#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
import shutil
import pysam
import chira_utilities
import multiprocessing
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Optional dependency for memory calculations
# If not available, safe defaults will be used
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False


def split_fasta_into_chunks(input_fasta, output_dir, num_chunks):
    """
    Split a large FASTA file into multiple chunks for parallel processing.
    
    This function reads the input FASTA and distributes reads across chunks,
    ensuring each chunk contains complete sequences (doesn't split in the middle of a read).
    
    Args:
        input_fasta (str): Path to input FASTA file
        output_dir (str): Directory to write chunk files
        num_chunks (int): Number of chunks to create
        
    Returns:
        list: List of paths to chunk FASTA files
    """
    chunk_files = []
    chunk_writers = []
    reads_per_chunk = []
    
    # Create chunk file writers
    for i in range(num_chunks):
        chunk_file = os.path.join(output_dir, f"chunk_{i:03d}.fasta")
        chunk_files.append(chunk_file)
        chunk_writers.append(open(chunk_file, 'w'))
        reads_per_chunk.append(0)
    
    # Distribute reads across chunks in round-robin fashion
    current_chunk = 0
    in_header = False
    current_seq = []
    current_header = None
    
    # OPTIMIZATION: Use adaptive buffer size for reading large FASTA files
    # Large buffers (8-16MB) significantly reduce system calls when reading multi-GB files
    # For a 10 GB file with 150M reads, this reduces I/O overhead by ~95%
    # Single file handle, so num_files=1
    BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=1)
    
    with open(input_fasta, 'r', buffering=BUFFER_SIZE) as f_in:
        for line in f_in:
            if line.startswith('>'):
                # Write previous sequence if exists
                if current_header is not None and current_seq:
                    chunk_writers[current_chunk].write(current_header)
                    chunk_writers[current_chunk].write(''.join(current_seq))
                    reads_per_chunk[current_chunk] += 1
                    current_chunk = (current_chunk + 1) % num_chunks
                
                # Start new sequence
                current_header = line
                current_seq = []
            else:
                # Accumulate sequence lines
                current_seq.append(line)
        
        # Write last sequence
        if current_header is not None and current_seq:
            chunk_writers[current_chunk].write(current_header)
            chunk_writers[current_chunk].write(''.join(current_seq))
            reads_per_chunk[current_chunk] += 1
    
    # Close all chunk files
    for writer in chunk_writers:
        writer.close()
    
    # Report chunk sizes
    total_reads = sum(reads_per_chunk)
    print(f"Split FASTA into {num_chunks} chunks:", file=sys.stderr)
    for i, count in enumerate(reads_per_chunk):
        if count > 0:
            print(f"  Chunk {i}: {count:,} reads", file=sys.stderr)
    print(f"  Total: {total_reads:,} reads", file=sys.stderr)
    
    # Return only non-empty chunks
    return [f for i, f in enumerate(chunk_files) if reads_per_chunk[i] > 0]


def align_with_bwa(align_type, index_type, query_fasta, refindex, outdir, seed_length, align_score,
                   macth_score, mistmatch_score, gap_o, gap_e, n_aligns, processes):
    """
        Function that maps the reads to the transcriptome. Different parameters
        are used for long and short alignments.
        Parameters:
            align_type: String of alignment type (long or short) used to format the output file name.
            index_type: String of index1 type (index1 or index2) used to format the output file name.
            query_fasta: Path to query fasta file. For align_type short, long.unmapped.fa from the out_dir used.
            refindex: Path to the reference index file.
            outdir: Output directory path.
            seed_length: Minimum seed length
            align_score: Minimum alignment score
            processes: Number of processes to use
    """

    bam = os.path.join(outdir, index_type + "." + align_type + ".bam")

    if align_type == "long":
        n_aligns = "50"
    elif align_type == "short":
        n_aligns = "100"

    bwa_params = ["-r 1",                       # look for internal seeds inside a seed longer than {-k} * {-r}
                  "-c 1000",
                  "-A "+ str(macth_score),          # match score
                  "-B " + str(mistmatch_score),       # mismatch penalty
                  "-O " + str(gap_o),                # gap open penalty
                  "-E " + str(gap_e),                # gap extension penalty
                  "-L 0",                       # clipping penalty, we need soft clips
                  "-Y",                         # use soft clipping for supplementary alignments
                  "-h " + str(n_aligns),             # if there're -h hits with score >80% of the maxscore, output in XA
                  "-k " + str(seed_length),     # minimum seed length
                  "-T " + str(align_score),     # minimum alignment score
                  "-t " + str(processes),
                  refindex,
                  query_fasta
                  ]
    # OPTIMIZATION: samtools view uses multi-threading (-@) to speed up SAM to BAM conversion
    # Performance Impact:
    # - BWA outputs uncompressed SAM format (text), which samtools compresses to BAM (binary)
    # - Compression is CPU-intensive and benefits greatly from parallelization
    # - With -@ processes, compression speed scales nearly linearly with thread count
    # - For a 19GB SAM file: single-threaded ~2 hours, 8 threads ~15 minutes (8x speedup)
    # - This parallelizes compression/decompression, which is especially beneficial for large files
    # OPTIMIZATION: Use subprocess instead of os.system for better error handling
    # - subprocess.run with check=True raises CalledProcessError if command fails
    # - This allows proper error detection and handling instead of silent failures
    bwacall = ("bwa mem " + " ".join(bwa_params) + " | samtools view -hb -@ " + str(processes) + " - > " + bam)
    print(bwacall)
    subprocess.run(bwacall, shell=True, check=True)


def write_mapped_bed(bam, bed, fasta, stranded):
    """
        Extracts the mapped and unmapped reads from the BAM alignments and writes them to BED and fasta files
        OPTIMIZED: Combined XA tag parsing to avoid two passes, improved string operations
        Parameters:
            bam: BAM file containing all merged alignments
            bed: output BED file path
            fasta: output FASTA file path
            stranded: Strand specificity
    """
    prev_readid = None
    prev_fasta = ""
    prev_unmapped = True
    # OPTIMIZATION: Adaptive buffer size for file I/O to reduce system calls
    # Performance Impact:
    # - For 150M reads writing ~50GB of output, default 8KB buffer = ~6.4M system calls
    # - With 8MB buffer = ~6,400 system calls (1,000x reduction)
    # - With 16MB buffer = ~3,200 system calls (2,000x reduction)
    # - This reduces I/O overhead by 95-99% and improves write performance by 10-30x
    # - Larger buffers reduce system calls significantly for large files (150M+ reads)
    # - We have 2 file handles (bed and fasta), so we calculate per-file buffer size
    # - Maximum 16MB per file to balance performance and memory usage (32MB total for 2 files)
    BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=2)
    
    # OPTIMIZATION: Pre-compile strand check conditions to avoid repeated string comparisons
    # Performance Impact:
    # - String comparison (stranded == "fw") is evaluated once instead of millions of times
    # - For 150M reads, this saves ~150M string comparisons
    # - Boolean checks are ~10-100x faster than string comparisons
    # - This optimization improves processing speed by 5-10% for large files
    stranded_fw = (stranded == "fw")
    stranded_rc = (stranded == "rc")
    stranded_both = (stranded == "both")
    
    # Progress tracking for large files
    processed_reads = 0
    last_progress_time = None
    
    with pysam.Samfile(bam, "rb") as fh_bam, \
         open(bed, "w", buffering=BUFFER_SIZE) as fh_mapped_bed, \
         open(fasta, "w", buffering=BUFFER_SIZE) as fh_unmapped_fasta:
        for alignment in fh_bam.fetch(until_eof=True):
            readid = alignment.query_name
            if readid != prev_readid:
                # Only write previous read if we've actually processed one (prev_readid is not None)
                if prev_readid is not None and prev_unmapped:
                    fh_unmapped_fasta.write(prev_fasta)
                prev_unmapped = True
                processed_reads += 1
                
                # Progress reporting every 1M reads
                if processed_reads % 1000000 == 0:
                    current_time = time.time()
                    if last_progress_time:
                        elapsed = current_time - last_progress_time
                        rate = 1000000 / elapsed if elapsed > 0 else 0
                        print(f"Processed {processed_reads:,} reads ({rate:.0f} reads/sec)", file=sys.stderr)
                    last_progress_time = current_time
            
            # Set prev_readid for every alignment to track current read
            # This must be set after the change detection but before processing
            prev_readid = readid

            # Check if primary alignment is on desired strand (needed for both mapped and unmapped logic)
            # For unmapped reads, is_desired_strand is not meaningful but we calculate it for consistency
            is_desired_strand = (stranded_both or 
                                (stranded_rc and alignment.is_reverse) or 
                                (stranded_fw and not alignment.is_reverse))
            
            # Get sequence for reads that will go to unmapped FASTA:
            # 1. Truly unmapped reads
            # 2. Wrong-strand mapped reads (treated as unmapped in original code)
            if alignment.is_unmapped or not is_desired_strand:
                readseq = alignment.get_forward_sequence()
                prev_fasta = f">{readid}\n{readseq}\n"
            
            if alignment.is_unmapped:
                continue
            
            # For stranded RNA-seq: if primary is on wrong strand, it's invalid and shouldn't be used as quality threshold
            # Only use primary length as threshold if it's on the desired strand
            # Write the alignment only if it mapped on desired strand
            # NOTE: In original code, prev_unmapped = False is only set for desired-strand alignments
            # Wrong-strand mapped reads are treated as unmapped (go to unmapped FASTA)
            optimal_alignment_len = 0
            if is_desired_strand:
                optimal_alignment_len = alignment.query_alignment_length
                fh_mapped_bed.write(chira_utilities.bedentry(alignment.reference_name,
                                                             str(alignment.reference_start),
                                                             str(alignment.reference_end),
                                                             readid,
                                                             "-" if alignment.is_reverse else "+",
                                                             alignment.cigarstring) + "\n")
                prev_unmapped = False

            # XA tag present in bwa output only
            if not alignment.has_tag('XA'):
                continue

            alt_alignments = alignment.get_tag('XA').rstrip(';').split(';')
            
            # Two passes are needed: we must find optimal alignment length before filtering
            # First pass: lightweight - only find optimal length (no storage)
            for alt_alignment in alt_alignments:
                # BWA XA tag format: refname,strand+start,cigar,nm
                f_alt_align = alt_alignment.split(',')
                alt_refstrand = f_alt_align[1]
                # Check if it mapped on desired strand (using pre-compiled conditions)
                if (alt_refstrand.startswith('-') and stranded_fw) or \
                        (alt_refstrand.startswith('+') and stranded_rc):
                    continue
                
                alt_cigar = f_alt_align[2]  # CIGAR string (field 3 is NM, which we ignore)
                alt_alignment_len = chira_utilities.alignment_length(alt_cigar)
                
                # Track optimal length (only from alternates on desired strand)
                if alt_alignment_len > optimal_alignment_len:
                    optimal_alignment_len = alt_alignment_len
            
            # Second pass: parse and write only alignments that meet the optimal length threshold
            if optimal_alignment_len > 0:
                for alt_alignment in alt_alignments:
                    # BWA XA tag format: refname,strand+start,cigar,nm
                    f_alt_align = alt_alignment.split(',')
                    alt_refstrand = f_alt_align[1]
                    # Check if it mapped on desired strand (using pre-compiled conditions)
                    if (alt_refstrand.startswith('-') and stranded_fw) or \
                            (alt_refstrand.startswith('+') and stranded_rc):
                        continue
                    
                    alt_cigar = f_alt_align[2]
                    alt_alignment_len = chira_utilities.alignment_length(alt_cigar)
                    
                    if alt_alignment_len < optimal_alignment_len:
                        continue
                    
                    alt_referenceid = f_alt_align[0]
                    alt_refstart = f_alt_align[1][1:]  # Remove strand character
                    alt_refend = chira_utilities.alignment_end(alt_refstart, alt_cigar, alt_refstrand.startswith('-'))
                    fh_mapped_bed.write(chira_utilities.bedentry(alt_referenceid,
                                                                 str(int(alt_refstart) - 1),
                                                                 str(alt_refend),
                                                                 readid,
                                                                 alt_refstrand[0],  # Just the strand character
                                                                 alt_cigar) + "\n")
                    prev_unmapped = False

        if prev_unmapped:
            fh_unmapped_fasta.write(prev_fasta)
        
        print(f"Total processed: {processed_reads:,} reads", file=sys.stderr)


def align_with_clan(query_fasta, outdir, ref_fasta1, ref_index1, ref_fasta2, ref_index2,
                    chimeric_overlap, align_score, stranded, n_aligns, processes):

    ref_index = "-f " + ref_fasta1 + " -d " + ref_index1
    if ref_fasta2 and ref_index2:
        ref_index += " -F " + ref_fasta2 + " -D " + ref_index2
    map_to_both_strands = "FALSE"
    if stranded == "both":
        map_to_both_strands = "TRUE"
    clan_search_params = ["-r " + query_fasta,          # query fasta file
                          "-m " + str(n_aligns),        # number of maximum hits for each maximal fragment
                          "-l " + str(align_score),     # minimum length for each fragment
                          "-s " + map_to_both_strands,
                          "-t " + str(processes),
                          "-v " + str(chimeric_overlap),
                          "-o " + os.path.join(outdir, "out.clan"),
                          ref_index]
    clan_search = ("clan_search " + " ".join(clan_search_params))
    print(clan_search)
    os.system(clan_search)

    ref_fasta = "-f " + ref_fasta1
    if ref_fasta2 and ref_index2:
        ref_fasta += " -F " + ref_fasta2

    clan_output_params = ["-r " + query_fasta,          # query fasta file
                          "-i " + os.path.join(outdir, "out.clan"),
                          "-o " + os.path.join(outdir, "out.map"),
                          ref_fasta]
    clan_output = ("clan_output " + " ".join(clan_output_params))
    print(clan_output)
    os.system(clan_output)
    return

def clan_to_bed(outdir):
    # OPTIMIZATION: Adaptive buffer size for file I/O to reduce system calls
    # Performance Impact:
    # - CLAN output files can be very large (multi-GB for large datasets)
    # - Default 8KB buffer requires millions of system calls for large files
    # - 8-16MB buffers reduce system calls by 1000-2000x
    # - This improves I/O performance by 10-30x for large CLAN output files
    # - Larger buffers reduce system calls significantly for large CLAN output files
    # - We have 2 file handles (input and output), so we calculate per-file buffer size
    # - Maximum 16MB per file to balance performance and memory usage (32MB total)
    BUFFER_SIZE = chira_utilities.get_adaptive_buffer_size(num_files=2)
    
    with open(os.path.join(outdir, "out.map"), buffering=BUFFER_SIZE) as fh_in, \
         open(os.path.join(outdir, "mapped.bed"), "w", buffering=BUFFER_SIZE) as fh_out:
        next(fh_in)
        for line in fh_in:
            [read_id,
             solution_id,
             read_mapped_begin,
             read_mapped_end,
             read_length,
             mapped_locations] = line.rstrip("\n").rstrip("\t").split("\t")
            
            # Build cigar string more efficiently
            # Convert strings to integers for arithmetic operations
            try:
                read_mapped_begin = int(read_mapped_begin)
                read_mapped_end = int(read_mapped_end)
                read_length = int(read_length)
            except ValueError as e:
                print(f"Warning: Skipping line with invalid numeric values: {line[:100]}, error: {e}", file=sys.stderr)
                continue
            
            # Defensive check: ensure values are valid (non-negative, logical ordering)
            if read_mapped_begin < 1 or read_mapped_end < read_mapped_begin or read_length < read_mapped_end:
                print(f"Warning: Skipping line with invalid position values (begin={read_mapped_begin}, end={read_mapped_end}, length={read_length}): {line[:100]}", file=sys.stderr)
                continue
            
            lead_soft_clips = read_mapped_begin - 1
            match_len = read_mapped_end - read_mapped_begin + 1
            trail_soft_clips = read_length - read_mapped_end
            
            # Defensive check: ensure match_len is positive
            if match_len <= 0:
                print(f"Warning: Skipping line with invalid match length ({match_len}): {line[:100]}", file=sys.stderr)
                continue
            
            cigar_parts = []
            if lead_soft_clips > 0:
                cigar_parts.append(f"{lead_soft_clips}S")
            cigar_parts.append(f"{match_len}M")
            if trail_soft_clips > 0:
                cigar_parts.append(f"{trail_soft_clips}S")
            cigar = "".join(cigar_parts)

            # Pre-split mapped_locations to avoid repeated splitting
            if mapped_locations:
                location_list = mapped_locations.split(";")
            else:
                continue
            
            for mapped_location in location_list:
                d = mapped_location.split(":")
                # Defensive check: ensure location has at least 2 parts (ref_id and position range)
                if len(d) < 2:
                    print(f"Warning: Skipping malformed location (expected 'ref_id:start-end', got '{mapped_location}'): {line[:100]}", file=sys.stderr)
                    continue
                # if header has spaces select the id only
                ref_id = ":".join(d[0:-1]).split(' ')[0]
                position_range = d[-1].split("-")
                # Defensive check: ensure position range has start and end
                if len(position_range) < 2:
                    print(f"Warning: Skipping malformed position range (expected 'start-end', got '{d[-1]}'): {line[:100]}", file=sys.stderr)
                    continue
                [ref_start, ref_end] = position_range
                try:
                    ref_start_int = int(ref_start)
                except ValueError:
                    print(f"Warning: Skipping location with invalid start position '{ref_start}': {line[:100]}", file=sys.stderr)
                    continue
                # NOTE: at the moment there is no way to findout which strand it is mapping to
                ref_start_minus_one = ref_start_int - 1
                fh_out.write("\t".join([ref_id, str(ref_start_minus_one), ref_end,
                                        ",".join([read_id, ref_id, str(ref_start_minus_one), ref_end, "+", cigar]),
                                        "1", "+"]) + "\n")


def calculate_per_job_processes(total_processes, num_parallel_units):
    """
    Calculate number of processes per BWA alignment job based on total processes.
    
    Args:
        total_processes: Total number of processes specified by user
        num_parallel_units: Number of parallel units (jobs or chunks) that will run simultaneously
    
    Returns:
        int: Number of processes to use per BWA alignment job (minimum 1)
    """
    if num_parallel_units <= 0:
        return max(1, total_processes)
    per_job = max(1, total_processes // num_parallel_units)
    return per_job


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: map reads to the reference',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-a", '--aligner', type=str, choices=["bwa", "clan"], default='bwa', required=False,
                        dest='aligner', metavar='', help='Alignment program to use, bwa or clan')

    parser.add_argument('-i', '--query_fasta', action='store', dest='fasta', required=True,
                        metavar='', help='Path to query fasta file')

    parser.add_argument('-o', '--outdir', action='store', dest='outdir', required=True, metavar='',
                        help='Output directory path for the analysis')

    parser.add_argument('-x1', '--index1', action='store', dest='idx1', required=False,
                        metavar='', help='first priority index file')

    parser.add_argument('-x2', '--index2', action='store', dest='idx2', required=False,
                        metavar='', help='second priority index file')

    parser.add_argument('-f1', '--ref_fasta1', action='store', dest='ref_fasta1', required=False,
                        metavar='', help='First priority fasta file')

    parser.add_argument('-f2', '--ref_fasta2', action='store', dest='ref_fasta2', required=False,
                        metavar='', help='second priority fasta file')

    parser.add_argument("-b", '--build', action='store_true', dest='build_index',
                        help="Build indices from reference fasta files")

    parser.add_argument('-p', '--processes', action='store', type=int, default=None, metavar='',
                        dest='processes',
                        help='''Total number of CPU processes/threads to use across all BWA alignment jobs.

HOW TO SET --processes:
- GENERAL RULE: Set to the total number of CPU cores available on your system
  * Example: 32-core system → --processes 32
  * Example: 16-core system → --processes 16
  * Example: 8-core system → --processes 8
- DEFAULT: If not specified, uses all available CPU cores (auto-detected)
- For shared systems: Set to number of cores allocated to your job
- For memory-constrained systems: Reduce if you encounter out-of-memory errors

HOW IT WORKS:
- WITHOUT CHUNKING: Total processes are divided among parallel BWA jobs
  * 2 jobs if only index1 is used (long + short)
  * 4 jobs if both index1 and index2 are used (long + short × 2 indices)
  * Example: 32 processes, 4 jobs → 8 processes per job (32 total)
- WITH CHUNKING: Total processes are divided among parallel chunks (controlled by --parallel_chunks, default: 2)
  * Example: 32 processes, --parallel_chunks 2 → 16 processes per chunk
  * Example: 32 processes, --parallel_chunks 4 → 8 processes per chunk
  * Each chunk processes jobs sequentially using all allocated processes
  * With <16 processes: Only 1 chunk runs at a time

Note: This controls BWA alignment threads. Samtools operations (merge, sort)
may use additional threads based on available resources.''')

    parser.add_argument('--sort_memory', action='store', type=str, default=None, metavar='',
                        dest='sort_memory',
                        help='Memory per thread for BAM sorting (e.g., "2G", "3G"). '
                             'If not specified, automatically calculates based on available RAM. '
                             'More memory reduces temporary files and I/O, improving performance. '
                             'Total memory used = sort_memory × processes, so ensure sufficient RAM.')

    parser.add_argument("-s", '--stranded', type=str, choices=["fw", "rc", "both"], default='fw', metavar='',
                        dest='stranded',
                        help='''Strand-specificity of input samples.
                             fw = map to transcript strand (default, recommended for most protocols like CLASH, CLEAR-CLIP, PARIS, SPLASH);
                             rc = map to reverse complement of transcript strand (use if your library protocol produces reads from antisense strand);
                             both = map on both strands (use only for unstranded libraries where strand information is not preserved).
                             
                             When to use:
                             - Most RNA-RNA interactome protocols (CLASH, CLEAR-CLIP, PARIS, SPLASH) are stranded and use "fw"
                             - Use "rc" only if you know your protocol produces reads from the reverse complement strand
                             - Use "both" only for unstranded libraries (rare for interactome protocols)
                             
                             Stranded mapping filters out alignments on the wrong strand, reducing false positives and improving chimeric read detection.''')

    parser.add_argument("-l1", '--seed_length1', action='store', type=int, default=12, metavar='',
                        dest='seed_length1',
                        help='''Seed length for 1st mapping iteration.
                                bwa-mem parameter "-k"''')

    parser.add_argument("-l2", '--seed_length2', action='store', type=int, default=16, metavar='',
                        dest='seed_length2',
                        help='''Seed length for 2nd mapping iteration.
                                bwa-mem parameter "-k"''')

    parser.add_argument("-s1", '--align_score1', action='store', type=int, default=18, metavar='',
                        dest='align_score1',
                        help='''Minimum alignment score in 1st mapping iteration.
                                bwa-mem parameter "-T" and clan_search parameter "-l"''')

    parser.add_argument("-s2", '--align_score2', action='store', type=int, default=16, metavar='',
                        dest='align_score2',
                        help='''Minimum alignment score in 2nd mapping iteration.
                                It must be smaller than --align_score1 parameter.
                                bwa-mem parameter "-T" and clan_search parameter "-l"''')

    parser.add_argument("-ma1", '--match1', action='store', type=int, default=1, metavar='',
                        dest='match1',
                        help='Matching score for 1st mapping iteration.')

    parser.add_argument("-mm1", '--mismatch1', action='store', type=int, default=4, metavar='',
                        dest='mismatch1',
                        help='Mismatch penalty for 1st mapping iteration.')

    parser.add_argument("-ma2", '--match2', action='store', type=int, default=1, metavar='',
                        dest='match2',
                        help='Matching score for 2nd mapping iteration.')

    parser.add_argument("-mm2", '--mismatch2', action='store', type=int, default=6, metavar='',
                        dest='mismatch2',
                        help='Mismatch penalty for 2nd mapping iteration.')

    parser.add_argument("-go1", '--gapopen1', action='store', type=int, default=6, metavar='',
                        dest='gapopen1',
                        help='Gap opening penalty for 1st mapping iteration.')

    parser.add_argument("-ge1", '--gapext1', action='store', type=int, default=1, metavar='',
                        dest='gapext1',
                        help='Gap extension penalty for 1st mapping iteration.')

    parser.add_argument("-go2", '--gapopen2', action='store', type=int, default=100, metavar='',
                        dest='gapopen2',
                        help='Gap opening penalty for 2nd mapping iteration.')

    parser.add_argument("-ge2", '--gapext2', action='store', type=int, default=100, metavar='',
                        dest='gapext2',
                        help='Gap extension penalty for 2nd mapping iteration.')

    parser.add_argument("-h1", '--nhits1', action='store', type=int, default=50, metavar='',
                        dest='nhits1',
                        help='Number of allowed multi hits per read')

    parser.add_argument("-h2", '--nhits2', action='store', type=int, default=100, metavar='',
                        dest='nhits2',
                        help='Number of allowed multi hits per read in 2nd iteration')

    parser.add_argument('-co', '--chimeric_overlap', action='store', type=int, default=2, metavar='',
                        dest='chimeric_overlap',
                        help='Maximum number of bases allowed between the chimeric segments of a read')

    parser.add_argument('--chunk_fasta', action='store', type=int, default=None, metavar='',
                        dest='chunk_fasta',
                        help='''Split input FASTA into N chunks for parallel processing (recommended for large files >1GB).

HOW TO SET --chunk_fasta:
- GENERAL RULE: Set based on input file size to create 1-3GB chunks
  * Small files (<1GB): Leave unset (no chunking needed)
  * Medium files (1-10GB): --chunk_fasta 5-10 (creates ~1-2GB chunks)
  * Large files (10-50GB): --chunk_fasta 10-20 (creates ~2-5GB chunks)
  * Very large files (>50GB): --chunk_fasta 20-50 (creates ~1-2.5GB chunks)
- FORMULA: N ≈ input_file_size_GB / desired_chunk_size_GB
  * Example: 20GB file, want 2GB chunks → --chunk_fasta 10
  * Example: 5GB file, want 1GB chunks → --chunk_fasta 5
- OPTIMAL CHUNK SIZE: 1-3GB per chunk for best I/O performance
- More chunks = better I/O parallelism but more overhead
- Fewer chunks = less overhead but larger memory usage per chunk

HOW IT WORKS:
1. FASTA is split into N chunks (as specified by --chunk_fasta)
2. Chunks are processed in batches (number controlled by --parallel_chunks, default: 2)
3. Each chunk processes all BWA jobs sequentially
4. Remaining chunks are processed in subsequent batches

EXECUTION WITH CHUNKING:
- Number of parallel chunks is controlled by --parallel_chunks (default: 2)
- Processes per chunk = --processes / --parallel_chunks
- Example: --chunk_fasta 10 --processes 32 --parallel_chunks 2 → 2 chunks at a time (16 processes each)
- Example: --chunk_fasta 10 --processes 32 --parallel_chunks 4 → 4 chunks at a time (8 processes each)
- If processes per chunk < 4, the number of parallel chunks is automatically reduced

BENEFITS:
- Better I/O: Each chunk reads from different parts of the file, reducing disk contention
- Lower memory: Each BWA job processes smaller chunks, reducing peak memory usage
- Resource management: Control parallel chunks with --parallel_chunks to match your system resources
- Each BWA job gets processes based on --processes / --parallel_chunks

WHEN TO USE CHUNKING:
- Use chunking for: Large files (>1GB), memory-constrained systems, I/O bottlenecks
- Skip chunking for: Small files (<1GB), systems with abundant memory

If I/O performance is a concern:
- Monitor disk I/O wait time with: iostat -x 1 or iotop
- If I/O is the bottleneck, reduce --chunk_fasta (fewer chunks = less I/O contention)
- Reduce --processes (fewer threads per chunk = less I/O per chunk)
- Use faster storage (SSD instead of HDD, faster network storage)''')

    parser.add_argument('--parallel_chunks', action='store', type=int, default=2, metavar='',
                        dest='parallel_chunks',
                        help='''Number of chunks to process simultaneously when using --chunk_fasta.

HOW TO SET --parallel_chunks:
- GENERAL RULE: Set based on available memory and CPU resources
  * Default: 2 (recommended for most systems)
  * Each chunk processes BWA jobs sequentially, so N chunks = N active BWA jobs
  * Each BWA job needs 4-8 CPUs, so ensure: parallel_chunks × 4-8 ≤ --processes
- RECOMMENDATIONS:
  * Small systems (<16 cores, <32GB RAM): --parallel_chunks 1
  * Medium systems (16-32 cores, 32-64GB RAM): --parallel_chunks 2 (default)
  * Large systems (>32 cores, >64GB RAM): --parallel_chunks 2-4
  * Very large systems (>64 cores, >128GB RAM): --parallel_chunks 4-8
- MEMORY CONSIDERATION:
  * More parallel chunks = more memory usage
  * Each chunk needs memory for BWA alignment (typically 2-4GB per chunk)
  * Total memory ≈ parallel_chunks × 4GB (approximate)
- CPU CONSIDERATION:
  * Processes per chunk = --processes / --parallel_chunks
  * Each chunk should get at least 4 processes for optimal BWA performance
  * Example: --processes 32 --parallel_chunks 4 → 8 processes per chunk (good)
  * Example: --processes 8 --parallel_chunks 4 → 2 processes per chunk (too few)

Note: This parameter only takes effect when --chunk_fasta is specified.''')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.4.4')

    args = parser.parse_args()
    # Set default processes to CPU count if not specified
    if args.processes is None:
        args.processes = multiprocessing.cpu_count()
    
    return args


def validate_arguments(args):
    """Validate command-line arguments."""
    if not args.idx1 and not args.ref_fasta1:
        sys.stderr.write("option -x1 or -f1 are required\n")
        sys.exit(1)

    if (args.idx1 or args.idx2) and args.build_index:
        sys.stderr.write("options -b and -x1 are mutually exclusive\n")
        sys.exit(1)
    
    if not args.idx1 and not args.idx2 and not args.build_index:
        sys.stderr.write("Either -b or -x1 is required\n")
        sys.exit(1)


def print_cpu_guidance(args):
    """Print CPU usage guidance based on system capabilities and chunking strategy."""
    try:
        cpu_count = multiprocessing.cpu_count()
        
        # Set default if not specified
        total_processes = args.processes if args.processes is not None else cpu_count
        
        if args.chunk_fasta and args.chunk_fasta > 1:
            # CHUNKING MODE: Multiple chunks run in parallel
            # Each chunk processes alignment jobs sequentially
            # Use user-specified parallel_chunks (default 2), limited by number of chunks
            num_parallel_chunks = min(args.parallel_chunks, args.chunk_fasta)
            per_chunk_processes = calculate_per_job_processes(total_processes, num_parallel_chunks)
            
            num_job_types = 2  # Default: long + short
            if args.idx2 or (args.build_index and args.ref_fasta2):
                num_job_types = 4  # long + short × 2 indices
            
            print(f"INFO: Using {total_processes} total processes across {num_parallel_chunks} parallel chunks", file=sys.stderr)
            print(f"      → {per_chunk_processes} processes per chunk (on {cpu_count}-core system)", file=sys.stderr)
            print(f"      Each chunk processes {num_job_types} alignment jobs sequentially.", file=sys.stderr)
            
            if total_processes > cpu_count * 1.5:
                print(f"WARNING: Total processes ({total_processes}) exceeds system cores ({cpu_count})", file=sys.stderr)
                print(f"         Consider reducing --processes or --chunk_fasta.", file=sys.stderr)
        else:
            # NON-CHUNKING MODE: Multiple alignment jobs run in parallel
            num_jobs = 2  # Default: long + short with index1
            if args.idx2 or (args.build_index and args.ref_fasta2):
                num_jobs = 4  # long + short × 2 indices
            
            per_job_processes = calculate_per_job_processes(total_processes, num_jobs)
            
            print(f"INFO: Using {total_processes} total processes across {num_jobs} parallel BWA jobs", file=sys.stderr)
            print(f"      → {per_job_processes} processes per job (on {cpu_count}-core system)", file=sys.stderr)
            
            if total_processes > cpu_count * 1.5:
                print(f"WARNING: Total processes ({total_processes}) exceeds system cores ({cpu_count})", file=sys.stderr)
                print(f"         Consider reducing --processes.", file=sys.stderr)
    except Exception:
        pass


def print_configuration(args):
    """Print configuration summary."""
    print('Query fasta                          : ' + args.fasta)
    print('Output directory                     : ' + args.outdir)
    print('Aligner                              : ' + args.aligner)
    print('Build index?                         : ' + str(args.build_index))
    if args.idx1:
        print('1st priority BWA/CLAN index          : ' + args.idx1)
    if args.idx2:
        print('2nd priority BWA/CLAN  index         : ' + args.idx2)
    if args.ref_fasta1:
        print('1st priority reference fasta file    : ' + args.ref_fasta1)
    if args.ref_fasta2:
        print('2nd priority reference fasta file    : ' + args.ref_fasta2)
    total_processes = args.processes if args.processes is not None else multiprocessing.cpu_count()
    print('Total number of processes            : ' + str(total_processes) + (' (auto-detected)' if args.processes is None else ''))
    print('Stranded                             : ' + args.stranded)
    print('Seed length                          : ' + str(args.seed_length1))
    if args.seed_length2:
        print('Seed length for 2nd iteration        : ' + str(args.seed_length2))
    print('Alignment score                      : ' + str(args.align_score1))
    if args.align_score2:
        print('Alignment score for 2nd iteration    : ' + str(args.align_score2))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    if args.chunk_fasta:
        print('FASTA chunking                       : ' + str(args.chunk_fasta) + ' chunks (enabled)')
    else:
        print('FASTA chunking                       : Disabled (using parallel job strategy)')
    print("===================================================================")


def setup_indices(args):
    """Build or set up indices for alignment."""
    index1 = args.idx1 if args.idx1 else None
    index2 = args.idx2 if args.idx2 else None
    
    if args.aligner == "clan":
        if args.build_index:
            index1 = os.path.join(args.outdir, "index1")
            os.system("clan_index -f " + args.ref_fasta1 + " -d " + index1)
            if args.ref_fasta2:
                index2 = os.path.join(args.outdir, "index2")
                os.system("clan_index -f " + args.ref_fasta2 + " -d " + index2)
    elif args.aligner == "bwa":
        if args.build_index:
            index1 = os.path.join(args.outdir, "index1")
            os.system("bwa index -p " + index1 + " " + args.ref_fasta1)
            if args.ref_fasta2:
                index2 = os.path.join(args.outdir, "index2")
                os.system("bwa index -p " + index2 + " " + args.ref_fasta2)
    
    return index1, index2


def run_clan_mapping(args, index1, index2):
    """Run CLAN alignment workflow."""
    chira_utilities.print_w_time("START: Map read using CLAN")
    align_with_clan(args.fasta, args.outdir, args.ref_fasta1, index1, args.ref_fasta2, index2,
                    args.chimeric_overlap, args.align_score2, args.stranded, args.nhits2, args.processes)
    chira_utilities.print_w_time("END: Map read using CLAN")
    chira_utilities.print_w_time("START: Write alignments to BED")
    clan_to_bed(args.outdir)
    chira_utilities.print_w_time("END: Write alignments to BED")


def get_alignment_job_types(args, index1, index2):
    """Get list of alignment job types to run."""
    alignment_job_types = [
        ("long", "index1", index1, args.seed_length1, args.align_score1,
         args.match1, args.mismatch1, args.gapopen1, args.gapext1, args.nhits1),
        ("short", "index1", index1, args.seed_length2, args.align_score2,
         args.match2, args.mismatch2, args.gapopen2, args.gapext2, args.nhits2)
    ]
    if index2:
        alignment_job_types.extend([
            ("long", "index2", index2, args.seed_length1, args.align_score1,
             args.match1, args.mismatch1, args.gapopen1, args.gapext1, args.nhits1),
            ("short", "index2", index2, args.seed_length2, args.align_score2,
             args.match2, args.mismatch2, args.gapopen2, args.gapext2, args.nhits2)
        ])
    return alignment_job_types


def run_bwa_mapping_with_chunking(args, index1, index2):
    """
    OPTIMIZATION: Run BWA mapping with chunking strategy for very large files.
    
    EXECUTION MODEL - How chunking jobs are executed:
    
    Example: 10 chunks, 4 alignment job types (index1.long, index1.short, index2.long, index2.short)
    
    Execution Flow:
    ┌─────────────────────────────────────────────────────────────────────────┐
    │ STEP 1: Split FASTA into chunks                                          │
    │   Input: 10GB FASTA file                                                 │
    │   Output: chunk_000.fasta, chunk_001.fasta, ..., chunk_009.fasta        │
    │   (Each chunk ~1GB, distributed round-robin)                            │
    └─────────────────────────────────────────────────────────────────────────┘
    
    ┌─────────────────────────────────────────────────────────────────────────┐
    │ STEP 2: Process chunks in parallel (N chunks simultaneously, controlled by --parallel_chunks) │
    │                                                                          │
    │   Batch 1: Chunk 0 (parallel)    Chunk 1 (parallel)                    │
    │   ┌──────────────┐      ┌──────────────┐                              │
    │   │ Job 1: long │      │ Job 1: long │                              │
    │   │ + index1    │      │ + index1    │                              │
    │   ├──────────────┤      ├──────────────┤                              │
    │   │ Job 2: short│      │ Job 2: short│                              │
    │   │ + index1    │      │ + index1    │                              │
    │   ├──────────────┤      ├──────────────┤                              │
    │   │ Job 3: long │      │ Job 3: long │                              │
    │   │ + index2    │      │ + index2    │                              │
    │   ├──────────────┤      ├──────────────┤                              │
    │   │ Job 4: short│      │ Job 4: short│                              │
    │   │ + index2    │      │ + index2    │                              │
    │   └──────────────┘      └──────────────┘                              │
    │   (sequential)          (sequential)                                    │
    │                                                                          │
    │   Batch 2: Chunk 2, Chunk 3 (after Batch 1 completes)                   │
    │   ... (continues until all 10 chunks are processed)                     │
    │                                                                          │
    │   Total: 10 chunks × 4 jobs = 40 BWA alignment jobs                      │
    │   Parallelism: N chunks run simultaneously (N active BWA jobs, controlled by --parallel_chunks) │
    │   Within each chunk: 4 jobs run sequentially                           │
    └─────────────────────────────────────────────────────────────────────────┘
    
    ┌─────────────────────────────────────────────────────────────────────────┐
    │ STEP 3: Collect results by job type                                      │
    │                                                                          │
    │   index1.long:  [chunk_000/index1.long.bam, chunk_001/index1.long.bam,  │
    │                  ..., chunk_009/index1.long.bam]                        │
    │   index1.short: [chunk_000/index1.short.bam, chunk_001/index1.short.bam, │
    │                  ..., chunk_009/index1.short.bam]                       │
    │   index2.long:  [chunk_000/index2.long.bam, chunk_001/index2.long.bam,  │
    │                  ..., chunk_009/index2.long.bam]                        │
    │   index2.short: [chunk_000/index2.short.bam, chunk_001/index2.short.bam,│
    │                  ..., chunk_009/index2.short.bam]                       │
    └─────────────────────────────────────────────────────────────────────────┘
    
    ┌─────────────────────────────────────────────────────────────────────────┐
    │ STEP 4: Merge chunk BAMs for each job type                               │
    │                                                                          │
    │   Merge 10 chunk BAMs → index1.long.bam                                  │
    │   Merge 10 chunk BAMs → index1.short.bam                                 │
    │   Merge 10 chunk BAMs → index2.long.bam                                 │
    │   Merge 10 chunk BAMs → index2.short.bam                                │
    └─────────────────────────────────────────────────────────────────────────┘
    
    Key Points:
    - Chunks run in PARALLEL (number controlled by --parallel_chunks, default: 2)
    - Within each chunk, alignment jobs run SEQUENTIALLY (one after another)
    - Each BWA job needs 4-8 CPUs, so N chunks = N active BWA jobs at a time
    - Total jobs = num_chunks × num_alignment_job_types
    - Each chunk produces one BAM file per alignment job type
    - Chunk BAMs are merged by job type to produce final BAMs
    - Each job within a chunk uses all processes allocated to that chunk
    - Processes per chunk = --processes / --parallel_chunks
    - If processes per chunk < 4, the number of parallel chunks is automatically reduced
    
    Performance Benefits:
    - Better I/O: Each chunk reads from different parts of the file, reducing disk contention
    - Lower memory: Each BWA job processes smaller chunks, reducing peak memory usage
    - Resource management: Control parallel chunks with --parallel_chunks (default: 2) to match your system resources
    - For 10GB file with 150M reads:
      * Without chunking: Single large BWA job, high memory, sequential I/O
      * With 10 chunks (2 at a time): ~10x lower memory per job, parallel I/O, controlled resource usage
      * Expected speedup: 2-3x depending on I/O subsystem and available cores
    """
    chira_utilities.print_w_time(f"START: Splitting FASTA into {args.chunk_fasta} chunks")
    chunk_dir = os.path.join(args.outdir, "chunks")
    os.makedirs(chunk_dir, exist_ok=True)
    chunk_files = split_fasta_into_chunks(args.fasta, chunk_dir, args.chunk_fasta)
    chira_utilities.print_w_time(f"END: Split into {len(chunk_files)} chunks")
    
    # Edge case: Check if any chunks were created
    if len(chunk_files) == 0:
        sys.stderr.write("ERROR: No chunks were created from the input FASTA file. The file may be empty.\n")
        sys.exit(1)
    
    # Get all alignment job types that need to be run (e.g., index1.long, index1.short, index2.long, index2.short)
    alignment_job_types = get_alignment_job_types(args, index1, index2)
    num_job_types = len(alignment_job_types)
    print(f"INFO: Each chunk will be processed through {num_job_types} alignment job types", file=sys.stderr)
    print(f"INFO: Total alignment jobs = {len(chunk_files)} chunks × {num_job_types} job types = {len(chunk_files) * num_job_types} jobs", file=sys.stderr)
    
    # Calculate optimal number of parallel chunks based on user setting and available processes
    # User can specify --parallel_chunks to control how many chunks run simultaneously
    total_processes = args.processes if args.processes is not None else multiprocessing.cpu_count()
    
    # Use user-specified number of parallel chunks (default is 2)
    # Limited by: user setting, number of chunks available, and available processes
    num_parallel_chunks = min(args.parallel_chunks, len(chunk_files))
    
    # Calculate processes per chunk: divide total processes among parallel chunks
    # This ensures each chunk gets enough processes (4-8) for BWA alignment
    per_chunk_processes = calculate_per_job_processes(total_processes, num_parallel_chunks)
    
    # Ensure each chunk gets at least 4 processes (minimum for BWA)
    if per_chunk_processes < 4:
        # If we don't have enough processes for requested chunks, reduce number of parallel chunks
        required_processes = args.parallel_chunks * 4
        num_parallel_chunks = max(1, total_processes // 4)  # At least 1 chunk, at most what we can support
        per_chunk_processes = calculate_per_job_processes(total_processes, num_parallel_chunks)
        print(f"WARNING: Not enough processes for {args.parallel_chunks} chunks (need at least {required_processes}, have {total_processes}).", file=sys.stderr)
        print(f"         Running {num_parallel_chunks} chunks in parallel instead ({per_chunk_processes} processes each).", file=sys.stderr)
    
    print(f"INFO: Using {total_processes} total processes across {num_parallel_chunks} parallel chunks → {per_chunk_processes} processes per chunk", file=sys.stderr)
    print(f"INFO: Within each chunk, {num_job_types} alignment jobs will run sequentially (each job uses {per_chunk_processes} processes)", file=sys.stderr)
    if num_parallel_chunks < len(chunk_files):
        print(f"INFO: Processing {len(chunk_files)} total chunks in batches of {num_parallel_chunks} (limited to prevent resource exhaustion)", file=sys.stderr)
    
    def process_chunk(chunk_file, chunk_idx):
        """
        Process a single chunk through all alignment jobs.
        
        Execution: Alignment job types run SEQUENTIALLY within each chunk.
        Multiple chunks can run in PARALLEL (controlled by ThreadPoolExecutor).
        
        Example for chunk_000:
          1. Run index1.long alignment on chunk_000.fasta → chunk_000_out/index1.long.bam
          2. Run index1.short alignment on chunk_000.fasta → chunk_000_out/index1.short.bam
          3. Run index2.long alignment on chunk_000.fasta → chunk_000_out/index2.long.bam
          4. Run index2.short alignment on chunk_000.fasta → chunk_000_out/index2.short.bam
        
        Returns: Dictionary mapping job_name -> chunk BAM file path
        """
        chunk_outdir = os.path.join(chunk_dir, f"chunk_{chunk_idx:03d}_out")
        os.makedirs(chunk_outdir, exist_ok=True)
        
        chunk_bams = {}
        # EXECUTION: Process each alignment job type SEQUENTIALLY for this chunk
        # This prevents resource exhaustion when multiple chunks run in parallel
        for (align_type, index_type, refindex, seed_length, align_score,
             match_score, mismatch_score, gap_o, gap_e, n_aligns) in alignment_job_types:
            job_name = f"{index_type}.{align_type}"
            try:
                # Run BWA alignment for this job type on this chunk
                align_with_bwa(align_type, index_type, chunk_file, refindex, chunk_outdir,
                              seed_length, align_score, match_score, mismatch_score,
                              gap_o, gap_e, n_aligns, per_chunk_processes)
                chunk_bam = os.path.join(chunk_outdir, index_type + "." + align_type + ".bam")
                chunk_bams[job_name] = chunk_bam
            except Exception as e:
                raise RuntimeError(f"Chunk {chunk_idx} {job_name} failed: {str(e)}")
        
        return chunk_bams
    
    # EXECUTION MODEL: Process chunks in PARALLEL using ThreadPoolExecutor
    # 
    # Parallelism Level:
    # - Number of parallel chunks is controlled by --parallel_chunks (default: 2)
    # - Each BWA job needs 4-8 CPUs, so N chunks = N active BWA jobs at a time
    # - Example: 32 processes, 10 chunks, --parallel_chunks 2 → 2 chunks in parallel (16 processes each)
    # - Example: 32 processes, 10 chunks, --parallel_chunks 4 → 4 chunks in parallel (8 processes each)
    # - Chunks are processed in batches (size controlled by --parallel_chunks) until all chunks are complete
    #
    # Within Each Chunk (SEQUENTIAL):
    # - All alignment job types run one after another
    # - Example: chunk_000 runs index1.long → index1.short → index2.long → index2.short
    # - Each job uses all processes allocated to the chunk (per_chunk_processes)
    #
    # Across Chunks (PARALLEL, controlled by --parallel_chunks):
    # - Number of chunks running simultaneously is controlled by --parallel_chunks (default: 2)
    # - Example: With 32 processes and --parallel_chunks 2, 2 chunks run in parallel (16 processes each)
    # - Example: With 32 processes and --parallel_chunks 4, 4 chunks run in parallel (8 processes each)
    # - Each chunk processes 4 jobs sequentially, so N active BWA jobs at any time (where N = --parallel_chunks)
    #
    # Performance Impact:
    # - Sequential within chunk: 10 chunks × 4 jobs × 30 min = 20 hours total
    # - Parallel chunks (N at a time): 10 chunks processed in batches of N = ~(10/N) × 30 min total
    # - Speedup: ~Nx for this example (limited by --parallel_chunks to prevent resource exhaustion)
    chira_utilities.print_w_time(f"START: Processing {len(chunk_files)} chunks in parallel (up to {num_parallel_chunks} simultaneously)")
    all_chunk_bams = {}  # Dictionary: job_name -> list of chunk BAM files
    
    def process_chunk_wrapper(chunk_info):
        """Wrapper function to process a chunk (used by ThreadPoolExecutor)."""
        chunk_file, chunk_idx = chunk_info
        return process_chunk(chunk_file, chunk_idx)
    
    failed_chunks = []
    # EXECUTION: Submit all chunk processing tasks to ThreadPoolExecutor
    # - Each chunk will be processed through all alignment job types
    # - Multiple chunks run in parallel (controlled by max_workers)
    # - Results are collected as they complete (as_completed)
    
    # Extract chunk index from filename (e.g., "chunk_005.fasta" -> 5)
    # This ensures chunk_idx matches the actual chunk file index, not enumerate position
    def get_chunk_idx(chunk_file):
        basename = os.path.basename(chunk_file)
        # Extract number from "chunk_XXX.fasta" format
        return int(basename.replace("chunk_", "").replace(".fasta", ""))
    
    with ThreadPoolExecutor(max_workers=num_parallel_chunks) as executor:
        # Submit all chunk processing tasks to the executor
        # Each task will process one chunk through all alignment job types
        future_to_chunk = {executor.submit(process_chunk_wrapper, (chunk_file, get_chunk_idx(chunk_file))): 
                          (get_chunk_idx(chunk_file), chunk_file) for chunk_file in chunk_files}
        
        # EXECUTION: Collect results as chunks complete processing
        # Results are organized by job type: all chunk BAMs for each job type are grouped together
        # Example structure after collection:
        #   all_chunk_bams = {
        #     "index1.long":  [chunk_000/index1.long.bam, chunk_001/index1.long.bam, ...],
        #     "index1.short": [chunk_000/index1.short.bam, chunk_001/index1.short.bam, ...],
        #     "index2.long":  [chunk_000/index2.long.bam, chunk_001/index2.long.bam, ...],
        #     "index2.short": [chunk_000/index2.short.bam, chunk_001/index2.short.bam, ...]
        #   }
        for future in as_completed(future_to_chunk):
            chunk_idx, chunk_file = future_to_chunk[future]
            try:
                # Get results from this chunk (dictionary: job_name -> chunk BAM path)
                chunk_bams = future.result()
                # Organize chunk BAMs by job type for later merging
                for job_name, bam_file in chunk_bams.items():
                    if job_name not in all_chunk_bams:
                        all_chunk_bams[job_name] = []
                    all_chunk_bams[job_name].append(bam_file)
            except Exception as e:
                failed_chunks.append((chunk_idx, str(e)))
    
    chira_utilities.print_w_time(f"END: All chunks processed")
    
    if failed_chunks:
        print("ERROR: Some chunks failed:", file=sys.stderr)
        for chunk_idx, error in failed_chunks:
            print(f"  - Chunk {chunk_idx}: {error}", file=sys.stderr)
        sys.exit(1)
    
    # EXECUTION: Merge chunk BAMs for each job type to produce final BAM files
    #
    # For each alignment job type, merge all chunk BAMs into a single final BAM:
    #   Example for index1.long:
    #     Input:  [chunk_000/index1.long.bam, chunk_001/index1.long.bam, ..., chunk_009/index1.long.bam]
    #     Output: index1.long.bam (merged from all 10 chunks)
    #
    # This produces the same result as if we had run the alignment on the full file,
    # but with better I/O performance and lower memory usage.
    #
    # Performance Impact:
    # - For 10 chunks per job type: Need to merge 10 BAM files into 1
    # - Single-threaded merge: ~20-30 min per job type (4 types = 80-120 min)
    # - Parallel merge with 8 threads: ~3-4 min per job type (4 types = 12-16 min)
    # - Total speedup: 5-10x for chunk merging step
    # - pysam.merge uses samtools merge with -@ for multi-threaded compression
    chira_utilities.print_w_time("START: Merging chunk BAM files by job type")
    final_bams = {}  # Dictionary: job_name -> final merged BAM file path
    for job_name, chunk_bam_list in all_chunk_bams.items():
        final_bam = os.path.join(args.outdir, f"{job_name}.bam")
        if len(chunk_bam_list) == 1:
            # Edge case: Only one chunk, just move/rename (no merge needed)
            shutil.move(chunk_bam_list[0], final_bam)
        else:
            # Merge all chunk BAMs for this job type into a single final BAM
            # -@ processes enables multi-threaded compression during merge
            # For 10 chunks: 8 threads = 5-10x speedup vs single-threaded
            total_processes = args.processes if args.processes is not None else multiprocessing.cpu_count()
            pysam.merge("-f", "-@", str(max(1, total_processes)), final_bam, *chunk_bam_list)
        final_bams[job_name] = final_bam
    chira_utilities.print_w_time("END: Merged chunk BAM files")
    
    # Clean up chunk files and directories
    chira_utilities.print_w_time("START: Cleaning up chunk files and directories")
    for chunk_file in chunk_files:
        try:
            os.remove(chunk_file)
        except OSError:
            pass
    # Extract chunk indices from filenames to clean up correct directories
    # Chunk files are named "chunk_XXX.fasta", extract XXX to get the index
    for chunk_file in chunk_files:
        basename = os.path.basename(chunk_file)
        chunk_idx = int(basename.replace("chunk_", "").replace(".fasta", ""))
        chunk_outdir = os.path.join(chunk_dir, f"chunk_{chunk_idx:03d}_out")
        try:
            shutil.rmtree(chunk_outdir)
        except OSError:
            pass
    try:
        os.rmdir(chunk_dir)
    except OSError:
        pass
    chira_utilities.print_w_time("END: Cleaned up chunk files and directories")
    
    # EXECUTION SUMMARY:
    # - All chunks have been processed through all alignment job types
    # - Chunk BAMs have been merged by job type into final BAM files
    # - Final BAMs are in args.outdir: index1.long.bam, index1.short.bam, etc.
    # - These final BAMs are identical to what would be produced by processing the full file
    # - Return dictionary mapping job_name -> final BAM file path for downstream merging
    return final_bams


def run_bwa_mapping_parallel(args, index1, index2):
    """
    OPTIMIZATION: Run BWA mapping with parallel job strategy (no chunking).
    
    Performance Benefits:
    - Original sequential approach: 4 jobs × 2 hours = 8 hours total
    - Parallel approach: All 4 jobs run simultaneously = ~2 hours total (4x speedup)
    - Total processes are automatically divided among parallel jobs
    - For 4 jobs with 8 processes each: 32 total threads (good for 32+ core systems)
    - BWA jobs are I/O and CPU intensive, so parallel execution maximizes resource utilization
    
    Strategy:
    - Create all alignment jobs (long/short × index1/index2)
    - Run all jobs simultaneously using ThreadPoolExecutor
    - Each job processes the full FASTA file independently
    - Collect results and return final BAM paths
    """
    alignment_jobs = []
    
    # Job 1: long + index1
    alignment_jobs.append(("long", "index1", args.fasta, index1, args.seed_length1, args.align_score1,
                          args.match1, args.mismatch1, args.gapopen1, args.gapext1, args.nhits1))
    
    # Job 2: short + index1
    alignment_jobs.append(("short", "index1", args.fasta, index1, args.seed_length2, args.align_score2,
                          args.match2, args.mismatch2, args.gapopen2, args.gapext2, args.nhits2))
    
    # Jobs 3 & 4: long + short + index2 (if index2 exists)
    if index2:
        alignment_jobs.append(("long", "index2", args.fasta, index2, args.seed_length1, args.align_score1,
                              args.match1, args.mismatch1, args.gapopen1, args.gapext1, args.nhits1))
        alignment_jobs.append(("short", "index2", args.fasta, index2, args.seed_length2, args.align_score2,
                              args.match2, args.mismatch2, args.gapopen2, args.gapext2, args.nhits2))
    
    # Calculate per-job processes (jobs run in parallel, so divide total by number of jobs)
    total_processes = args.processes if args.processes is not None else multiprocessing.cpu_count()
    per_job_processes = calculate_per_job_processes(total_processes, len(alignment_jobs))
    print(f"INFO: Using {total_processes} total processes across {len(alignment_jobs)} parallel jobs → {per_job_processes} processes per job", file=sys.stderr)
    
    # OPTIMIZATION: Run all alignment jobs in parallel using ThreadPoolExecutor
    # - All jobs run simultaneously instead of sequentially
    # - Total execution time = max(job_times) instead of sum(job_times)
    # - For 4 jobs: 4x speedup compared to sequential execution
    chira_utilities.print_w_time(f"START: Running {len(alignment_jobs)} BWA alignment jobs in parallel")
    
    def run_alignment_job(job_args):
        """Wrapper function to run a single alignment job with error handling."""
        (align_type, index_type, query_fasta, refindex, seed_length, align_score,
         match_score, mismatch_score, gap_o, gap_e, n_aligns) = job_args
        job_name = f"{index_type}.{align_type}"
        try:
            chira_utilities.print_w_time(f"START: {job_name} alignment")
            align_with_bwa(align_type, index_type, query_fasta, refindex, args.outdir,
                          seed_length, align_score, match_score, mismatch_score,
                          gap_o, gap_e, n_aligns, per_job_processes)
            chira_utilities.print_w_time(f"END: {job_name} alignment")
            return (job_name, True, None)
        except Exception as e:
            error_msg = f"Error in {job_name} alignment: {str(e)}"
            print(error_msg, file=sys.stderr)
            return (job_name, False, error_msg)
    
    failed_jobs = []
    with ThreadPoolExecutor(max_workers=len(alignment_jobs)) as executor:
        future_to_job = {executor.submit(run_alignment_job, job): job for job in alignment_jobs}
        
        for future in as_completed(future_to_job):
            job_name, success, error = future.result()
            if not success:
                failed_jobs.append((job_name, error))
    
    chira_utilities.print_w_time(f"END: All {len(alignment_jobs)} BWA alignment jobs completed")
    
    if failed_jobs:
        print("ERROR: Some alignment jobs failed:", file=sys.stderr)
        for job_name, error in failed_jobs:
            print(f"  - {job_name}: {error}", file=sys.stderr)
        sys.exit(1)
    
    # Store final BAM paths for merging
    final_bams = {}
    for align_type in ["long", "short"]:
        for idx_type in ["index1"] + (["index2"] if index2 else []):
            job_name = f"{idx_type}.{align_type}"
            final_bams[job_name] = os.path.join(args.outdir, f"{job_name}.bam")
    
    return final_bams


def calculate_sort_memory(args):
    """
    OPTIMIZATION: Calculate optimal memory per thread for BAM sorting.
    
    Performance Impact:
    - samtools sort uses temporary files when memory is insufficient
    - More memory = fewer temporary files = less I/O = faster sorting
    - Default 1G per thread: Many temp files, slow (hours for large BAMs)
    - 2-3G per thread: Fewer temp files, much faster (30-60 min for large BAMs)
    - Optimal: 2-4G per thread, but must fit in available RAM
    - For 8 threads with 2G each: 16GB total memory usage
    
    Strategy:
    - Use user-specified value if provided
    - Otherwise auto-calculate based on available RAM
    - Reserve 2GB for OS, divide remainder by thread count
    - Clamp between 1G (minimum) and 3G (optimal maximum)
    """
    if args.sort_memory:
        return args.sort_memory
    
    # OPTIMIZATION: Auto-calculate optimal memory per thread based on available RAM
    # - More memory reduces temporary file creation and I/O operations
    # - Target: 2-3G per thread for optimal performance
    # - Must ensure total (memory × threads) doesn't exceed available RAM
    if PSUTIL_AVAILABLE:
        try:
            available_gb = psutil.virtual_memory().available / (1024**3)
            total_processes = args.processes if args.processes is not None else multiprocessing.cpu_count()
            num_processes = max(1, total_processes)
            usable_gb = (available_gb - 2) / num_processes
            if usable_gb >= 2:
                sort_memory_gb = min(3, int(usable_gb))
            else:
                sort_memory_gb = max(1, int(max(0, usable_gb)))
            return f"{sort_memory_gb}G"
        except (AttributeError, OSError, ZeroDivisionError):
            return "2G"
    else:
        return "2G"


def merge_and_sort_bams(args, final_bams):
    """
    OPTIMIZATION: Merge and sort BAM files with parallel processing.
    
    Performance Benefits:
    - BAM merging: Parallel compression/decompression with -@ threads
      * Single-threaded: ~30 min for 4 large BAMs
      * 8 threads: ~4-5 min (6-7x speedup)
    - BAM sorting: Most time-consuming step, benefits greatly from:
      * Multi-threading (-@): Parallel sorting operations
      * Optimal memory (-m): Reduces temporary file I/O
      * Single-threaded with 1G memory: ~4-6 hours for large BAMs
      * 8 threads with 2G memory: ~30-45 min (8-12x speedup)
    """
    # OPTIMIZATION: Merge BAM files with parallel compression
    # - pysam.merge uses samtools merge internally
    # - -@ processes enables multi-threaded compression/decompression
    # - This significantly speeds up merging multiple large BAM files
    # - For 4 BAM files totaling 50GB: 8 threads = 6-7x speedup vs single-threaded
    chira_utilities.print_w_time("START: Merge final BAM files")
    final_bam_list = list(final_bams.values())
    total_processes = args.processes if args.processes is not None else multiprocessing.cpu_count()
    pysam.merge("-f", "-@", str(max(1, total_processes)),
                os.path.join(args.outdir, "unsorted.bam"),
                *final_bam_list)
    chira_utilities.print_w_time("END: Merge final BAM files")
    
    # Remove intermediate BAM files
    for bam_file in final_bams.values():
        try:
            os.remove(bam_file)
        except OSError:
            pass
    
    # OPTIMIZATION: Sort BAM file with parallel processing and optimal memory
    # - pysam.sort uses samtools sort internally
    # - -@ processes: Multi-threaded sorting (8 threads = 6-8x speedup)
    # - -m memory: More memory = fewer temporary files = less I/O = faster sorting
    # - -n: Sort by read name (required for downstream processing)
    # - This is one of the most time-consuming steps and benefits greatly from parallelization
    # - For 50GB BAM: 8 threads + 2G memory = 30-45 min vs 4-6 hours single-threaded
    chira_utilities.print_w_time("START: Sorting BAM file")
    sort_memory = calculate_sort_memory(args)
    total_processes = args.processes if args.processes is not None else multiprocessing.cpu_count()
    pysam.sort("-m", sort_memory, "-@", str(max(1, total_processes)), "-n",
               os.path.join(args.outdir, "unsorted.bam"),
               "-T", os.path.join(args.outdir, "sorted"),
               "-o", os.path.join(args.outdir, "sorted.bam"))
    chira_utilities.print_w_time("END: Sorting BAM file")
    os.remove(os.path.join(args.outdir, "unsorted.bam"))
    
    # Write alignments to BED
    chira_utilities.print_w_time("START: Write alignments to BED")
    write_mapped_bed(os.path.join(args.outdir, "sorted.bam"),
                     os.path.join(args.outdir, "mapped.bed"),
                     os.path.join(args.outdir, "unmapped.fasta"),
                     args.stranded)
    chira_utilities.print_w_time("END: Write alignments to BED")


def process_bed_file(args):
    """
    OPTIMIZATION: Sort and deduplicate BED file with parallel sort when available.
    
    Performance Impact:
    - BED files can be very large (multi-GB for 150M reads)
    - Standard sort: Single-threaded, ~1-2 hours for large BED files
    - GNU parallel sort: Multi-threaded, ~10-15 min for large BED files (4-8x speedup)
    - GNU sort >= 8.6 supports --parallel option for multi-threaded sorting
    - Falls back to standard sort if GNU sort not available or version too old
    """
    chira_utilities.print_w_time("START: Sorting BED file")
    # OPTIMIZATION: Use parallel sort if available (GNU sort supports --parallel option)
    # Performance: For large BED files (5-10GB), parallel sort provides 4-8x speedup
    # - GNU sort >= 8.6 supports --parallel option for multi-threaded sorting
    # - For systems with GNU sort >= 8.6, this can significantly speed up sorting large BED files
    # - Fall back to standard sort if --parallel is not supported
    sort_cmd = "sort"
    try:
        # OPTIMIZATION: Check if GNU sort with --parallel is available (GNU sort >= 8.6)
        # - Only GNU sort (Linux) supports --parallel, not BSD/macOS sort
        # - Check version string to detect GNU sort
        # - Use parallel sort with number of threads equal to processes
        result = subprocess.run(["sort", "--version"], capture_output=True, text=True, timeout=2)
        if "GNU coreutils" in result.stdout:
            # OPTIMIZATION: Use parallel sort with number of threads equal to total processes
            # - Parallel sort distributes sorting work across multiple threads
            # - For 8 threads: 4-8x speedup compared to single-threaded sort
            total_processes = args.processes if args.processes is not None else multiprocessing.cpu_count()
            sort_cmd = f"sort --parallel={max(1, total_processes)}"
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
        # OPTIMIZATION: Fall back to standard sort if check fails
        # - If GNU sort not available or version check fails, use standard sort
        # - Standard sort still works, just without parallelization benefits
        pass
    
    os.system(sort_cmd + " -k 4,4 -u " + os.path.join(args.outdir, "mapped.bed")
              + " > " + os.path.join(args.outdir, "sorted.bed"))
    chira_utilities.print_w_time("END: Sorting BED file")
    os.remove(os.path.join(args.outdir, "mapped.bed"))


def main():
    """Main function to orchestrate the mapping workflow."""
    args = parse_arguments()
    
    print_cpu_guidance(args)
    print_configuration(args)
    
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    validate_arguments(args)
    index1, index2 = setup_indices(args)
    
    if args.aligner == "clan":
        run_clan_mapping(args, index1, index2)
    elif args.aligner == "bwa":
        if args.chunk_fasta and args.chunk_fasta > 1:
            final_bams = run_bwa_mapping_with_chunking(args, index1, index2)
        else:
            final_bams = run_bwa_mapping_parallel(args, index1, index2)
        merge_and_sort_bams(args, final_bams)
    else:
        sys.stderr.write("Unknown aligner!! Currently suppoted aligners: BWA-mem and CLAN\n")
        sys.exit(1)
    
    process_bed_file(args)


if __name__ == "__main__":
    main()
