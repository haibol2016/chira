#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
import pysam
import chira_utilities

# Optional dependency for automatic memory detection
# If not available, memory calculation will use a safe default
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False


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
    # This parallelizes compression/decompression, which is especially beneficial for large files
    bwacall = ("bwa mem " + " ".join(bwa_params) + " | samtools view -hb -@ " + str(processes) + " - > " + bam)
    print(bwacall)
    os.system(bwacall)


def write_mapped_bed(bam, bed, fasta, stranded):
    """
        Extracts the mapped and unmapped reads from the BAM alignments and writes them to BED and fasta files
        Parameters:
            bam: BAM file containing all merged alignments
            bed: output BED file path
            fasta: output FASTA file path
            stranded: Strand specificity
    """
    prev_readid = None
    prev_fasta = ""
    prev_unmapped = True
    # OPTIMIZATION: Use larger buffer size (2MB) for file I/O to reduce system calls
    # This significantly improves performance for large BAM files with many small writes
    # 2MB is a good balance between memory usage and I/O efficiency
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
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

            # Only get sequence when needed (for unmapped reads)
            readseq = None
            prev_readid = readid

            if alignment.is_unmapped:
                # Only get sequence for unmapped reads
                readseq = alignment.get_forward_sequence()
                prev_fasta = f">{readid}\n{readseq}\n"
                continue
            
            # Check if primary alignment is on desired strand (cache the result)
            is_desired_strand = (stranded == "both" or 
                                (stranded == "rc" and alignment.is_reverse) or 
                                (stranded == "fw" and not alignment.is_reverse))
            
            # For stranded RNA-seq: if primary is on wrong strand, it's invalid and shouldn't be used as quality threshold
            # Only use primary length as threshold if it's on the desired strand
            # Write the alignment only if it mapped on desired strand
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
                # Check if it mapped on desired strand
                if (alt_refstrand.startswith('-') and stranded == "fw") or \
                        (alt_refstrand.startswith('+') and stranded == "rc"):
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
                    # Check if it mapped on desired strand
                    if (alt_refstrand.startswith('-') and stranded == "fw") or \
                            (alt_refstrand.startswith('+') and stranded == "rc"):
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
    # OPTIMIZATION: Use larger buffer size (2MB) for file I/O to reduce system calls
    # This significantly improves performance for large CLAN output files
    # 2MB is a good balance between memory usage and I/O efficiency
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
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


if __name__ == "__main__":

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

    parser.add_argument('-p', '--processes', action='store', type=int, default=1, metavar='',
                        dest='processes',
                        help='Number of processes to use')

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

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.4.4')

    args = parser.parse_args()
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
    print('Number of processes                  : ' + str(args.processes))
    print('Stranded                             : ' + args.stranded)
    print('Seed length                          : ' + str(args.seed_length1))
    if args.seed_length2:
        print('Seed length for 2nd iteration        : ' + str(args.seed_length2))
    print('Alignment score                      : ' + str(args.align_score1))
    if args.align_score2:
        print('Alignment score for 2nd iteration    : ' + str(args.align_score2))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    print("===================================================================")

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if not args.idx1 and not args.ref_fasta1:
        sys.stderr.write("option -x1 or -f1 are required")
        sys.exit(1)

    if (args.idx1 or args.idx2) and args.build_index:
        sys.stderr.write("options -b and -x1 are mutually exclusive")
        sys.exit(1)
    if not args.idx1 and not args.idx2 and not args.build_index:
        sys.stderr.write("Either -b or -x1 is required")
        sys.exit(1)

    index1 = index2 = None
    if args.idx1:
        index1 = args.idx1
    if args.idx2:
        index2 = args.idx2
    if args.aligner == "clan":
        if args.build_index:
            index1 = os.path.join(args.outdir, "index1")
            os.system("clan_index -f " + args.ref_fasta1 + " -d " + index1)
            if args.ref_fasta2:
                index2 = os.path.join(args.outdir, "index2")
                os.system("clan_index -f " + args.ref_fasta2 + " -d " + index2)

        # use only align_score2 for mapping with CLAN
        chira_utilities.print_w_time("START: Map read using CLAN")
        align_with_clan(args.fasta, args.outdir, args.ref_fasta1, index1, args.ref_fasta2, index2,
                        args.chimeric_overlap, args.align_score2, args.stranded, args.nhits2, args.processes)
        chira_utilities.print_w_time("END: Map read using CLAN")
        chira_utilities.print_w_time("START: Write alignments to BED")
        clan_to_bed(args.outdir)
        chira_utilities.print_w_time("END: Write alignments to BED")

    elif args.aligner == "bwa":
        # build indices
        if args.build_index:
            index1 = os.path.join(args.outdir, "index1")
            os.system("bwa index -p " + index1 + " " + args.ref_fasta1)
            if args.ref_fasta2:
                index2 = os.path.join(args.outdir, "index2")
                os.system("bwa index -p " + index2 + " " + args.ref_fasta2)

        # align with bwa
        chira_utilities.print_w_time("START: Map long read segments to index1 at " + index1)
        align_with_bwa("long", "index1", args.fasta, index1, args.outdir, args.seed_length1, args.align_score1,
                       args.match1, args.mismatch1, args.gapopen1, args.gapext1, args.nhits1, args.processes)
        chira_utilities.print_w_time("END: Map long read segments to index1 at " + index1)

        chira_utilities.print_w_time("START: Map short read segments to index1 at " + index1)
        align_with_bwa("short", "index1", args.fasta, index1, args.outdir, args.seed_length2, args.align_score2,
                       args.match2, args.mismatch2, args.gapopen2, args.gapext2, args.nhits2, args.processes)
        chira_utilities.print_w_time("END: Map short read segments to index1 at " + index1)

        if index2:
            chira_utilities.print_w_time("START: Map long read segments to index2 at " + index2)
            align_with_bwa("long", "index2", args.fasta, index2, args.outdir, args.seed_length1, args.align_score1,
                           args.match1, args.mismatch1, args.gapopen1, args.gapext1, args.nhits1, args.processes)
            chira_utilities.print_w_time("END: Map long read segments to index2 at " + index2)
            chira_utilities.print_w_time("START: Map short read segments to index2 at " + index2)
            align_with_bwa("short", "index2", args.fasta, index2, args.outdir, args.seed_length2, args.align_score2,
                           args.match2, args.mismatch2, args.gapopen2, args.gapext2, args.nhits2, args.processes)
            chira_utilities.print_w_time("END: Map short read segments to index2 at " + index2)

            chira_utilities.print_w_time("START: Merge BAM files")
            # OPTIMIZATION: pysam.merge uses multi-threading (-@) to parallelize BAM merging
            # This significantly speeds up merging multiple BAM files, especially for large datasets
            # -f to force if file already exists
            pysam.merge("-f", "-@", str(max(1, args.processes)),
                        os.path.join(args.outdir, "unsorted.bam"),
                        os.path.join(args.outdir, "index1.long.bam"),
                        os.path.join(args.outdir, "index1.short.bam"),
                        os.path.join(args.outdir, "index2.long.bam"),
                        os.path.join(args.outdir, "index2.short.bam"))
            chira_utilities.print_w_time("END: Merge BAM files")
        else:
            chira_utilities.print_w_time("START: Merge BAM files")
            # OPTIMIZATION: pysam.merge uses multi-threading (-@) to parallelize BAM merging
            # -f to force if file already exists
            pysam.merge("-f", "-@", str(max(1, args.processes)),
                        os.path.join(args.outdir, "unsorted.bam"),
                        os.path.join(args.outdir, "index1.long.bam"),
                        os.path.join(args.outdir, "index1.short.bam"))
            chira_utilities.print_w_time("END: Merge BAM files")
        # Remove intermediate BAM files
        try:
            os.remove(os.path.join(args.outdir, "index1.long.bam"))
            os.remove(os.path.join(args.outdir, "index1.short.bam"))
            if index2:
                os.remove(os.path.join(args.outdir, "index2.long.bam"))
                os.remove(os.path.join(args.outdir, "index2.short.bam"))
        except OSError:
            pass

        chira_utilities.print_w_time("START: Sorting BAM file")
        # OPTIMIZATION: pysam.sort uses multi-threading (-@) to parallelize BAM sorting
        # This is one of the most time-consuming steps and benefits greatly from parallelization
        # OPTIMIZATION: Calculate optimal memory per thread for sorting
        # More memory reduces temporary files and I/O operations, significantly improving performance
        # Rule of thumb: Use 2-4G per thread, but ensure total (memory × processes) < available RAM
        if args.sort_memory:
            sort_memory = args.sort_memory
        else:
            # Auto-calculate: try to use 2-3G per thread, but cap at available RAM
            if PSUTIL_AVAILABLE:
                try:
                    available_gb = psutil.virtual_memory().available / (1024**3)
                    # Reserve 2GB for OS and other processes, then divide by number of processes
                    # Defensive check: ensure processes > 0 to avoid division by zero
                    # Use max(1, args.processes) to ensure we always divide by at least 1
                    num_processes = max(1, args.processes)
                    usable_gb = (available_gb - 2) / num_processes
                    # Use 2-3G per thread if possible, but don't exceed available memory
                    # If usable_gb < 2, use what's available (minimum 1G) to avoid exceeding RAM
                    if usable_gb >= 2:
                        sort_memory_gb = min(3, int(usable_gb))
                    else:
                        # Use available memory, but at least 1G (samtools sort minimum)
                        # Ensure usable_gb is positive before converting to int
                        sort_memory_gb = max(1, int(max(0, usable_gb)))
                    sort_memory = f"{sort_memory_gb}G"
                except (AttributeError, OSError, ZeroDivisionError):
                    # Fallback if psutil fails at runtime or division error (shouldn't happen, but be safe)
                    sort_memory = "2G"
            else:
                # Fallback: use 2G per thread (better default than 1G)
                # This is a safe default that works well for most systems
                # Install psutil (pip install psutil) for automatic memory detection
                sort_memory = "2G"
        
        # -m sets memory per thread, -n sorts by read name
        pysam.sort("-m", sort_memory, "-@", str(max(1, args.processes)), "-n",
                   os.path.join(args.outdir, "unsorted.bam"),
                   "-T", os.path.join(args.outdir, "sorted"),
                   "-o", os.path.join(args.outdir, "sorted.bam"))
        chira_utilities.print_w_time("END: Sorting BAM file")
        os.remove(os.path.join(args.outdir, "unsorted.bam"))

        chira_utilities.print_w_time("START: Write alignments to BED")
        write_mapped_bed(os.path.join(args.outdir, "sorted.bam"),
                         os.path.join(args.outdir, "mapped.bed"),
                         os.path.join(args.outdir, "unmapped.fasta"),
                         args.stranded)
        chira_utilities.print_w_time("END: Write alignments to BED")
    else:
        sys.stderr.write("Unknown aligner!! Currently suppoted aligners: BWA-mem and CLAN")
        sys.exit(1)

    chira_utilities.print_w_time("START: Sorting BED file")
    # OPTIMIZATION: Use parallel sort if available (GNU sort supports --parallel option)
    # For systems with GNU sort >= 8.6, this can significantly speed up sorting large BED files
    # Fall back to standard sort if --parallel is not supported
    sort_cmd = "sort"
    try:
        # Check if GNU sort with --parallel is available (GNU sort >= 8.6)
        result = subprocess.run(["sort", "--version"], capture_output=True, text=True, timeout=2)
        if "GNU coreutils" in result.stdout and args.processes > 0:
            # Use parallel sort with number of threads equal to processes
            # Defensive check: ensure processes > 0 (default is 1, so this should always be true)
            sort_cmd = f"sort --parallel={max(1, args.processes)}"
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
        # Fall back to standard sort if check fails
        pass
    
    os.system(sort_cmd + " -k 4,4 -u " + os.path.join(args.outdir, "mapped.bed")
              + " > " + os.path.join(args.outdir, "sorted.bed"))
    chira_utilities.print_w_time("END: Sorting BED file")
    os.remove(os.path.join(args.outdir, "mapped.bed"))
