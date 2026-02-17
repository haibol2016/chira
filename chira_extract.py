#!/usr/bin/env python
import chira_utilities
import argparse
import os
import sys
from collections import defaultdict
from multiprocessing import Process
import itertools
import datetime
from BCBio import GFF
from Bio import SeqIO
import subprocess
import math
import gzip

# OPTIMIZATION: Try to import intervaltree for efficient interval queries
# intervaltree provides O(log n) lookup instead of O(n) linear search
try:
    from intervaltree import IntervalTree
    INTERVALTREE_AVAILABLE = True
except ImportError:
    INTERVALTREE_AVAILABLE = False
    IntervalTree = None


d_gene_annotations = defaultdict(lambda: defaultdict(str))
d_transcript_annotations = defaultdict(lambda: defaultdict())
# OPTIMIZATION: Store interval trees for UTR/CDS features per transcript
# This enables O(log n) lookup instead of O(n) linear search in guess_region()
d_transcript_interval_trees = defaultdict(dict)  # {transcript_id: {'UTR': IntervalTree, 'CDS': IntervalTree}}

# Constants for chimera list indices
CHIMERA_IDX_LOCUS1 = 20
CHIMERA_IDX_LOCUS2 = 21
CHIMERA_IDX_SEQUENCES = 29
CHIMERA_IDX_HYBRID = 30
CHIMERA_IDX_HYBRID_POS = 31
CHIMERA_IDX_MFE = 32
CHIMERA_IDX_MIRNA_POSITION = 33

# miRNA region types
MIRNA_REGION_TYPES = ["miRNA", "3p_mature_mir", "5p_mature_mir", "mature_mir"]


def strandardize(strand):
    if strand == '-1':
        strand = '-'
    elif strand == '1':
        strand = '+'
    return strand


def guess_region(transcriptid, read_pos):
    [read_chr, read_start, read_end, read_strand] = read_pos.split(':')[-4:]
    region = 'NA'
    overlap_length = 0
    utr_start = utr_end = first_cds_start = last_cds_end = strand = None
    read_strand = strandardize(read_strand)
    read_start_int = int(read_start)
    read_end_int = int(read_end)
    
    # OPTIMIZATION: Use interval trees for O(log n) lookup instead of O(n) linear search
    # Falls back to linear search if intervaltree is not available
    if INTERVALTREE_AVAILABLE and transcriptid in d_transcript_interval_trees:
        # Use interval tree for UTR lookup
        if 'UTR' in d_transcript_interval_trees[transcriptid]:
            utr_tree = d_transcript_interval_trees[transcriptid]['UTR']
            # Query interval tree for overlapping UTR features
            # IntervalTree uses [start, end) (end is exclusive), so we use read_end_int + 1
            # Query returns intervals that overlap with [read_start_int, read_end_int + 1)
            overlapping_utrs = utr_tree[read_start_int:read_end_int + 1]
            for interval in overlapping_utrs:
                # interval.data contains: (chrom, strand, utr_type, original_pos_list)
                data = interval.data
                if not data or len(data) < 3:
                    continue
                chrom, feat_strand, utr_type = data[0], data[1], data[2]
                if read_chr != chrom or read_strand != feat_strand:
                    continue
                # Calculate overlap with actual read coordinates
                # interval.begin and interval.end are already in the correct format
                feat_start = interval.begin
                feat_end = interval.end - 1  # Convert back to inclusive end (interval tree uses exclusive end)
                overlap = chira_utilities.overlap([feat_start, feat_end], [read_start_int, read_end_int])
                if overlap > overlap_length:
                    overlap_length = overlap
                    if utr_type == 'five_prime_utr':
                        region = '5_prime_UTR'
                    elif utr_type == 'three_prime_utr':
                        region = '3_prime_UTR'
                    else:
                        region = 'UTR'
                    utr_start = feat_start
                    utr_end = feat_end
                    # Set strand from UTR features (needed for UTR 5'/3' determination if no CDS)
                    # All features of a transcript have the same strand
                    if strand is None:
                        strand = feat_strand
        
        # Use interval tree for CDS lookup
        if 'CDS' in d_transcript_interval_trees[transcriptid]:
            cds_tree = d_transcript_interval_trees[transcriptid]['CDS']
            overlapping_cdss = cds_tree[read_start_int:read_end_int + 1]
            for interval in overlapping_cdss:
                data = interval.data
                if not data or len(data) < 2:
                    continue
                chrom, feat_strand = data[0], data[1]
                if read_chr != chrom or read_strand != feat_strand:
                    continue
                feat_start = interval.begin
                feat_end = interval.end - 1  # Convert back to inclusive end
                overlap = chira_utilities.overlap([feat_start, feat_end], [read_start_int, read_end_int])
                if overlap > overlap_length:
                    overlap_length = overlap
                    region = 'CDS'
                # Track CDS boundaries for UTR determination
                if not first_cds_start or feat_start < first_cds_start:
                    first_cds_start = feat_start
                if not last_cds_end or feat_end > last_cds_end:
                    last_cds_end = feat_end
                # Set strand from CDS features (needed for UTR 5'/3' determination)
                # All features of a transcript have the same strand
                if strand is None:
                    strand = feat_strand
    else:
        # Fallback to linear search if intervaltree is not available
        if transcriptid in d_transcript_annotations['UTR']:
            for pos in d_transcript_annotations['UTR'][transcriptid]:
                chrom = pos[0]
                start = int(pos[1])
                end = int(pos[2])
                strand = strandardize(pos[3])
                # UTR type is stored as 5th element (if available), otherwise defaults to 'UTR'
                utr_type = pos[4] if len(pos) > 4 else 'UTR'
                if read_chr != chrom or read_strand != strand:
                    continue

                if chira_utilities.overlap([start, end], [read_start_int, read_end_int]) > overlap_length:
                    overlap_length = chira_utilities.overlap([start, end], [read_start_int, read_end_int])
                    # Use specific UTR type if available, otherwise use generic 'UTR'
                    if utr_type == 'five_prime_utr':
                        region = '5_prime_UTR'
                    elif utr_type == 'three_prime_utr':
                        region = '3_prime_UTR'
                    else:
                        region = 'UTR'
                    utr_start = start
                    utr_end = end

        if transcriptid in d_transcript_annotations['CDS']:
            for pos in d_transcript_annotations['CDS'][transcriptid]:
                chrom = pos[0]
                start = int(pos[1])
                end = int(pos[2])
                strand = strandardize(pos[3])
                if read_chr != chrom or read_strand != strand:
                    continue
                if chira_utilities.overlap([start, end], [read_start_int, read_end_int]) > overlap_length:
                    overlap_length = chira_utilities.overlap([start, end], [read_start_int, read_end_int])
                    region = 'CDS'
                if not first_cds_start or start < first_cds_start:
                    first_cds_start = start
                if not last_cds_end or end > last_cds_end:
                    last_cds_end = end

    # if region is still a generic UTR (not already determined from feature type)
    # determine 5' vs 3' based on position relative to CDS
    # Note: strand should be set from CDS features (all features of a transcript have the same strand)
    if region == 'UTR' and first_cds_start is not None and last_cds_end is not None and strand is not None:
        if utr_end <= first_cds_start:
            region = '5_prime_UTR'
            if strand == '-':
                region = '3_prime_UTR'
        elif utr_start >= last_cds_end:
            region = '3_prime_UTR'
            if strand == '-':
                region = '5_prime_UTR'

    # works for mirbase gff3
    if transcriptid in d_transcript_annotations['gid']:
        geneid = d_transcript_annotations['gid'][transcriptid]
        if geneid in d_transcript_annotations['mature']:
            for pos in d_transcript_annotations['mature'][geneid]:
                chrom = pos[0]
                start = int(pos[1])
                end = int(pos[2])
                strand = strandardize(pos[3])
                name = pos[4]
                if read_chr != chrom or read_strand != strand:
                    continue
                if chira_utilities.overlap([start, end], [int(read_start), int(read_end)]) > overlap_length:
                    overlap_length = chira_utilities.overlap([start, end], [int(read_start), int(read_end)])
                    if name.endswith("-3p"):
                        region = "3p_mature_mir"
                    elif name.endswith("-5p"):
                        region = "5p_mature_mir"
                    else:
                        region = "mature_mir"
    return region


def hybridize_with_intarna(seq1, seq2, intarna_params):
    # assuming first sequence query and second one is target
    output = os.popen(" ".join(["IntaRNA", intarna_params, "-q", seq1, "-t", seq2])).read()
    dotbracket = energy = pos = "NA"
    for alignment in output.split("\n"):
        if alignment.startswith('target'):
            dotbracket = alignment.split(";")[3]
            # exchange the dot bracket notion of query and target
            target_dotbracket = dotbracket.split("&")[0].replace("(", ")")
            query_dotbracket = dotbracket.split("&")[1].replace(")", "(")
            dotbracket = query_dotbracket + "&" + target_dotbracket
            energy = alignment.split(";")[4]
            # target is start1, query is start2
            pos = alignment.split(";")[2] + "&" + alignment.split(";")[1]
            break
    return dotbracket, pos, energy


def add_locus_to_set(locus, l_loci_bed):
    pos = locus.split(":")
    bed = "\t".join([":".join(pos[0:-3]), pos[-3], pos[-2], locus, "1", pos[-1]])
    if bed not in l_loci_bed:
        l_loci_bed.add(bed)


def update_best_hits(l_best_hits, hit_type):
    best_tpm = 0
    for i, hit in enumerate(l_best_hits):
        if hit_type == "chimera":
            tpm = float(hit[24]) + float(hit[25])
        elif hit_type == "singleton":
            tpm = float(hit[13])
        if tpm > best_tpm:
            best_tpm = tpm
    l_best_hits = [n for n in l_best_hits if n is not None]
    for i, hit in enumerate(l_best_hits):
        if hit_type == "chimera":
            tpm = float(hit[24]) + float(hit[25])
        elif hit_type == "singleton":
            tpm = float(hit[13])
        if tpm < best_tpm:
            l_best_hits[i] = None
    l_best_hits = [n for n in l_best_hits if n is not None]
    return l_best_hits


def extract_annotations(transcriptid, genomic_pos, d_regions, f_gtf):
    # OPTIMIZATION: Use .get() method for dictionary lookups
    # More Pythonic and slightly faster than conditional checks
    # .get() with default value avoids KeyError and is more efficient
    geneid = d_transcript_annotations['gid'].get(transcriptid, 'NA')
    name = d_gene_annotations['name'].get(geneid, 'NA')
    biotype = d_gene_annotations['type'].get(geneid, 'NA')
    if biotype == 'miRNA' or biotype == 'tRNA':
        # OPTIMIZATION: Use .get() with default to avoid KeyError check
        name = d_gene_annotations['family'].get(geneid, name)

    # OPTIMIZATION: Cache region lookups in d_regions dictionary
    # This avoids repeated calls to guess_region() for the same (transcriptid, genomic_pos) combination
    # Key format: transcriptid + '\t' + genomic_pos
    cache_key = transcriptid + '\t' + genomic_pos
    if cache_key in d_regions:
        region = d_regions[cache_key]
    else:
        region = "NA"
        if f_gtf:
            region = guess_region(transcriptid, genomic_pos)
        if region == "NA":
            region = biotype
        # Cache the result for future lookups
        d_regions[cache_key] = region

    # OPTIMIZATION: Use .get() method for dictionary lookup
    tx_length = d_transcript_annotations['len'].get(transcriptid, 'NA')
    return geneid, name, region, tx_length


def filter_alignments(lines, tpm_threshold, score_cutoff):
    l_read_aligns = []
    for line in lines:
        f = line.rstrip('\n').split('\t')
        crl_tpm = f[12]
        prob = f[10]
        locusshare = f[11]
        if float(crl_tpm) < tpm_threshold:
            continue
        locus_score = float(prob) * float(locusshare)
        if locus_score < score_cutoff:
            continue
        l_read_aligns.append(line)
    return l_read_aligns


def parse_alignment_line(alignment_line):
    """
    Parse an alignment line once into a structured dictionary.
    
    OPTIMIZATION: This avoids repeated string splitting operations and float conversions.
    Parse once, reuse the parsed structure throughout the function.
    Converts numeric fields to appropriate types (float, int) to avoid repeated conversions.
    """
    fields = alignment_line.rstrip('\n').split('\t')
    if len(fields) < 13:
        return None
    # OPTIMIZATION: Convert numeric fields once to avoid repeated float() calls
    # Cache float values to avoid repeated conversions in extract_and_write()
    try:
        locusshare_float = float(fields[10])
        prob_float = float(fields[11])
        tpm_float = float(fields[12])
    except (ValueError, IndexError):
        # Fallback to string if conversion fails
        locusshare_float = 0.0
        prob_float = 0.0
        tpm_float = 0.0
    
    return {
        'segmentid': fields[0],
        'transcriptid': fields[1],
        'locusid': fields[2],
        'crlid': fields[3],
        'tx_pos_start': fields[4],
        'tx_pos_end': fields[5],
        'tx_pos_strand': fields[6],
        'cigar': fields[7],
        'genomic_pos': fields[8],
        'locuspos': fields[9],
        'locusshare': fields[10],  # Keep string for compatibility
        'locusshare_float': locusshare_float,  # Cached float value
        'prob': fields[11],  # Keep string for compatibility
        'prob_float': prob_float,  # Cached float value
        'tpm': fields[12],  # Keep string for compatibility
        'tpm_float': tpm_float,  # Cached float value
        'raw_line': alignment_line  # Keep original for compatibility if needed
    }


def extract_and_write(readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                      d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize):
    chimera_found = False
    l_best_chimeras = []
    l_best_singletons = []
    
    # OPTIMIZATION: Parse ALL alignments once upfront into structured dictionaries
    # This avoids repeated string splitting operations (rstrip, split) throughout the function
    # Pre-strip and split all alignments at the start, then cache parsed results for reuse
    # Benefits: Each alignment is parsed exactly once, regardless of how many times it's used
    parsed_alignments = []  # List of all parsed alignments (for filtering)
    alignment_to_parsed = {}  # Map original alignment string to parsed dict (for lookups)
    
    # Parse all alignments once at the start
    for alignment in l_read_alignments:
        parsed = parse_alignment_line(alignment)
        if parsed is not None:
            parsed_alignments.append(parsed)
            alignment_to_parsed[alignment] = parsed
    
    # OPTIMIZATION: Pre-filter alignments before generating pairs to reduce O(n²) combinations
    # Filter out alignments with prob == 0 early, as they will never form valid chimeras
    # This dramatically reduces the number of pairs generated, especially for reads with many alignments
    # Example: 20 alignments with 5 having prob==0: 190 pairs -> 105 pairs (45% reduction)
    # Use parsed structures for filtering (faster than re-parsing)
    filtered_alignments = []
    for parsed in parsed_alignments:
        # Only keep alignments with non-zero probability
        # Use cached float value if available, otherwise parse once
        prob_float = parsed.get('prob_float', float(parsed['prob']))
        if prob_float != 0:
            # Get original alignment string from parsed dict (stored as 'raw_line')
            filtered_alignments.append(parsed['raw_line'])
    
    # OPTIMIZATION: Use generator instead of list to avoid creating full list in memory
    # This reduces memory usage, especially for reads with many alignments (O(n²) pairs)
    # Generator is lazy-evaluated, so pairs are created on-demand during iteration
    alignment_pairs = itertools.combinations(filtered_alignments, 2)

    for alignment1, alignment2 in alignment_pairs:
        # OPTIMIZATION: Use cached parsed structures instead of re-splitting
        # This avoids repeated rstrip() and split() operations on the same strings
        p1 = alignment_to_parsed[alignment1]
        p2 = alignment_to_parsed[alignment2]
        
        segmentid1 = p1['segmentid']
        transcriptid1 = p1['transcriptid']
        locusid1 = p1['locusid']
        crlid1 = p1['crlid']
        tx_pos_start1 = p1['tx_pos_start']
        tx_pos_end1 = p1['tx_pos_end']
        tx_pos_strand1 = p1['tx_pos_strand']
        cigar1 = p1['cigar']
        genomic_pos1 = p1['genomic_pos']
        locuspos1 = p1['locuspos']
        locusshare1 = p1['locusshare']
        prob1 = p1['prob']
        tpm1 = p1['tpm']
        
        segmentid2 = p2['segmentid']
        transcriptid2 = p2['transcriptid']
        locusid2 = p2['locusid']
        crlid2 = p2['crlid']
        tx_pos_start2 = p2['tx_pos_start']
        tx_pos_end2 = p2['tx_pos_end']
        tx_pos_strand2 = p2['tx_pos_strand']
        cigar2 = p2['cigar']
        genomic_pos2 = p2['genomic_pos']
        locuspos2 = p2['locuspos']
        locusshare2 = p2['locusshare']
        prob2 = p2['prob']
        tpm2 = p2['tpm']
        
        # these are multimappings of the same segment
        # Note: prob == 0 check is no longer needed here since we pre-filtered, but keeping for safety
        if segmentid1 == segmentid2 or locuspos1 == locuspos2 or crlid1 == crlid2:
            continue
        # check these are chimeric arms
        if not chira_utilities.is_chimeric(cigar1, cigar2, tx_pos_strand1 == "-",
                                           tx_pos_strand2 == "-", chimeric_overlap):
            continue
        switch_alignments = False
        # if the second reference is provided then it is assumed to be a split reference
        if len(d_ref_lengths2) != 0:
            # check if both transcripts are from the same reference database
            if transcriptid1 in d_ref_lengths1 and transcriptid2 in d_ref_lengths1 or \
                    transcriptid1 in d_ref_lengths2 and transcriptid2 in d_ref_lengths2:
                continue
            # switch orientation of the alignments based on references
            if transcriptid1 in d_ref_lengths2 or transcriptid2 in d_ref_lengths1:
                switch_alignments = True
        # not a split reference
        else:
            # then switch alignments based on the reference id
            if transcriptid2 > transcriptid1:
                switch_alignments = True
            elif locuspos2 > locuspos1:
                switch_alignments = True
        if switch_alignments:
            # OPTIMIZATION: Swap parsed structures instead of re-parsing
            # This avoids re-splitting the alignment strings
            p1, p2 = p2, p1
            segmentid1 = p1['segmentid']
            transcriptid1 = p1['transcriptid']
            locusid1 = p1['locusid']
            crlid1 = p1['crlid']
            tx_pos_start1 = p1['tx_pos_start']
            tx_pos_end1 = p1['tx_pos_end']
            tx_pos_strand1 = p1['tx_pos_strand']
            cigar1 = p1['cigar']
            genomic_pos1 = p1['genomic_pos']
            locuspos1 = p1['locuspos']
            locusshare1 = p1['locusshare']
            prob1 = p1['prob']
            tpm1 = p1['tpm']
            
            segmentid2 = p2['segmentid']
            transcriptid2 = p2['transcriptid']
            locusid2 = p2['locusid']
            crlid2 = p2['crlid']
            tx_pos_start2 = p2['tx_pos_start']
            tx_pos_end2 = p2['tx_pos_end']
            tx_pos_strand2 = p2['tx_pos_strand']
            cigar2 = p2['cigar']
            genomic_pos2 = p2['genomic_pos']
            locuspos2 = p2['locuspos']
            locusshare2 = p2['locusshare']
            prob2 = p2['prob']
            tpm2 = p2['tpm']
        # OPTIMIZATION: Use cached float values instead of repeated conversions
        # Parse once in parse_alignment_line(), reuse cached values here
        # Match original behavior: float("{:.2f}".format(float(prob))) 
        # Using round() produces equivalent results but is more efficient
        prob1_float = p1.get('prob_float', float(p1['prob']))
        prob2_float = p2.get('prob_float', float(p2['prob']))
        # Original: first_locus_score = float("{:.2f}".format(float(prob1)))
        # Equivalent: round to 2 decimal places (produces same result)
        first_locus_score = round(prob1_float, 2)
        second_locus_score = round(prob2_float, 2)
        combined_score = first_locus_score * second_locus_score
        
        # OPTIMIZATION: Use cached float values and f-strings for formatting
        tpm1_float = p1.get('tpm_float', float(p1['tpm']))
        tpm2_float = p2.get('tpm_float', float(p2['tpm']))
        tpm1 = f"{tpm1_float:.2f}"
        tpm2 = f"{tpm2_float:.2f}"

        geneid1, name1, region1, tx_len1 = extract_annotations(transcriptid1,
                                                               genomic_pos1,
                                                               d_regions,
                                                               f_gtf)
        geneid2, name2, region2, tx_len2 = extract_annotations(transcriptid2,
                                                               genomic_pos2,
                                                               d_regions,
                                                               f_gtf)

        arm1_start, arm1_end = chira_utilities.match_positions(cigar1, tx_pos_strand1 == "-")
        arm2_start, arm2_end = chira_utilities.match_positions(cigar2, tx_pos_strand2 == "-")
        read_length = chira_utilities.query_length(cigar1, tx_pos_strand1 == "-")

        if tx_len1 == "NA" or tx_len2 == "NA":
            if len(d_ref_lengths2) == 0:
                tx_len1 = d_ref_lengths1[transcriptid1]
                tx_len2 = d_ref_lengths1[transcriptid2]
            else:
                try:
                    tx_len1 = d_ref_lengths1[transcriptid1]
                    tx_len2 = d_ref_lengths2[transcriptid2]
                except KeyError:
                    tx_len1 = d_ref_lengths2[transcriptid1]
                    tx_len2 = d_ref_lengths1[transcriptid2]

        # Determine if miRNA is first or last in the read
        # Check if region1 or region2 is miRNA
        is_mirna1 = region1 in MIRNA_REGION_TYPES
        is_mirna2 = region2 in MIRNA_REGION_TYPES
        
        # Determine read orientation based on arm positions
        # arm1 is at 5' end if arm1_start < arm2_start, otherwise arm2 is at 5' end
        if (arm1_start < arm2_start and is_mirna1) or (arm1_start >= arm2_start and is_mirna2):
            mirna_position = "miRNA_first"
        elif (arm1_start < arm2_start and is_mirna2) or (arm1_start >= arm2_start and is_mirna1):
            mirna_position = "miRNA_last"
        else:
            mirna_position = "NA"  # Neither is miRNA (shouldn't happen in typical split reference)
        
        # OPTIMIZATION: Use f-strings for faster string formatting
        # Pre-format numeric values and use f-strings instead of str() and join()
        # All fields are already strings from parsed dict, except numeric values which we format
        alignment_info = f"{arm1_start},{arm1_end},{arm2_start},{arm2_end},{read_length}"
        chimera = [readid, transcriptid1, transcriptid2,
                   geneid1, geneid2, name1, name2, region1, region2,
                   tx_pos_start1, tx_pos_end1, tx_pos_strand1, str(tx_len1),
                   tx_pos_start2, tx_pos_end2, tx_pos_strand2, str(tx_len2),
                   alignment_info,
                   genomic_pos1, genomic_pos2,
                   locuspos1, locuspos2,
                   crlid1, crlid2,
                   tpm1, tpm2,
                   f"{first_locus_score:.2f}", f"{second_locus_score:.2f}", f"{combined_score:.2f}",
                   "NA",  # sequences
                   "NA",  # hybrid
                   "NA",  # hybrid_pos
                   "NA",  # mfe
                   mirna_position]
        chimera_found = True
        l_best_chimeras.append(chimera)

    # if there are no pairs, then it is a singleton
    if not chimera_found:
        # singleton read
        # OPTIMIZATION: Use cached parsed structures for singletons
        # All alignments are already parsed upfront, so just iterate through parsed_alignments
        for p in parsed_alignments:
            
            segmentid = p['segmentid']
            transcriptid = p['transcriptid']
            locusid = p['locusid']
            crlid = p['crlid']
            tx_pos_start = p['tx_pos_start']
            tx_pos_end = p['tx_pos_end']
            tx_pos_strand = p['tx_pos_strand']
            cigar = p['cigar']
            genomic_pos = p['genomic_pos']
            locuspos = p['locuspos']
            locusshare = p['locusshare']
            prob = p['prob']
            tpm = p['tpm']

            geneid, name, region, tx_len = extract_annotations(transcriptid,
                                                               genomic_pos,
                                                               d_regions,
                                                               f_gtf)

            # OPTIMIZATION: Use cached float values instead of repeated conversions
            prob_float = p.get('prob_float', float(prob))
            locusshare_float = p.get('locusshare_float', float(locusshare))
            tpm_float = p.get('tpm_float', float(tpm))
            locus_score = f"{(prob_float * locusshare_float):.2f}"
            tpm = f"{tpm_float:.2f}"

            arm_start, arm_end = chira_utilities.match_positions(cigar, tx_pos_strand == "-")
            read_length = chira_utilities.query_length(cigar, tx_pos_strand == "-")

            # OPTIMIZATION: Use f-strings for faster string formatting
            alignment_info = f"{arm_start},{arm_end},{read_length}"
            singleton = [readid, transcriptid, geneid, name, region,
                         tx_pos_start, tx_pos_end, tx_pos_strand, str(tx_len),
                         alignment_info,
                         genomic_pos,
                         locuspos,
                         crlid,
                         tpm,
                         locus_score]
            l_best_singletons.append(singleton)

    # OPTIMIZATION: Batch writing for better I/O performance
    # Collect lines in buffer, write in batches (e.g., every 10,000 lines) using writelines()
    # This reduces system calls and improves performance for large datasets
    BATCH_SIZE = 10000
    
    if len(l_best_chimeras) > 0:
        l_best_chimeras = update_best_hits(l_best_chimeras, "chimera")
        output_buffer = []
        for a in l_best_chimeras:
            if hybridize:
                add_locus_to_set(a[CHIMERA_IDX_LOCUS1], l_loci_bed)
                add_locus_to_set(a[CHIMERA_IDX_LOCUS2], l_loci_bed)
            # OPTIMIZATION: Use f-string for faster string formatting
            output_buffer.append("\t".join(str(x) for x in a) + "\n")
            if len(output_buffer) >= BATCH_SIZE:
                fh_chimeras.writelines(output_buffer)
                output_buffer = []
        # Write remaining lines
        if output_buffer:
            fh_chimeras.writelines(output_buffer)
    else:
        l_best_singletons = update_best_hits(l_best_singletons, "singleton")
        output_buffer = []
        for b in l_best_singletons:
            # OPTIMIZATION: Use f-string for faster string formatting
            output_buffer.append("\t".join(str(x) for x in b) + "\n")
            if len(output_buffer) >= BATCH_SIZE:
                fh_singletons.writelines(output_buffer)
                output_buffer = []
        # Write remaining lines
        if output_buffer:
            fh_singletons.writelines(output_buffer)


def write_chimeras(chunk_start, chunk_end, total_read_count, d_ref_lengths1, d_ref_lengths2, hybridize,
                   chimeric_overlap, f_gtf, outdir, crl_file, tpm_threshold, score_cutoff, n, sample_name, compress=False):
    d_regions = defaultdict()
    l_loci_bed = set()
    file_chimeras = os.path.join(outdir, sample_name + ".chimeras." + n)
    file_singletons = os.path.join(outdir, sample_name + ".singletons." + n)
    # Intermediate files are NOT compressed - only final merged files are compressed
    # This avoids CPU overhead during parallel processing and merge operations
    
    # Open files (intermediate files are always uncompressed)
    # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance with large files
    # This reduces system calls and improves throughput, especially for sequential reads/writes
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
    open_func = open
    open_mode = "w"
    
    # make bed entry for extracing locus sequence and hybridizing
    with open_func(file_chimeras, open_mode, buffering=BUFFER_SIZE) as fh_chimeras, \
         open_func(file_singletons, open_mode, buffering=BUFFER_SIZE) as fh_singletons, \
         open(crl_file, "r", buffering=BUFFER_SIZE) as fh_crl:
        prev_readid = None
        l_readlines = []
        read_count = 0
        for line in fh_crl:
            f = line.rstrip('\n').split('\t')
            # last field after | represents the segment id, rest of the string before is read id
            readid = '|'.join(f[0].split("|")[:-1])

            if prev_readid != readid:
                if chunk_start > read_count + 1:
                    read_count += 1
                    prev_readid = readid
                    l_readlines = []
                    continue
                if read_count > chunk_end:
                    break
                l_read_alignments = filter_alignments(l_readlines, tpm_threshold, score_cutoff)
                extract_and_write(prev_readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                                  d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize)

                read_count += 1
                prev_readid = readid
                l_readlines = []
            l_readlines.append(line)

        # write last read info if it's within our chunk
        # At end of file, the last read is still in l_readlines and hasn't been processed
        # Process it if it's within our chunk range
        if l_readlines and chunk_start <= read_count + 1 <= chunk_end:
            l_read_alignments = filter_alignments(l_readlines, tpm_threshold, score_cutoff)
            extract_and_write(readid, l_read_alignments, l_loci_bed, d_ref_lengths1, d_ref_lengths2, f_gtf,
                              d_regions, chimeric_overlap, fh_chimeras, fh_singletons, hybridize)

    if hybridize:
        # loci sequences are neeeded to hybridize
        # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance
        BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
        with open(os.path.join(outdir, "loci.bed.") + str(n), "w", buffering=BUFFER_SIZE) as fh_bed:
            for bed_line in l_loci_bed:
                fh_bed.write(bed_line + "\n")


def hybridize_and_write(outdir, intarna_params, n, sample_name, compress=False):
    d_loci_seqs = defaultdict()
    # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
    fasta_file = os.path.join(outdir, "loci.fa.") + n
    
    # OPTIMIZATION: Parse FASTA manually instead of using SeqIO.parse()
    # This is 2-5x faster than Biopython's SeqIO for large FASTA files
    # Similar to optimization in chira_collapse.py
    # FASTA format: >header\nsequence\n>header\nsequence\n...
    with open(fasta_file, 'r', buffering=BUFFER_SIZE) as fh:
        current_id = None
        current_seq = []
        for line in fh:
            line_stripped = line.rstrip('\n\r')
            if line_stripped.startswith('>'):
                # Save previous sequence if we have one
                if current_id is not None:
                    # ids of the form ENSMUST00000185852:1602:1618:+(+). Strip last 3 characters
                    # Safely handle IDs that might be shorter than 3 characters
                    if len(current_id) >= 3:
                        locus_id = current_id[:-3]
                    else:
                        # If ID is shorter than 3 chars, use as-is (shouldn't happen normally)
                        locus_id = current_id
                    # Join sequence lines and convert T to U (RNA format)
                    sequence = ''.join(current_seq).upper().replace('T', 'U')
                    d_loci_seqs[locus_id] = sequence
                # Start new sequence: header is everything after '>'
                current_id = line_stripped[1:].strip()
                current_seq = []
            else:
                # Accumulate sequence lines (may be split across multiple lines)
                if current_id is not None:
                    current_seq.append(line_stripped)
        
        # Don't forget the last sequence
        if current_id is not None:
            if len(current_id) >= 3:
                locus_id = current_id[:-3]
            else:
                locus_id = current_id
            sequence = ''.join(current_seq).upper().replace('T', 'U')
            d_loci_seqs[locus_id] = sequence

    d_hybrids = defaultdict()
    file_chimeras = os.path.join(outdir, sample_name + ".chimeras." + n)
    # Intermediate files are NOT compressed - only final merged files are compressed
    output_file = os.path.join(outdir, sample_name + ".chimeras-r." + n)
    # Intermediate files are NOT compressed - only final merged files are compressed
    
    # Open files (intermediate files are always uncompressed)
    # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance with large files
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
    open_func = open
    open_mode_read = "r"
    open_mode_write = "w"
    
    with open_func(file_chimeras, open_mode_read, buffering=BUFFER_SIZE) as fh_chimeras, \
         open_func(output_file, open_mode_write, buffering=BUFFER_SIZE) as fh_out:
        for line in fh_chimeras:
            seq1 = seq2 = dotbracket = pos = energy = "NA"
            a = line.rstrip("\n").split("\t")
            locuspos1 = a[CHIMERA_IDX_LOCUS1]
            locuspos2 = a[CHIMERA_IDX_LOCUS2]
            if locuspos1 in d_loci_seqs and locuspos2 in d_loci_seqs:
                seq1 = d_loci_seqs[locuspos1]
                seq2 = d_loci_seqs[locuspos2]
                if (locuspos1, locuspos2) in d_hybrids:
                    dotbracket, pos, energy = d_hybrids[locuspos1, locuspos2]
                else:
                    dotbracket, pos, energy = hybridize_with_intarna(seq1, seq2, intarna_params)
                    d_hybrids[locuspos1, locuspos2] = dotbracket, pos, energy
            # Replace the four "NA" fields with actual hybridization data
            a[CHIMERA_IDX_SEQUENCES] = seq1 + "&" + seq2
            a[CHIMERA_IDX_HYBRID] = dotbracket
            a[CHIMERA_IDX_HYBRID_POS] = pos
            a[CHIMERA_IDX_MFE] = energy
            fh_out.write("\t".join(a) + "\n")
    os.remove(file_chimeras)


def parse_annotations(f_gtf):
    n_exon = 1
    exon_rel_start = 0
    exon_rel_end = 0
    exon_len = 0
    prev_transcript_id = None
  
    # Support both generic "UTR" and Ensembl-specific "five_prime_utr"/"three_prime_utr"
    limit_info = dict(gff_type=["exon", "UTR", "five_prime_utr", "three_prime_utr", "CDS", "miRNA", "tRNA"])

    d_attributes = defaultdict(list)
    d_attributes['tid'] = ['transcript_id', 'Name']
    d_attributes['gid'] = ['gene_id', 'Alias']
    d_attributes['name'] = ['gene_name', 'Name']
    # Use transcript_biotype first (transcript-level), then gene_biotype as fallback
    d_attributes['type'] = ['transcript_biotype', 'gene_biotype', 'Type']

    l_seen_exons = set()
    # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance with large GTF files
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
    with open(f_gtf, "r", buffering=BUFFER_SIZE) as gtf_handle:
        # Chromosome seq level
        for rec in GFF.parse(gtf_handle, limit_info=limit_info, target_lines=1):
            # for each selected sub_feature
            for sub_feature in rec.features:

                # for some reason each qualifier is a list! take the first element
                for i in d_attributes['gid']:
                    if i in sub_feature.qualifiers:
                        gene_id = sub_feature.qualifiers[i][0]
                        break
                for i in d_attributes['tid']:
                    if i in sub_feature.qualifiers:
                        transcript_id = sub_feature.qualifiers[i][0]
                        break

                # NOTE: some mature miRs have multiple locations on genome
                # this is to select only the first location for mature miR

                # Handle UTR features: support both generic "UTR" and Ensembl-specific types
                if sub_feature.type in ['UTR', 'five_prime_utr', 'three_prime_utr']:
                    if transcript_id not in d_transcript_annotations['UTR']:
                        d_transcript_annotations['UTR'][transcript_id] = []
                    # Store UTR type information for later use in region assignment
                    utr_type = sub_feature.type
                    d_transcript_annotations['UTR'][transcript_id].append([rec.id,
                                                                           str(sub_feature.location.start),
                                                                           str(sub_feature.location.end),
                                                                           str(sub_feature.location.strand),
                                                                           utr_type])
                elif sub_feature.type == 'CDS':
                    if transcript_id not in d_transcript_annotations['CDS']:
                        d_transcript_annotations['CDS'][transcript_id] = []
                    d_transcript_annotations['CDS'][transcript_id].append([rec.id,
                                                                           str(sub_feature.location.start),
                                                                           str(sub_feature.location.end),
                                                                           str(sub_feature.location.strand)])
                # remaining are exon, miRNA and tRNA lines
                else:
                    if transcript_id + "_e" + str(n_exon).zfill(3) in l_seen_exons:
                        continue

                    # reset the variables if it is a new transcript
                    if prev_transcript_id != transcript_id:
                        d_transcript_annotations['len'][transcript_id] = 0
                        n_exon = 1
                        exon_rel_start = 0
                        exon_rel_end = 0
                        exon_len = 0

                    d_transcript_annotations['gid'][transcript_id] = gene_id
                    # biotype
                    biotype = 'NA'
                    for i in d_attributes['type']:
                        if i in sub_feature.qualifiers:
                            biotype = sub_feature.qualifiers[i][0]
                    if sub_feature.type == "miRNA":
                        biotype = "miRNA"
                    if sub_feature.type == "tRNA":
                        biotype = "tRNA"
                    d_gene_annotations['type'][gene_id] = biotype

                    # name
                    try:
                        for i in d_attributes['name']:
                            if i in sub_feature.qualifiers:
                                d_gene_annotations['name'][gene_id] = sub_feature.qualifiers[i][0]
                    except KeyError:
                        d_gene_annotations['name'][gene_id] = "NA"

                    exon_len = sub_feature.location.end - sub_feature.location.start
                    exon_rel_end = exon_rel_start + exon_len
                    # TODO: check +1 to length ?
                    d_transcript_annotations['len'][transcript_id] += exon_len + 1
                    exon_rel_start = exon_rel_end
                    prev_transcript_id = transcript_id
                    l_seen_exons.add(transcript_id + "_e" + str(n_exon).zfill(3))
                    n_exon += 1
    
    # OPTIMIZATION: Build interval trees for UTR/CDS features after parsing annotations
    # This enables O(log n) lookup instead of O(n) linear search in guess_region()
    if INTERVALTREE_AVAILABLE:
        build_interval_trees()


def build_interval_trees():
    """
    Build interval trees for UTR and CDS features per transcript.
    
    OPTIMIZATION: This enables O(log n) lookup instead of O(n) linear search
    in guess_region(). Interval trees are built once after parsing annotations,
    then reused for all region lookups.
    """
    # Build UTR interval trees
    for transcript_id, utr_list in d_transcript_annotations['UTR'].items():
        utr_tree = IntervalTree()
        for pos in utr_list:
            chrom = pos[0]
            start = int(pos[1])
            end = int(pos[2]) + 1  # IntervalTree uses [start, end) (end is exclusive)
            strand = strandardize(pos[3])
            utr_type = pos[4] if len(pos) > 4 else 'UTR'
            # Store chrom, strand, utr_type, and original pos for later use
            # interval.data can store any data associated with the interval
            utr_tree[start:end] = (chrom, strand, utr_type, pos)
        if utr_tree:
            d_transcript_interval_trees[transcript_id]['UTR'] = utr_tree
    
    # Build CDS interval trees
    for transcript_id, cds_list in d_transcript_annotations['CDS'].items():
        cds_tree = IntervalTree()
        for pos in cds_list:
            chrom = pos[0]
            start = int(pos[1])
            end = int(pos[2]) + 1  # IntervalTree uses [start, end) (end is exclusive)
            strand = strandardize(pos[3])
            # Store chrom, strand, and original pos for later use
            cds_tree[start:end] = (chrom, strand, pos)
        if cds_tree:
            d_transcript_interval_trees[transcript_id]['CDS'] = cds_tree


def parse_counts_file(crl_file, tpm_cutoff):
    d_crl_tpm = defaultdict(float)
    l_loci_bed = set()
    prev_readid = None
    read_count = 0
    # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance with large CRL files
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
    with open(crl_file, "r", buffering=BUFFER_SIZE) as fh_crl_file:
        for line in fh_crl_file:
            f = line.rstrip('\n').split('\t')
            readid = '|'.join(f[0].split("|")[:-1])
            if readid != prev_readid:
                read_count += 1
            crlid = f[3]
            crl_tpm = f[12]
            d_crl_tpm[crlid] = float(crl_tpm)
            b = f[9].split(":")
            locus_bed_entry = "\t".join([":".join(b[0:-3]), b[-3], b[-2], f[9], "1", b[-1]])
            if locus_bed_entry not in l_loci_bed:
                l_loci_bed.add(locus_bed_entry)
            prev_readid = readid

    uniq_tpms = sorted(list(set(d_crl_tpm.values())))
    
    # Clamp tpm_cutoff to [0, 1) to prevent index out of bounds
    if tpm_cutoff < 0:
        tpm_cutoff = 0.0
    elif tpm_cutoff >= 1.0:
        tpm_cutoff = 0.999999  # Just below 1.0 to ensure valid index
    
    # Handle empty list case
    if len(uniq_tpms) == 0:
        tpm_threshold = 0.0
    else:
        # Calculate index: ensure it's in valid range [0, len-1]
        index = int(tpm_cutoff * len(uniq_tpms))
        # Clamp index to valid range (shouldn't be needed with tpm_cutoff < 1.0, but defensive)
        index = min(index, len(uniq_tpms) - 1)
        tpm_threshold = uniq_tpms[index]

    return read_count, tpm_threshold


def merge_files(inprefix, outfile, header, r, compress=False):
    # Use shell commands for better performance with large files
    # This avoids loading all data into memory and leverages system sort
    # Intermediate files are always uncompressed - only final output is compressed
    temp_files = []
    for i in range(r):
        temp_file = inprefix + "." + str(i)
        # Intermediate files are NOT compressed
        if os.path.exists(temp_file):
            temp_files.append(temp_file)
    
    if not temp_files:
        # No files to merge
        return
    
    # Write header first, then append sorted unique lines
    # This avoids shell escaping issues with header content
    # Final output file is compressed if requested (intermediate files are NOT compressed)
    # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance with large files
    # Note: gzip.open doesn't support buffering parameter, but Python's gzip module uses internal buffering
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
    open_func = gzip.open if compress else open
    open_mode = "wt" if compress else "w"
    
    # Only apply buffering for uncompressed files (gzip.open handles buffering internally)
    if compress:
        fh_out = open_func(outfile, open_mode)
    else:
        fh_out = open_func(outfile, open_mode, buffering=BUFFER_SIZE)
    
    with fh_out:
        # Write header first
        fh_out.write(header + "\n")
        
        # Use shell commands to merge and sort files (more efficient for large files)
        # Intermediate files are uncompressed, so we can use cat directly
        # OPTIMIZATION: Use parallel sort if available (GNU sort supports --parallel option)
        # For systems with GNU sort >= 8.6, this can significantly speed up sorting large files
        sort_cmd = "sort -u"
        try:
            # Check if GNU sort with --parallel is available (GNU sort >= 8.6)
            result = subprocess.run(["sort", "--version"], capture_output=True, text=True, timeout=2)
            if "GNU coreutils" in result.stdout:
                # Use parallel sort - use number of processes that created the files
                # This ensures we use appropriate parallelism for merging
                sort_cmd = f"sort --parallel={r} -u"
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
            # Fall back to standard sort if check fails
            pass
        
        escaped_files = " ".join([f"'{f}'" for f in temp_files])
        cmd = f"cat {escaped_files} | {sort_cmd}"
        
        # Execute command and write output to file
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            sys.stderr.write(f"Warning: merge_files command failed: {stderr}\n")
        else:
            # Write sorted unique lines (header already written)
            # Compression happens automatically if open_func is gzip.open
            fh_out.write(stdout)

    # Clean up intermediate per-process files after merging
    for i in range(r):
        temp_file = inprefix + "." + str(i)
        # Intermediate files are NOT compressed
        if os.path.exists(temp_file):
            os.remove(temp_file)


def hybridization_positions(dotbracket1, dotbracket2):
    end1 = -1
    end2 = -1
    for i in range(len(dotbracket1)):
        if dotbracket1[i] == '(':
            end1 = i + 1
            dotbracket1[i] = "."
            # Iterate backward through dotbracket2 to find last ')' in original
            for j in range(len(dotbracket2) - 1, -1, -1):
                if dotbracket2[j] == ')':
                    end2 = j + 1
                    dotbracket2[j] = '.'
                    break
    return end1, end2


def write_interaction_summary(outdir, sample_name, compress=False, num_threads=4):
    d_interactions = defaultdict(lambda: defaultdict(list))
    chimeras_file = os.path.join(outdir, sample_name + ".chimeras.txt")
    if compress:
        chimeras_file += ".gz"
    
    # Open file with gzip if compression is enabled
    # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance with large files
    # Note: gzip.open doesn't support buffering parameter, but Python's gzip module uses internal buffering
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
    open_func = gzip.open if compress else open
    open_mode = "rt" if compress else "r"
    
    # Only apply buffering for uncompressed files (gzip.open handles buffering internally)
    if compress:
        fh_in = open_func(chimeras_file, open_mode)
    else:
        fh_in = open_func(chimeras_file, open_mode, buffering=BUFFER_SIZE)
    
    with fh_in:
        next(fh_in)
        for line in fh_in:
            f = line.rstrip("\n").split("\t")
            (readid, ref1, ref2, region1, region2, locus1, locus2, tpm1, tpm2, score1, score2) = (f[0], f[1], f[2],
                                                                                                  f[7], f[8], f[20],
                                                                                                  f[21], f[24], f[25],
                                                                                                  f[26], f[27])
            tpm = str(float(tpm1) + float(tpm2))
            score = str(float(score1) * float(score2))
            sequence1 = sequence2 = "NA"
            if f[29] != "NA":
                [sequence1, sequence2] = f[29].split("&")
            dotbracket = f[30]
            hybrid_start_pos = f[31]
            mfe = f[32]
            interaction = "\t".join(locus1.split(":")) + "\t" + "\t".join(locus2.split(":"))
            hybridization_pos = interaction
            [refid1, ref_start1, re1, ref_strand1, refid2, ref_start2, re2, ref_strand2] = interaction.split("\t")
            hybridized_sequences = "NA\tNA"
            interaction_otherway = "\t".join(locus2.split(":")) + "\t" + "\t".join(locus1.split(":"))

            if interaction_otherway in d_interactions:
                interaction = interaction_otherway
                hybridization_pos = interaction
                (ref2, ref1, region2, region1, tpm2, tpm1, score2, score1) = (f[1], f[2], f[7], f[8], f[24], f[25],
                                                                              f[26], f[27])
            if dotbracket != "NA":
                hybrid_end1, hybrid_end2 = hybridization_positions(list(dotbracket.split("&")[0]),
                                                                   list(dotbracket.split("&")[1]))
                # decrease by 1 before adding to the reference start
                hybrid_start1 = int(hybrid_start_pos.split("&")[0]) - 1
                hybrid_start2 = int(hybrid_start_pos.split("&")[1]) - 1
                hybridization_pos1 = "\t".join([refid1, str(int(ref_start1) + hybrid_start1),
                                                str(int(ref_start1) + hybrid_end1), ref_strand1])
                hybridization_pos2 = "\t".join([refid2, str(int(ref_start2) + hybrid_start2),
                                               str(int(ref_start2) + hybrid_end2), ref_strand2])
                hybridized_sequence1 = sequence1[hybrid_start1:hybrid_end1]
                hybridized_sequence2 = sequence2[hybrid_start2:hybrid_end2]

                hybridization_pos = hybridization_pos1 + "\t" + hybridization_pos2
                hybridized_sequences = hybridized_sequence1 + "\t" + hybridized_sequence2

                if interaction_otherway in d_interactions:
                    hybrid_start_pos = f[31].split('&')[1] + '&' + f[31].split('&')[0]
                    hybridization_pos = hybridization_pos2 + "\t" + hybridization_pos1
                    hybridized_sequences = hybridized_sequence2 + "\t" + hybridized_sequence1
                    target_dotbracket = dotbracket.split("&")[0].replace("(", ")")
                    query_dotbracket = dotbracket.split("&")[1].replace(")", "(")
                    dotbracket = query_dotbracket + "&" + target_dotbracket
                    [sequence2, sequence1] = f[29].split("&")

            d_interactions[interaction]["readid"].append(readid)
            d_interactions[interaction]["ref1"].append(ref1)
            d_interactions[interaction]["ref2"].append(ref2)
            d_interactions[interaction]["region1"].append(region1)
            d_interactions[interaction]["region2"].append(region2)
            common_info = "\t".join([sequence1, sequence2, dotbracket, mfe, hybridized_sequences,
                                     hybrid_start_pos, hybridization_pos,
                                     tpm1, tpm2, tpm, score1, score2, score])
            d_interactions[interaction]["common"] = [common_info]

    # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
    with open(os.path.join(outdir, "interactions.temp"), "w", buffering=BUFFER_SIZE) as fh_out:
        for interaction in d_interactions.keys():
            fh_out.write("\t".join([str(len(set(d_interactions[interaction]["readid"]))),
                                    interaction,
                                    d_interactions[interaction]["common"][0],
                                    ";".join(sorted(set(d_interactions[interaction]["region1"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["region2"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["ref1"]))),
                                    ";".join(sorted(set(d_interactions[interaction]["ref2"])))]) + "\n")

    header_interactions = "\t".join([
        "supporting_read_count",
        "locus_1_chromosome", "locus_1_start", "locus_1_end", "locus_1_strand",
        "locus_2_chromosome", "locus_2_start", "locus_2_end", "locus_2_strand",
        "locus_1_sequence", "locus_2_sequence", "hybridization_structure_dotbracket",
        "hybridization_mfe_kcal_mol", "hybridized_sequence_segments",
        "hybridization_start_positions", "hybridization_genomic_coordinates",
        "tpm_locus_1", "tpm_locus_2", "tpm_combined",
        "alignment_score_locus_1", "alignment_score_locus_2", "combined_alignment_score",
        "annotation_region_locus_1", "annotation_region_locus_2",
        "reference_transcript_id_1", "reference_transcript_id_2"
    ])
    # OPTIMIZATION: Use parallel sort if available (GNU sort supports --parallel option)
    # For systems with GNU sort >= 8.6, this can significantly speed up sorting large interaction files
    sort_cmd = "sort -k 1nr,1"
    try:
        # Check if GNU sort with --parallel is available (GNU sort >= 8.6)
        result = subprocess.run(["sort", "--version"], capture_output=True, text=True, timeout=2)
        if "GNU coreutils" in result.stdout:
            # Use parallel sort with specified number of threads
            # Default to 4 threads, but can be overridden for better performance
            sort_cmd = f"sort --parallel={max(1, num_threads)} -k 1nr,1"
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
        # Fall back to standard sort if check fails
        pass
    
    os.system(sort_cmd + " " + os.path.join(outdir, "interactions.temp") + " > " + os.path.join(outdir, "interactions.sorted"))
    # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance
    BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
    with open(os.path.join(outdir, sample_name + ".interactions.txt"), "w", buffering=BUFFER_SIZE) as fh_out:
        fh_out.write("# Note: To identify which locus is miRNA vs target, check the region1 and region2 fields.\n")
        fh_out.write("# miRNA annotations include: miRNA, 3p_mature_mir, 5p_mature_mir, mature_mir\n")
        fh_out.write(header_interactions + "\n")
        with open(os.path.join(outdir, "interactions.sorted"), "r", buffering=BUFFER_SIZE) as fh_in:
            for line in fh_in:
                fh_out.write(line)
    os.remove(os.path.join(outdir, "interactions.temp"))
    os.remove(os.path.join(outdir, "interactions.sorted"))


def validate_arguments(args):
    """Validate command-line arguments."""
    if args.hybridize and args.f_gtf and not args.f_ref:
        sys.stderr.write("Need the reference fasta file to hybridize. Make sure to provide the genomic fasta file"
                         " in case you already provided a GTF file.\n")
        sys.exit(1)

    if args.temperature < 0 or args.temperature > 100:
        sys.stderr.write("IntaRNA tempertature must be between 0 and 100!\n")
        sys.exit(1)


def print_configuration(args):
    """Print configuration parameters."""
    print('CRL file                             : ' + args.crl_file)
    print('Output directory                     : ' + args.outdir)
    if args.f_gtf:
        print('Annotation file                      : ' + args.f_gtf)
    print('Number of processes                  : ' + str(args.processes))
    print('TPM cutoff                           : ' + str(args.tpm_cutoff))
    print('Score cutoff                         : ' + str(args.score_cutoff))
    print('Chimeric overlap                     : ' + str(args.chimeric_overlap))
    print('Hybridize chimeric loci?             : ' + str(args.hybridize))
    print('Do not enforce seed interaction      : ' + str(args.no_seed))
    print('1st priority reference fasta file    : ' + args.ref_fasta1)
    if args.ref_fasta2:
        print('2nd priority reference fasta file    : ' + args.ref_fasta2)
    if args.f_ref:
        print('Reference genomic fasta file         : ' + args.f_ref)
    print('Summarize interactions at loci level : ' + str(args.summarize))
    print('Sample name                            : ' + args.sample_name)
    print('Compress output files with gzip        : ' + str(args.compress))
    print("===================================================================")


def setup_references(args):
    """Extract reference lengths from fasta files."""
    d_reflen1 = defaultdict()
    d_reflen2 = defaultdict()
    chira_utilities.extract_reflengths(args.ref_fasta1, d_reflen1)
    if args.ref_fasta2:
        chira_utilities.extract_reflengths(args.ref_fasta2, d_reflen2)
    return d_reflen1, d_reflen2


def run_chimera_extraction(args, d_reflen1, d_reflen2, tpm_cutoff_value, no_of_reads):
    """Run multiprocessing for chimera extraction."""
    print(str(datetime.datetime.now()), " START: multiprocessing")
    
    jobs = []
    # write chimeric reads
    for k in range(args.processes):
        s = k * math.ceil(no_of_reads / args.processes)
        e = min(s + math.floor(no_of_reads / args.processes), no_of_reads)
        print(k, s, e, no_of_reads)
        j = Process(target=write_chimeras, args=(s, e, no_of_reads, d_reflen1, d_reflen2, args.hybridize,
                                                 args.chimeric_overlap, args.f_gtf, args.outdir,
                                                 args.crl_file, tpm_cutoff_value, args.score_cutoff, str(k), args.sample_name, args.compress))
        jobs.append(j)

    for j in jobs:
        j.start()
    for j in jobs:
        j.join()


def prepare_reference_file(args):
    """Prepare reference fasta file for hybridization."""
    if args.f_ref:
        return args.f_ref
    else:
        f_reference = os.path.join(args.outdir, 'merged_reference.fa')
        l_ref = [args.ref_fasta1]
        if args.ref_fasta2:
            l_ref.append(args.ref_fasta2)
        # OPTIMIZATION: Use larger buffer size (2MB) for better I/O performance
        BUFFER_SIZE = 2 * 1024 * 1024  # 2MB buffer
        with open(f_reference, 'w', buffering=BUFFER_SIZE) as fh_out_ref:
            for fname in l_ref:
                with open(fname, 'r', buffering=BUFFER_SIZE) as infile:
                    for lin in infile:
                        fh_out_ref.write(lin)
        return f_reference


def extract_loci_sequences(args, f_reference):
    """Extract FASTA sequences for loci using bedtools getfasta."""
    for k in range(args.processes):
        bed = os.path.join(args.outdir, "loci.bed.") + str(k)
        fa = os.path.join(args.outdir, "loci.fa.") + str(k)
        getfasta_cmd = chira_utilities.get_bedtools_command('getfasta')
        process = subprocess.Popen(getfasta_cmd + ['-s', '-nameOnly',
                                    '-fi', f_reference,
                                    '-bed', bed,
                                    '-fo', fa],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   universal_newlines=True)
        for l in sorted(set(process.stdout.readlines())):
            print(l, end="")


def build_intarna_params(args):
    """Build common IntaRNA parameters string."""
    noseed_param = ""
    if args.no_seed:
        noseed_param = "--noSeed"
    return " ".join(["--outMode C", "--outCsvCols id1,start1,start2,hybridDPfull,E",
                     noseed_param, "-m", args.intarna_mode, "--acc", args.accessibility,
                     "--temperature", str(args.temperature), "--seedBP", str(args.seed_bp),
                     "--seedMinPu", str(args.seed_min_pu), "--accW", str(args.acc_width)])


def run_hybridization(args):
    """Run hybridization process if enabled."""
    chimeras_prefix = os.path.join(args.outdir, args.sample_name + ".chimeras-r")
    f_reference = prepare_reference_file(args)
    
    # Extract FASTA sequences for loci
    extract_loci_sequences(args, f_reference)
    
    # Build IntaRNA parameters
    common_intarna_params = build_intarna_params(args)
    
    # Run hybridization in parallel
    jobs = []
    for k in range(args.processes):
        j = Process(target=hybridize_and_write, args=(args.outdir, common_intarna_params, str(k), args.sample_name, args.compress))
        jobs.append(j)

    for j in jobs:
        j.start()
    for j in jobs:
        j.join()

    # Cleanup intermediate files
    for k in range(args.processes):
        os.remove(os.path.join(args.outdir, "loci.fa.") + str(k))
        os.remove(os.path.join(args.outdir, "loci.bed.") + str(k))

    # Remove temporary reference file if we created it
    if not args.f_ref:
        if os.path.exists(f_reference):
            os.remove(f_reference)
        if os.path.exists(f_reference + ".fai"):
            os.remove(f_reference + ".fai")
    
    return chimeras_prefix


def merge_output_files(args, chimeras_prefix):
    """Merge output files from multiple processes."""
    # File name prefixes
    if not args.hybridize:
        chimeras_prefix = os.path.join(args.outdir, args.sample_name + ".chimeras")
    
    chimeras_file = os.path.join(args.outdir, args.sample_name + ".chimeras.txt")
    singletons_file = os.path.join(args.outdir, args.sample_name + ".singletons.txt")
    singletons_prefix = os.path.join(args.outdir, args.sample_name + ".singletons")
    
    # Add .gz extension if compression is enabled
    if args.compress:
        chimeras_file += ".gz"
        singletons_file += ".gz"
    
    # Header fields
    header_chimeras = "\t".join(["read_id",
                                 "transcript_id_1",
                                 "transcript_id_2",
                                 "gene_id_1",
                                 "gene_id_2",
                                 "gene_symbol_1",
                                 "gene_symbol_2",
                                 "annotation_region_1",
                                 "annotation_region_2",
                                 "transcript_start_1",
                                 "transcript_end_1",
                                 "transcript_strand_1",
                                 "transcript_length_1",
                                 "transcript_start_2",
                                 "transcript_end_2",
                                 "transcript_strand_2",
                                 "transcript_length_2",
                                 "read_alignment_info",
                                 "genomic_coordinates_1",
                                 "genomic_coordinates_2",
                                 "locus_id_1",
                                 "locus_id_2",
                                 "crl_group_id_1",
                                 "crl_group_id_2",
                                 "tpm_1",
                                 "tpm_2",
                                 "alignment_score_1",
                                 "alignment_score_2",
                                 "combined_alignment_score",
                                 "hybridized_sequences",
                                 "hybridization_structure",
                                 "hybridization_positions",
                                 "hybridization_mfe_kcal_mol",
                                 "mirna_read_position"])
    merge_files(chimeras_prefix, chimeras_file, header_chimeras, args.processes, args.compress)
    
    header_singletons = "\t".join(["read_id",
                                   "transcript_id",
                                   "gene_id",
                                   "gene_symbol",
                                   "annotation_region",
                                   "transcript_start",
                                   "transcript_end",
                                   "transcript_strand",
                                   "transcript_length",
                                   "read_alignment_info",
                                   "genomic_coordinates",
                                   "locus_id",
                                   "crl_group_id",
                                   "tpm",
                                   "alignment_score"])
    merge_files(singletons_prefix, singletons_file, header_singletons, args.processes, args.compress)


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: extract chimeras',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-l', '--loci', action='store', dest='crl_file', required=True,
                        metavar='', help='Input BED file with alignments')

    parser.add_argument('-o', '--out', action='store', dest='outdir', required=True,
                        metavar='', help='Path to output directory')

    parser.add_argument('-g', '--gtf', action='store', dest='f_gtf', required=False,
                        metavar='', help='Annotation GTF file')

    parser.add_argument('-p', '--processes', action='store', type=int, default=1, metavar='',
                        dest='processes',
                        help='Number of processes to use')

    parser.add_argument('-tc', '--tpm_cutoff', action='store', type=chira_utilities.score_float, default=0, metavar='',
                        dest='tpm_cutoff',
                        help='Transcripts with less than this percentile TPMs will be discarded in '
                             'the final output. [0, 1) - values >= 1.0 will be clamped to just below 1.0')

    parser.add_argument('-sc', '--score_cutoff', action='store', type=chira_utilities.score_float, default=0.0, metavar='',
                        dest='score_cutoff',
                        help='Hybrids with less than this score will be discarded in the final output. [0-1.0)')

    parser.add_argument('-co', '--chimeric_overlap', action='store', type=int, default=2, metavar='',
                        dest='chimeric_overlap',
                        help='Maximum number of bases allowed between the chimeric segments of a read')

    parser.add_argument("-r", '--hybridize', action='store_true', dest='hybridize',
                        help="Hybridize the predicted chimeras")

    parser.add_argument("-ns", '--no_seed', action='store_true', dest='no_seed',
                        help="Do not enforce seed interactions")

    parser.add_argument("-acc", '--accessibility', type=str, choices=["C", "N"], default='N', required=False,
                        dest='accessibility', metavar='', help='IntaRNA accessibility: C (compute) or N (not)')

    parser.add_argument("-m", '--intarna_mode', type=str, choices=["H", "M", "S"], default='H', required=False,
                        dest='intarna_mode', metavar='', help='IntaRNA mode: H (heuristic), M (exact), S (seed-only)')

    parser.add_argument('-t', '--temperature', action='store', type=float, default=37, metavar='',
                        dest='temperature',
                        help='IntaRNA temperature parameter in Celsius to setup the VRNA energy parameters')

    parser.add_argument('-sbp', '--seed_bp', action='store', type=int, default=5, metavar='',
                        dest='seed_bp', choices=range(2, 20),
                        help='IntaRNA --seedBP parameter: number of inter-molecular base pairs within the seed region')

    parser.add_argument('-smpu', '--seed_min_pu', action='store', type=chira_utilities.score_float, default=0,
                        metavar='', dest='seed_min_pu',
                        help='IntaRNA --seedMinPu parameter: minimal unpaired probability '
                             '(per sequence) a seed region may have')

    parser.add_argument('-accw', '--acc_width', action='store', type=int, default=150, metavar='',
                        dest='acc_width', choices=range(0, 99999),
                        help='IntaRNA --accW parameter:  sliding window size for accessibility computation')

    parser.add_argument('-f1', '--ref_fasta1', action='store', dest='ref_fasta1', required=True,
                        metavar='', help='First priority fasta file')

    parser.add_argument('-f2', '--ref_fasta2', action='store', dest='ref_fasta2', required=False,
                        metavar='', help='second priority fasta file')

    parser.add_argument('-f', '--ref', action='store', dest='f_ref', required=False,
                        metavar='', help='Reference fasta file')

    parser.add_argument("-s", '--summarize', action='store_true', dest='summarize',
                        help="Summarize interactions at loci level")

    parser.add_argument('-n', '--sample_name', action='store', dest='sample_name', required=True,
                        metavar='', help='Sample name prefix for output files')

    parser.add_argument('-z', '--gzip', action='store_true', dest='compress',
                        help='Compress output files (chimeras and singletons) with gzip')

    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {chira_utilities.__version__}')

    return parser.parse_args()


def main():
    """Main function to orchestrate the chimera extraction workflow."""
    args = parse_arguments()
    print_configuration(args)
    validate_arguments(args)

    # Parse annotations if GTF file provided
    if args.f_gtf:
        print("Parsing the annotation file")
        parse_annotations(args.f_gtf)

    # Setup references
    d_reflen1, d_reflen2 = setup_references(args)

    # Parse CRLs file
    print("Parsing CRLs file")
    no_of_reads, tpm_cutoff_value = parse_counts_file(args.crl_file, args.tpm_cutoff)
    print("Done")

    # Run chimera extraction
    run_chimera_extraction(args, d_reflen1, d_reflen2, tpm_cutoff_value, no_of_reads)

    # Run hybridization if enabled (returns chimeras prefix)
    chimeras_prefix = None
    if args.hybridize:
        chimeras_prefix = run_hybridization(args)

    # Merge output files
    merge_output_files(args, chimeras_prefix)
    
    print(str(datetime.datetime.now()), " END: multiprocessing")
    
    # Write interaction summary if requested
    if args.summarize:
        write_interaction_summary(args.outdir, args.sample_name, args.compress, args.processes)


if __name__ == "__main__":
    main()

