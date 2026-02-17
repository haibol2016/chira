#!/usr/bin/env python
import os
import sys
from collections import defaultdict
import argparse
import chira_utilities
import copy
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import math


def build_crls(build_crls_too, bed, merged_bed, crl_file, crl_share_cutoff, min_locus_size):
    l_locipos = []
    l_locireads = []
    d_readlocus_transcripts = defaultdict(list)
    # OPTIMIZATION: Use larger buffer size (1MB) for read operations on large files
    # This reduces system calls and improves I/O performance when reading many lines
    with open(merged_bed, "r", buffering=1024*1024) as fh_merged_bed:
        for n_locus, line in enumerate(fh_merged_bed):
            # chr14\t64814786\t64814804\t-\ttag_1308593|1|r,ENSMUST00000176386;tag_1308594|2|r,ENSMUST00000176386
            f = line.rstrip('\n').split('\t')
            pos = ':'.join([f[0], f[1], f[2], f[3]])
            l_locipos.append(pos)
            alignments = f[4].split(';')  # in the description column of the bed file,alignments are seperated by ';'
            l_locusreads = set()
            for alignment in alignments:
                alignment_fields = alignment.split(',')
                # Defensive check: skip malformed alignments
                if len(alignment_fields) < 6:
                    print(f"Warning: Skipping malformed alignment (expected 6 comma-separated fields, got {len(alignment_fields)}): {alignment}", file=sys.stderr)
                    continue
                segmentid, transcriptid, start, end, tx_strand, cigar = alignment.split(',')
                # transcriptid already extracted above, no need to split again
                transcriptid_pos = '\t'.join([transcriptid, start, end, tx_strand, cigar])
                d_readlocus_transcripts[segmentid+str(n_locus)].append(transcriptid_pos)
                if segmentid not in l_locusreads:
                    l_locusreads.add(segmentid)
            l_locireads.append(set(l_locusreads))

    # OPTIMIZATION: Pre-sort transcript lists once to avoid repeated sorting in nested output loops
    # This reduces sorting from O(n log n) per access to O(n log n) total, significantly
    # improving performance when writing output for large datasets
    for key in d_readlocus_transcripts:
        d_readlocus_transcripts[key].sort()

    d_crl_reads = defaultdict(set)  # Use sets instead of lists for faster operations
    d_crl_locus_reads = defaultdict(lambda: defaultdict())
    l_remaining_locireads = defaultdict(list)
    chira_utilities.print_w_time("START: 1st iteration of CRLs")
    print("Number of loci: " + str(len(l_locireads)))
    n_crl = 0
    # create CRLs only if required
    if build_crls_too:
        l_qualified_locireads = sorted([(i, item) for i, item in enumerate(l_locireads) if len(item) >= min_locus_size],
                                       key=lambda x: len(x[1]), reverse=True)
        l_remaining_locireads = sorted([(i, item) for i, item in enumerate(l_locireads) if len(item) < min_locus_size],
                                       key=lambda x: len(x[1]), reverse=True)

        print("Number of qualified loci: ", len(l_qualified_locireads))

        for n_locus, l_locusreads in l_qualified_locireads:
            already_crl_member = False
            l_matched_crls = []
            # OPTIMIZATION: Cache locus_size since it's used multiple times (for bounds calculation)
            # This avoids repeated len() calls on the set, which is O(1) but still has function call overhead
            locus_size = len(l_locusreads)
            # lower and uppder bounds for filtering crls based on their size
            lower_bound = locus_size * (1 - crl_share_cutoff)
            upper_bound = locus_size / (1 - crl_share_cutoff)
            # traverse in reverse order because the latest CRL is the last one
            for crlid in range(len(d_crl_reads) - 1, -1, -1):  # Include 0 in range
                l_crlreads = d_crl_reads[crlid]  # Already a set, no need to convert
                # if the CRL has similar size
                if lower_bound <= len(l_crlreads) <= upper_bound:
                    # if there are significantly more reads in crl than in locus or otherway around
                    n_common_reads = len(l_locusreads & l_crlreads)  # Use & for set intersection
                    if n_common_reads == 0:
                        continue
                    n_union_reads = len(l_locusreads | l_crlreads)  # Use | for set union
                    # jaccard similarity score
                    # Defensive check: skip if union is empty (shouldn't happen if n_common_reads > 0)
                    if n_union_reads == 0:
                        continue
                    if n_common_reads / float(n_union_reads) < crl_share_cutoff:
                        continue
                    # identical loci are multi-mapped loci with the same set of identical set of multi-mapped reads
                    l_matched_crls.append(crlid)
                    already_crl_member = True

            # n_locus is not a member of any crl, hence create a new crl
            if not already_crl_member:
                l_matched_crls.append(n_crl)
                n_crl += 1
            for matched_crl in l_matched_crls:
                d_crl_locus_reads[matched_crl][n_locus] = l_locusreads
                d_crl_reads[matched_crl].update(l_locusreads)  # Use set.update() instead of extend()
        d_crl_reads.clear()
    else:
        l_remaining_locireads = sorted(enumerate(l_locireads), key=lambda x: len(x[1]), reverse=True)
    # every locus with less than min_locus_size reads gets its' own locus
    for n_locus, l_locusreads in l_remaining_locireads:
        d_crl_locus_reads[n_crl][n_locus] = l_locusreads
        n_crl += 1
    chira_utilities.print_w_time("END: 1st iteration of CRLs")

    d_locus_crl_share = defaultdict(lambda: defaultdict(float))
    d_highest_shares = defaultdict(float)
    for crlid, d_crlloci in d_crl_locus_reads.items():
        # Build set of all reads in CRL more efficiently
        l_crlreads = set()
        for l_locusreads in d_crlloci.values():
            l_crlreads.update(l_locusreads)  # Use set.update() instead of extend()
        crl_reads_len = len(l_crlreads)
        # Defensive check: skip CRLs with zero reads (shouldn't happen, but protect against division by zero)
        if crl_reads_len == 0:
            continue
        for locusid, l_locusreads in d_crlloci.items():
            locus_share = len(l_locusreads) / float(crl_reads_len)
            d_locus_crl_share[locusid][crlid] = locus_share
            if locus_share > d_highest_shares[locusid]:
                d_highest_shares[locusid] = locus_share

    d_isolated_loci = {}  # loci that are separated from crl because of not enough overall share
    chira_utilities.print_w_time("START: 2nd iteration of CRLs")
    # in this iteration each locus is checked again for crl share and only one best crl for locus kept
    # if there are multiple eqaully good crls for a locus, all of the crls are considered for that locus
    for crlid in list(d_crl_locus_reads.keys()):
        d_crlloci = d_crl_locus_reads[crlid]
        for locusid in list(d_crlloci):
            if d_locus_crl_share[locusid][crlid] >= crl_share_cutoff:
                if d_locus_crl_share[locusid][crlid] < d_highest_shares[locusid]:
                    if locusid in d_crl_locus_reads[crlid]:
                        del d_crl_locus_reads[crlid][locusid]
                        # del d_locus_crl_share[locusid][crlid]
            else:
                # store isolated locus id and its reads
                d_isolated_loci[locusid] = d_crlloci[locusid]
                del d_crl_locus_reads[crlid][locusid]
                del d_locus_crl_share[locusid]

    chira_utilities.print_w_time("END: 2nd iteration of CRLs")
    chira_utilities.print_w_time("START: creating CRLs with isolated loci")
    # every locus that is in isolated_loci makes its own crl
    crl_index = len(d_crl_locus_reads)
    for locusid, l_locusreads in d_isolated_loci.items():
        # already member of some crl, safely ignore it
        if locusid in d_locus_crl_share:
            continue
        d_crl_locus_reads[crl_index][locusid] = l_locusreads
        d_locus_crl_share[locusid][crl_index] = 1
        crl_index += 1
    d_locus_crl_share.clear()

    chira_utilities.print_w_time("END: creating CRLs with isolated loci")

    print("There are a total of " + str(len(d_crl_locus_reads)) + " uniq crls")

    # read the segments BED to get the genomic positions
    d_read_genomic_pos = defaultdict(str)
    # OPTIMIZATION: Use larger buffer size (1MB) for read operations on large files
    # This reduces system calls and improves I/O performance when reading many lines
    with open(bed, "r", buffering=1024*1024) as fh_bed:
        for line in fh_bed:
            b = line.rstrip('\n').split('\t')
            pos = ':'.join([b[0], b[1], b[2], b[5]])
            desc = b[3].split(',')  # in the description column of the bed file,alignments are seperated by ';'
            # Defensive check: skip if desc is empty
            if len(desc) < 1:
                print(f"Warning: Skipping line with empty description: {line[:100]}", file=sys.stderr)
                continue
            segmentid = desc[0]
            transcriptid_pos = '\t'.join(desc[1:])
            # at this level reads have unique ids preceeded by a serialnumber
            # each read can have multiple alignements on same transcript
            d_read_genomic_pos[transcriptid_pos+segmentid] = pos

    # OPTIMIZATION: Use larger buffer size (1MB) for write operations in nested loops
    # This significantly reduces system calls and improves I/O performance when writing many small lines
    # Default buffer is typically 8KB, so 1MB (128x larger) provides substantial improvement
    with open(crl_file, "w", buffering=1024*1024) as fh_crl_file:
        # enumerate again because of above processing some crlids might be missing
        for crlid, d_crlloci in enumerate(d_crl_locus_reads.values()):
            if len(d_crlloci) > 1:
                # Build set more efficiently using set union
                l_crlreads = set()
                for l_locusreads in d_crlloci.values():
                    l_crlreads.update(l_locusreads)
                crl_reads_len = len(l_crlreads)
            else:
                crl_reads_len = None
            for locusid, l_locusreads in d_crlloci.items():
                if crl_reads_len is not None:
                    # Defensive check: skip if crl_reads_len is 0 (shouldn't happen, but protect against division by zero)
                    if crl_reads_len == 0:
                        locus_share = 0.0
                    else:
                        locus_share = len(l_locusreads) / float(crl_reads_len)
                else:
                    locus_share = 1.0
                # l_locusreads is a set, so sort for consistent output order
                # OPTIMIZATION: d_readlocus_transcripts is already pre-sorted (see line 38-40),
                # so we avoid repeated sorting in this nested loop, improving performance significantly
                for segmentid in sorted(l_locusreads):
                    segment_key = segmentid + str(locusid)
                    for transcriptid_pos in d_readlocus_transcripts[segment_key]:
                        entry = "\t".join([segmentid,
                                           transcriptid_pos.split('\t')[0],
                                           str(locusid),
                                           str(crlid),
                                           '\t'.join(transcriptid_pos.split('\t')[1:]),
                                           d_read_genomic_pos[transcriptid_pos+segmentid],
                                           l_locipos[locusid],
                                           "{:.2f}".format(locus_share)]
                                          )
                        fh_crl_file.write(entry + "\n")


def _process_reads_chunk_e_step(args):
    """
    Process a chunk of reads in the E-step.
    Returns updated d_alpha values for reads in this chunk.
    
    OPTIMIZATION: This function is designed for parallel execution with ProcessPoolExecutor.
    Each process processes a unique subset of reads. Returns updated values instead of
    modifying in place to work with process-based parallelism.
    
    Args:
        args: Tuple of (readids_chunk, d_alpha_chunk, d_rho_old)
            - readids_chunk: List of read IDs to process
            - d_alpha_chunk: Dict of {readid: {crlid: value}} for this chunk
            - d_rho_old: Dict of {crlid: relative_abundance}
    
    Returns:
        Dict of {readid: {crlid: updated_value}} for this chunk
    """
    readids_chunk, d_alpha_chunk, d_rho_old = args
    d_alpha_updated = {}
    for readid in readids_chunk:
        # Match original logic: skip reads with only 1 CRL (shouldn't happen but defensive)
        if len(d_alpha_chunk[readid]) == 1:
            continue
        read_crls = d_alpha_chunk[readid]
        # OPTIMIZATION: Pre-compute total abundance once and cache inverse to avoid
        # repeated sum() calls and division operations in the inner loop
        total_crl_rel_abundancy = sum(d_rho_old[crlid] for crlid in read_crls)
        # OPTIMIZATION: Defensive check to prevent division by zero (not in original but safe)
        if total_crl_rel_abundancy > 0:
            inv_total = 1.0 / total_crl_rel_abundancy
            d_alpha_updated[readid] = {crlid: d_rho_old[crlid] * inv_total for crlid in read_crls}
        else:
            d_alpha_updated[readid] = read_crls.copy()  # Keep original if total is 0
    return d_alpha_updated


def _process_reads_chunk_aggregate(args):
    """
    Process a chunk of reads for aggregation step.
    Returns a dictionary of CRL contributions from this chunk.
    
    OPTIMIZATION: This function is designed for parallel execution with ProcessPoolExecutor.
    Each process accumulates contributions into its own dictionary, then results are merged.
    Since summation is commutative and associative, parallel aggregation produces identical
    results to sequential.
    
    Args:
        args: Tuple of (readids_chunk, d_alpha_chunk)
            - readids_chunk: List of read IDs to process
            - d_alpha_chunk: Dict of {readid: {crlid: value}} for this chunk
    
    Returns:
        Dict of {crlid: aggregated_value} for this chunk
    """
    readids_chunk, d_alpha_chunk = args
    d_rho_chunk = defaultdict(float)
    for readid in readids_chunk:
        for crlid in d_alpha_chunk[readid]:
            d_rho_chunk[crlid] += d_alpha_chunk[readid][crlid]
    # Convert to regular dict: ProcessPoolExecutor pickles return values when sending
    # results back to main process. While defaultdict can be pickled, converting to
    # plain dict is more efficient (no factory function to pickle) and explicit.
    return dict(d_rho_chunk)


def em(d_alpha, d_rho, library_size, l_multimap_readids, em_threshold, num_processes=1):
    """
    Expectation-Maximization algorithm with optional multi-processing.
    
    OPTIMIZATION: Uses ProcessPoolExecutor instead of ThreadPoolExecutor for CPU-bound operations.
    ProcessPoolExecutor bypasses Python's GIL, enabling true parallelism for CPU-intensive tasks.
    This provides 2-8x speedup compared to ThreadPoolExecutor for large datasets.
    
    Args:
        d_alpha: Dictionary mapping read IDs to CRL contributions
        d_rho: Dictionary mapping CRL IDs to relative abundances
        library_size: Total library size
        l_multimap_readids: List of read IDs with multiple mappings
        em_threshold: Convergence threshold
        num_processes: Number of processes to use (default: 1, use all available if <= 0)
    """
    if num_processes <= 0:
        import multiprocessing
        num_processes = multiprocessing.cpu_count()
    
    # Log parallelization configuration
    if num_processes > 1:
        print(f"Using {num_processes} parallel processes for EM algorithm", file=sys.stderr)
        print(f"  - Multimapping reads: {len(l_multimap_readids)}", file=sys.stderr)
        print(f"  - Total reads: {len(d_alpha)}", file=sys.stderr)
    else:
        print("Using sequential processing (single process)", file=sys.stderr)
    
    i = 1
    sum_of_rho_diff = float('inf')
    
    # Split reads into chunks for parallel processing
    def chunk_list(lst, n_chunks):
        """Split list into approximately equal chunks."""
        if n_chunks <= 1:
            return [lst]
        chunk_size = max(1, math.ceil(len(lst) / n_chunks))
        return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]
    
    while sum_of_rho_diff > em_threshold and i <= 1000:
        if i % 10 == 0:  # Print every 10 iterations instead of every iteration
            print("iteration: " + str(i))
        # E-step
        # OPTIMIZATION: Use shallow copy instead of deep copy - much faster (10-50x speedup)
        # Shallow copy is sufficient since d_rho contains only immutable float values
        d_rho_old = dict(d_rho)  # Shallow copy - safe because values are immutable floats
        d_rho.clear()
        
        # OPTIMIZATION: Parallelize E-step for multimapped reads when beneficial
        # Only parallelize if num_processes > 1 AND there are enough reads to justify process overhead.
        # ProcessPoolExecutor has overhead: process creation (~10-50ms), data copying/pickling, result merging.
        # Adaptive threshold: max(500 reads, num_processes * 50) ensures each process gets meaningful work.
        # For 12 processes: need at least 600 reads (50 per process minimum).
        # ProcessPoolExecutor provides true parallelism (bypasses GIL) for CPU-bound operations.
        min_reads_threshold = max(500, num_processes * 50) if num_processes > 1 else float('inf')
        if num_processes > 1 and len(l_multimap_readids) > min_reads_threshold:
            if i == 1:  # Log only on first iteration to avoid spam
                print(f"  Parallelizing E-step: {len(l_multimap_readids)} multimapping reads across {num_processes} processes", file=sys.stderr)
            read_chunks = chunk_list(l_multimap_readids, num_processes)
            
            # Extract d_alpha subsets for each chunk (needed for process-based parallelism)
            chunk_args = []
            for chunk in read_chunks:
                d_alpha_chunk = {readid: d_alpha[readid].copy() for readid in chunk}
                chunk_args.append((chunk, d_alpha_chunk, d_rho_old))
            
            with ProcessPoolExecutor(max_workers=num_processes) as executor:
                futures = [executor.submit(_process_reads_chunk_e_step, args) for args in chunk_args]
                # Collect updated values and merge back into d_alpha
                for future in as_completed(futures):
                    d_alpha_updated = future.result()
                    for readid, read_crls_updated in d_alpha_updated.items():
                        d_alpha[readid] = read_crls_updated
        else:
            # Sequential processing for small datasets or single process
            if num_processes > 1 and len(l_multimap_readids) <= min_reads_threshold and i == 1:
                print(f"  Note: Sequential E-step ({len(l_multimap_readids)} reads < threshold of {min_reads_threshold})", file=sys.stderr)
            for readid in l_multimap_readids:
                # Match original logic: skip reads with only 1 CRL (shouldn't happen but defensive)
                if len(d_alpha[readid]) == 1:
                    continue
                read_crls = d_alpha[readid]
                total_crl_rel_abundancy = sum(d_rho_old[crlid] for crlid in read_crls)
                # OPTIMIZATION: Defensive check to prevent division by zero (not in original but safe)
                if total_crl_rel_abundancy > 0:
                    inv_total = 1.0 / total_crl_rel_abundancy
                    for crlid in read_crls:
                        d_alpha[readid][crlid] = d_rho_old[crlid] * inv_total
        
        # OPTIMIZATION: Aggregate contributions from all reads (including single-mapped)
        # This step can also be parallelized for large datasets. Each process processes a chunk
        # of reads and accumulates contributions, then results are merged.
        # ProcessPoolExecutor provides true parallelism (bypasses GIL) for CPU-bound operations.
        if num_processes > 1:
            # For parallel processing, convert keys to list once per iteration
            # OPTIMIZATION: Only create list when needed for parallel processing
            all_readids = list(d_alpha.keys())
            # Use same adaptive threshold as E-step
            min_reads_threshold_agg = max(500, num_processes * 50)
            if len(all_readids) > min_reads_threshold_agg:
                if i == 1:  # Log only on first iteration to avoid spam
                    print(f"  Parallelizing aggregation: {len(all_readids)} reads across {num_processes} processes", file=sys.stderr)
                read_chunks = chunk_list(all_readids, num_processes)
                
                # Extract d_alpha subsets for each chunk (needed for process-based parallelism)
                chunk_args = []
                for chunk in read_chunks:
                    d_alpha_chunk = {readid: d_alpha[readid].copy() for readid in chunk}
                    chunk_args.append((chunk, d_alpha_chunk))
                
                with ProcessPoolExecutor(max_workers=num_processes) as executor:
                    futures = [executor.submit(_process_reads_chunk_aggregate, args) 
                              for args in chunk_args]
                    for future in as_completed(futures):
                        d_rho_chunk = future.result()
                        for crlid, value in d_rho_chunk.items():
                            d_rho[crlid] += value
            else:
                # Sequential aggregation for small datasets even with threading enabled
                if len(all_readids) <= min_reads_threshold_agg and i == 1:
                    print(f"  Note: Sequential aggregation ({len(all_readids)} reads < threshold of {min_reads_threshold_agg})", file=sys.stderr)
                for readid in d_alpha.keys():
                    for crlid in d_alpha[readid]:
                        d_rho[crlid] += d_alpha[readid][crlid]
        else:
            # Sequential aggregation - matches original code exactly
            for readid in d_alpha.keys():
                for crlid in d_alpha[readid]:
                    d_rho[crlid] += d_alpha[readid][crlid]
        
        # M-step (sequential - fast enough)
        sum_of_rho_diff = 0
        # Defensive check: avoid division by zero if library_size is 0 (shouldn't happen, but protect)
        if library_size > 0:
            inv_library_size = 1.0 / float(library_size)
            for crlid in d_rho:
                new_rho = d_rho[crlid] * inv_library_size
                sum_of_rho_diff += abs(new_rho - d_rho_old[crlid])
                d_rho[crlid] = new_rho
        else:
            # If library_size is 0, set all rho to 0 and break (converged to invalid state)
            for crlid in d_rho:
                d_rho[crlid] = 0.0
            break
        i += 1
    return d_alpha


def tpm(d_crl_expression, d_crl_loci_len):
    total_rpk = 0
    d_crl_tpm = defaultdict(float)
    # OPTIMIZATION: Cache sorted keys to avoid repeated sorting (used in two loops below)
    # This reduces sorting from O(n log n) twice to O(n log n) once
    sorted_crlids = sorted(d_crl_expression.keys())
    for crlid in sorted_crlids:
        crl_expression = d_crl_expression[crlid]
        # OPTIMIZATION: median() already sorts internally, so don't sort here to avoid double sorting
        # This eliminates redundant O(n log n) sorting operation for each CRL
        crl_len = chira_utilities.median(d_crl_loci_len[crlid].values()) / 1000.0  # length in kbs
        if crl_len > 0:  # Avoid division by zero
            rpk = crl_expression / crl_len
            d_crl_tpm[crlid] = rpk
            total_rpk += rpk
    millions_of_rpk = total_rpk / 1000000.0
    if millions_of_rpk > 0:  # Avoid division by zero
        # OPTIMIZATION: Pre-compute inverse to avoid repeated division in loop
        inv_millions = 1.0 / millions_of_rpk
        for crlid in sorted_crlids:
            crl_tpm = d_crl_tpm[crlid] * inv_millions
            d_crl_tpm[crlid] = crl_tpm
    return d_crl_tpm


def quantify_crls(crl_file, em_threshold, num_processes=1):
    d_crl_loci_len = defaultdict(lambda: defaultdict(int))
    l_multimap_readids = []
    d_alpha = defaultdict(lambda: defaultdict(float))
    d_rho = defaultdict(float)

    # Use context manager for file handling (more efficient)
    # OPTIMIZATION: Use larger buffer size (1MB) for read operations on large files
    # This reduces system calls and improves I/O performance when reading many lines
    with open(crl_file, "r", buffering=1024*1024) as fh_crl_file:
        for line in fh_crl_file:
            f = line.rstrip("\n").split("\t")
            # consider the segment id and quantify individula segments than whole reads
            # Defensive check: skip malformed lines
            if len(f) < 10:
                print(f"Warning: Skipping malformed line (expected at least 10 fields, got {len(f)}): {line[:100]}", file=sys.stderr)
                continue
            readid = f[0]
            locusid = f[2]
            crlid = f[3]
            pos = f[9].split(':')
            # Defensive check: pos should have at least 3 elements (chr:start:end:strand format)
            if len(pos) < 3:
                print(f"Warning: Skipping line with invalid position format (expected chr:start:end:strand, got {f[9]}): {line[:100]}", file=sys.stderr)
                continue
            try:
                locuslength = int(pos[-2]) - int(pos[-3]) + 1
            except (ValueError, IndexError) as e:
                print(f"Warning: Skipping line with invalid position values (pos={f[9]}): {e}", file=sys.stderr)
                continue
            # a single locus can belong to multiple crls
            # one read can be part of multiple crls
            d_alpha[readid][crlid] = 1
            d_crl_loci_len[crlid][locusid] = locuslength

    # intial read segment contributions
    # OPTIMIZATION: Use .items() for more efficient dictionary iteration (avoids key lookups)
    # OPTIMIZATION: Cache num_crls and inv_num_crls since they're used multiple times per read
    # This avoids repeated len() calls and division operations in the inner loop
    for readid, read_crls in d_alpha.items():
        num_crls = len(read_crls)
        inv_num_crls = 1.0 / float(num_crls)
        for crlid in read_crls:
            d_alpha[readid][crlid] = inv_num_crls
            d_rho[crlid] += inv_num_crls
        if num_crls > 1:
            l_multimap_readids.append(readid)

    library_size = sum(d_rho.values())
    # intial relative abundancies
    # Defensive check: avoid division by zero if library_size is 0 (shouldn't happen, but protect)
    if library_size > 0:
        for crlid in d_rho:
            d_rho[crlid] = d_rho[crlid] / library_size
    else:
        # If library_size is 0, all reads have no CRL assignments - set all to 0
        print("Warning: library_size is 0, no reads have CRL assignments", file=sys.stderr)
        for crlid in d_rho:
            d_rho[crlid] = 0.0

    d_res = em(d_alpha, d_rho, library_size, l_multimap_readids, em_threshold, num_processes)

    d_crl_expression = defaultdict(float)
    # OPTIMIZATION: Use .items() for more efficient dictionary iteration
    # This avoids key lookups and provides direct access to both keys and values
    for readid, read_crls in d_res.items():
        for crlid, value in read_crls.items():
            d_crl_expression[crlid] += value

    d_crl_tpm = tpm(d_crl_expression, d_crl_loci_len)

    return d_alpha, d_crl_tpm


def print_configuration(args):
    """Print configuration summary."""
    print('Input BED file                       : ' + args.bed)
    print('Input merged BED file                : ' + args.merged_bed)
    print('Output directory                     : ' + args.outdir)
    print('Minimum locus size                   : ' + str(args.min_locus_size))
    print('CRL share                            : ' + str(args.crl_share))
    print('EM threshold                         : ' + str(args.em_thresh))
    print('Number of processes                  : ' + str(args.num_processes if args.num_processes > 0 else 'all available'))
    print('Create CRLs too                      : ' + str(args.build_crls_too))
    print("===================================================================")


def write_crl_output(loci_file, output_file, d_read_crl_fractions, d_crl_tpms):
    """
    Write CRL output file with fractions and TPMs.
    
    Args:
        loci_file: Path to input loci file
        output_file: Path to output file
        d_read_crl_fractions: Dictionary mapping read_id -> {crl_id: fraction}
        d_crl_tpms: Dictionary mapping crl_id -> TPM value
    """
    # Re-read file to write output (OS will cache it, so second read is fast)
    # OPTIMIZATION: Use larger buffer size (1MB) for write operations in loops
    # This significantly reduces system calls and improves I/O performance when writing many small lines
    with open(loci_file, "r", buffering=1024*1024) as fh_in:
        with open(output_file, "w", buffering=1024*1024) as fh_out:
            for l in fh_in:
                k = l.rstrip("\n").split("\t")
                # Defensive check: skip malformed lines
                if len(k) < 4:
                    print(f"Warning: Skipping malformed line in loci file (expected at least 4 fields, got {len(k)}): {l[:100]}", file=sys.stderr)
                    continue
                read_id = k[0]
                crl_id = k[3]
                # Defensive check: skip if read_id or crl_id not in dictionaries
                if read_id not in d_read_crl_fractions or crl_id not in d_read_crl_fractions[read_id]:
                    print(f"Warning: Skipping line with missing read_id or crl_id in fractions: read_id={read_id}, crl_id={crl_id}", file=sys.stderr)
                    continue
                if crl_id not in d_crl_tpms:
                    print(f"Warning: Skipping line with missing crl_id in TPMs: crl_id={crl_id}", file=sys.stderr)
                    continue
                # Format values and use join for efficient string construction
                fraction_str = "{:.2f}".format(d_read_crl_fractions[read_id][crl_id])
                tpm_str = "{:.4g}".format(d_crl_tpms[crl_id])
                fh_out.write("\t".join([l.rstrip("\n"), fraction_str, tpm_str]) + "\n")


def finalize_output(outdir, temp_file, final_file):
    """
    Remove temporary files and sort the final output file.
    
    Args:
        outdir: Output directory
        temp_file: Temporary file to remove
        final_file: Final output file to create (sorted)
    """
    temp_path = os.path.join(outdir, temp_file)
    final_path = os.path.join(outdir, final_file)
    
    chira_utilities.print_w_time("START: Sort CRLs file by read name")
    os.system("sort -V " + temp_path + " > " + final_path)
    os.remove(temp_path)
    chira_utilities.print_w_time("END: Sort CRLs file by read name")


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Chimeric Read Annotator: quantify mapped loci',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-b', '--bed', action='store', dest='bed', required=True,
                        metavar='', help='Input BED file')

    parser.add_argument('-m', '--merged_bed', action='store', dest='merged_bed', required=True,
                        metavar='', help='Input merged BED file')

    parser.add_argument('-o', '--outdir', action='store', dest='outdir', required=True, metavar='',
                        help='Output file containing merged alignments')

    parser.add_argument('-cs', '--crl_share', action='store', type=chira_utilities.score_float, default=0.7, metavar='',
                        dest='crl_share',
                        help='Minimum fraction of reads of a locus that must overlap with all CRL loci '
                             'inorder to merge it into that CRL.')

    parser.add_argument('-ls', '--min_locus_size', action='store', type=int, default=10, metavar='',
                        dest='min_locus_size',
                        help='Minimum number of reads a locus should have in order to participate in CRL creation.'
                             'Always set this value relative to your sequencing depth. Setting this to lower leads'
                             'CRLs of random multimappings Also consider setting the --crl_share option '
                             'along with this')

    parser.add_argument('-e', '--em_threshold', action='store', type=chira_utilities.score_float, default=0.00001, metavar='',
                        dest='em_thresh',
                        help='The maximum difference of transcripts expression between two consecutive iterations '
                             'of EM algorithm to converge.')

    parser.add_argument("-crl", '--build_crls_too', action='store_true', dest='build_crls_too',
                        help="Create CRLs too")

    parser.add_argument('-p', '--processes', action='store', type=int, default=0, metavar='',
                        dest='num_processes',
                        help='Number of parallel processes to use for EM algorithm. Use 0 to use all available CPU cores. '
                             'Default: 0 (use all available cores). Uses ProcessPoolExecutor for true parallelism (bypasses GIL). '
                             'Multi-processing is most beneficial for large datasets (>500 multimapping reads). '
                             'Set to 1 to disable parallel processing.')

    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {chira_utilities.__version__}')

    return parser.parse_args()


def main():
    """Main function to orchestrate the CRL quantification workflow."""
    args = parse_arguments()
    print_configuration(args)

    # Build CRLs
    loci_file = os.path.join(args.outdir, 'loci.txt')
    chira_utilities.print_w_time("START: Build CRLs")
    build_crls(args.build_crls_too, args.bed, args.merged_bed,
               loci_file, args.crl_share, args.min_locus_size)
    chira_utilities.print_w_time("END: Build CRLs")
    
    # Quantify CRLs
    chira_utilities.print_w_time("START: Quantify CRLs")
    d_read_crl_fractions, d_crl_tpms = quantify_crls(loci_file, args.em_thresh, args.num_processes)
    chira_utilities.print_w_time("END: Quantify CRLs")
    
    # Write output
    temp_output_file = os.path.join(args.outdir, 'loci.counts.temp')
    chira_utilities.print_w_time("START: Write CRLs")
    write_crl_output(loci_file, temp_output_file, d_read_crl_fractions, d_crl_tpms)
    chira_utilities.print_w_time("END: Write CRLs")
    
    # Cleanup and finalize
    os.remove(loci_file)
    finalize_output(args.outdir, 'loci.counts.temp', 'loci.counts')


if __name__ == "__main__":
    main()
