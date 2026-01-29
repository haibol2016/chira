#!/usr/bin/env python
import os
import sys
from collections import defaultdict
import argparse
import chira_utilities
import copy


def build_crls(build_crls_too, bed, merged_bed, crl_file, crl_share_cutoff, min_locus_size):
    l_locipos = []
    l_locireads = []
    d_readlocus_transcripts = defaultdict(list)
    with open(merged_bed) as fh_merged_bed:
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
            # lower and uppder bounds for filtering crls based on their size
            lower_bound = len(l_locusreads) * (1 - crl_share_cutoff)
            upper_bound = len(l_locusreads) / (1 - crl_share_cutoff)
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
    with open(bed) as fh_bed:
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

    with open(crl_file, "w") as fh_crl_file:
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
                for segmentid in sorted(l_locusreads):
                    segment_key = segmentid + str(locusid)
                    for transcriptid_pos in sorted(d_readlocus_transcripts[segment_key]):
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


def em(d_alpha, d_rho, library_size, l_multimap_readids, em_threshold):
    i = 1
    sum_of_rho_diff = float('inf')
    while sum_of_rho_diff > em_threshold and i <= 1000:
        if i % 10 == 0:  # Print every 10 iterations instead of every iteration
            print("iteration: " + str(i))
        # E-step
        # Shallow copy is sufficient since d_rho contains only immutable float values
        # Using copy() preserves defaultdict behavior, but dict() is slightly faster
        d_rho_old = dict(d_rho)  # Shallow copy - safe because values are immutable floats
        d_rho.clear()
        # iterate through multi-mapped read segments
        # Note: l_multimap_readids only contains reads with multiple CRLs (len > 1)
        for readid in l_multimap_readids:
            read_crls = d_alpha[readid]
            # Pre-compute total abundance to avoid repeated lookups
            total_crl_rel_abundancy = sum(d_rho_old[crlid] for crlid in read_crls)
            if total_crl_rel_abundancy > 0:
                inv_total = 1.0 / total_crl_rel_abundancy
                for crlid in read_crls:
                    d_alpha[readid][crlid] = d_rho_old[crlid] * inv_total

        for readid in d_alpha.keys():
            for crlid in d_alpha[readid]:
                d_rho[crlid] += d_alpha[readid][crlid]
        # relative abundancies
        # M-step
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
    # Cache sorted keys to avoid repeated sorting (used in two loops)
    sorted_crlids = sorted(d_crl_expression.keys())
    for crlid in sorted_crlids:
        crl_expression = d_crl_expression[crlid]
        crl_len = chira_utilities.median(sorted(d_crl_loci_len[crlid].values())) / 1000.0  # length in kbs
        if crl_len > 0:  # Avoid division by zero
            rpk = crl_expression / crl_len
            d_crl_tpm[crlid] = rpk
            total_rpk += rpk
    millions_of_rpk = total_rpk / 1000000.0
    if millions_of_rpk > 0:  # Avoid division by zero
        inv_millions = 1.0 / millions_of_rpk
        for crlid in sorted_crlids:
            crl_tpm = d_crl_tpm[crlid] * inv_millions
            d_crl_tpm[crlid] = crl_tpm
    return d_crl_tpm


def quantify_crls(crl_file, em_threshold):
    d_crl_loci_len = defaultdict(lambda: defaultdict(int))
    l_multimap_readids = []
    d_alpha = defaultdict(lambda: defaultdict(float))
    d_rho = defaultdict(float)

    # Use context manager for file handling (more efficient)
    with open(crl_file, "r") as fh_crl_file:
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
    for readid in d_alpha.keys():
        for crlid in d_alpha[readid].keys():
            d_alpha[readid][crlid] = 1 / float(len(d_alpha[readid]))
            d_rho[crlid] += d_alpha[readid][crlid]
        if len(d_alpha[readid]) > 1:
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

    d_res = em(d_alpha, d_rho, library_size, l_multimap_readids, em_threshold)

    d_crl_expression = defaultdict(float)
    for readid in d_res.keys():
        for crlid in d_res[readid]:
            d_crl_expression[crlid] += d_res[readid][crlid]

    d_crl_tpm = tpm(d_crl_expression, d_crl_loci_len)

    return d_alpha, d_crl_tpm


if __name__ == "__main__":

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

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.4.4')

    args = parser.parse_args()

    print('Input BED file                       : ' + args.bed)
    print('Input merged BED file                : ' + args.merged_bed)
    print('Output directory                     : ' + args.outdir)
    print('Minimum locus size                   : ' + str(args.min_locus_size))
    print('CRL share                            : ' + str(args.crl_share))
    print('EM threshold                         : ' + str(args.em_thresh))
    print('Create CRLs too                      : ' + str(args.build_crls_too))
    print("===================================================================")

    chira_utilities.print_w_time("START: Build CRLs")
    build_crls(args.build_crls_too, args.bed, args.merged_bed,
               os.path.join(args.outdir, 'loci.txt'), args.crl_share, args.min_locus_size)
    chira_utilities.print_w_time("END: Build CRLs")
    chira_utilities.print_w_time("START: Quantify CRLs")
    d_read_crl_fractions, d_crl_tpms = quantify_crls(os.path.join(args.outdir, 'loci.txt'), args.em_thresh)
    chira_utilities.print_w_time("END: Quantify CRLs")
    chira_utilities.print_w_time("START: Write CRLs")
    # Re-read file to write output (OS will cache it, so second read is fast)
    with open(os.path.join(args.outdir, 'loci.txt'), "r") as fh_in:
        with open(os.path.join(args.outdir, 'loci.counts.temp'), "w") as fh_out:
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
                fh_out.write("\t".join([l.strip("\n"),
                                        "{:.2f}".format(d_read_crl_fractions[read_id][crl_id]),
                                        "{:.4g}".format(d_crl_tpms[crl_id])]) + "\n")
    chira_utilities.print_w_time("END: Write CRLs")
    os.remove(os.path.join(args.outdir, 'loci.txt'))
    chira_utilities.print_w_time("START: Sort CRLs file by read name")
    os.system("sort -V " + os.path.join(args.outdir, 'loci.counts.temp') + " > " + os.path.join(args.outdir, 'loci.counts'))
    os.remove(os.path.join(args.outdir, 'loci.counts.temp'))
    chira_utilities.print_w_time("END: Sort CRLs file by read name")
