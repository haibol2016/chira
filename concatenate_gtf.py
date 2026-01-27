#!/usr/bin/env python
"""
Concatenate mature miRNA GTF file with target transcriptome GTF file.

This script combines:
1. Mature miRNA GTF file (output from gff3_to_gtf.py) - with comment lines removed
2. Target transcriptome GTF file (output from remove_mirna_hairpin_from_gtf.py) - with miRNA hairpins removed

The output is a combined GTF file suitable for use as ref2.fasta annotation in split-reference analysis.
"""

import argparse
import sys
import re


def concatenate_gtf_files(mirna_gtf, target_gtf, output_gtf, keep_target_comments=True):
    """
    Concatenate mature miRNA GTF with target transcriptome GTF.
    
    Args:
        mirna_gtf: Path to mature miRNA GTF file (from gff3_to_gtf.py)
        target_gtf: Path to target transcriptome GTF file (with miRNA hairpins removed)
        output_gtf: Path to output combined GTF file
        keep_target_comments: If True, keep comment lines from target GTF (default: True)
    """
    mirna_count = 0
    target_count = 0
    
    version_pattern = re.compile(r'transcript_id\s+"([^".]+)\.\d+"')
    
    try:
        with open(output_gtf, 'w') as fh_out:
            # First, write target transcriptome GTF (with optional comments)
            print(f"Reading target transcriptome GTF: {target_gtf}")
            with open(target_gtf, 'r') as fh_target:
                for line in fh_target:
                    if line.startswith('#'):
                        if keep_target_comments:
                            fh_out.write(line)
                    else:
                        # remove version number from transcript_id
                        line = version_pattern.sub(r'transcript_id "\1"', line)
                        fh_out.write(line)
                        target_count += 1
            
            # Then, append mature miRNA GTF (without comment lines)
            print(f"Reading mature miRNA GTF: {mirna_gtf}")
            with open(mirna_gtf, 'r') as fh_mirna:
                for line in fh_mirna:
                    # Skip comment lines from miRNA GTF
                    if line.startswith('#'):
                        continue
                    # remove version number from transcript_id
                    line = version_pattern.sub(r'transcript_id "\1"', line)
                    fh_out.write(line)
                    mirna_count += 1
        
        print(f"\nCombined GTF file created:")
        print(f"  - {target_count} lines from target transcriptome")
        print(f"  - {mirna_count} lines from mature miRNAs")
        print(f"  - Total: {target_count + mirna_count} lines")
        print(f"Output written to: {output_gtf}")
        
    except FileNotFoundError as e:
        print(f"Error: File not found: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing GTF files: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Concatenate mature miRNA GTF file with target transcriptome GTF file',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-m', '--mirna-gtf', action='store', dest='mirna_gtf', required=True,
                        metavar='', help='Mature miRNA GTF file (output from gff3_to_gtf.py)')
    parser.add_argument('-t', '--target-gtf', action='store', dest='target_gtf', required=True,
                        metavar='', help='Target transcriptome GTF file (output from remove_mirna_hairpin_from_gtf.py)')
    parser.add_argument('-o', '--output', action='store', dest='output_gtf', required=True,
                        metavar='', help='Output combined GTF file')
    parser.add_argument('--remove-target-comments', action='store_true', dest='remove_target_comments',
                        help='Remove comment lines from target GTF as well (default: keep target comments)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    concatenate_gtf_files(
        args.mirna_gtf,
        args.target_gtf,
        args.output_gtf,
        keep_target_comments=not args.remove_target_comments
    )

