#!/usr/bin/env python
"""
Remove all included microRNA entries from an Ensembl GTF file.

This script filters out all lines in a GTF file that are related to microRNAs,
including miRNA genes, miRNA transcripts, and miRNA features.
"""

import argparse
import sys
import re


def is_mirna_line(line, mirna_pattern=None):
    """
    Check if a GTF line is related to microRNA.
    
    Args:
        line: A line from a GTF file
        mirna_pattern: Optional compiled regex pattern for matching miRNA in attributes.
                      If None, uses default checks (feature type and biotype only).
    
    Returns:
        True if the line is related to microRNA, False otherwise
    """
    # Skip comment lines
    if line.startswith('#'):
        return False
    
    # Split the line into fields
    fields = line.strip().split('\t')
    if len(fields) < 9:
        return False
    
    # Check feature type
    feature_type = fields[2]
    if feature_type in ['miRNA', 'miRNA_primary_transcript']:
        return True
    
    # Check attributes for miRNA biotype
    attributes = fields[8]
    
    # Check for gene_biotype "miRNA"
    if re.search(r'gene_biotype\s+"miRNA"', attributes):
        return True
    
    # Check for transcript_biotype "miRNA" (if present)
    if re.search(r'transcript_biotype\s+"miRNA"', attributes):
        return True
    
    # Check for miRNA using user-provided pattern if available
    if mirna_pattern:
        if mirna_pattern.search(attributes):
            return True
    
    return False


def remove_mirna_from_gtf(input_gtf, output_gtf, keep_comments=True, mirna_pattern=None):
    """
    Remove all microRNA entries from a GTF file.
    
    Args:
        input_gtf: Path to input GTF file
        output_gtf: Path to output GTF file
        keep_comments: If True, keep comment lines in output (default: True)
        mirna_pattern: Optional compiled regex pattern for matching miRNA in attributes.
                      If None, only checks feature type and biotype fields.
    """
    removed_count = 0
    kept_count = 0
    
    try:
        with open(input_gtf, 'r') as fh_in, open(output_gtf, 'w') as fh_out:
            for line in fh_in:
                # Keep comment lines if requested
                if line.startswith('#'):
                    if keep_comments:
                        fh_out.write(line)
                    continue
                
                # Check if line is miRNA-related
                if is_mirna_line(line, mirna_pattern):
                    removed_count += 1
                else:
                    fh_out.write(line)
                    kept_count += 1
        
        print(f"Removed {removed_count} microRNA-related lines")
        print(f"Kept {kept_count} non-microRNA lines")
        print(f"Output written to: {output_gtf}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_gtf}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing GTF file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Remove all microRNA entries from an Ensembl GTF file',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-i', '--input', action='store', dest='input_gtf', required=True,
                        metavar='', help='Input GTF file')
    parser.add_argument('-o', '--output', action='store', dest='output_gtf', required=True,
                        metavar='', help='Output GTF file (without microRNA entries)')
    parser.add_argument('-p', '--pattern', action='store', dest='mirna_pattern', default=None,
                        metavar='', help='Regular expression pattern for matching miRNA in GTF attributes. '
                        'If not provided, only feature type and biotype fields are checked. '
                        'Example: \'gene_name\\s+"[^"]*(?:[Mm][Ii][Rr][_-]|[^"]*[-_][Mm][Ii][Rr][_-])\'')
    parser.add_argument('--remove-comments', action='store_true', dest='remove_comments',
                        help='Remove comment lines from output (default: keep comments)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Compile regex pattern if provided
    mirna_pattern = None
    if args.mirna_pattern:
        try:
            mirna_pattern = re.compile(args.mirna_pattern)
        except re.error as e:
            print(f"Error: Invalid regular expression pattern: {e}", file=sys.stderr)
            sys.exit(1)
    
    remove_mirna_from_gtf(
        args.input_gtf,
        args.output_gtf,
        keep_comments=not args.remove_comments,
        mirna_pattern=mirna_pattern
    )

