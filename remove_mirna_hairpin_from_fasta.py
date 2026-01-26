#!/usr/bin/env python
"""
Remove microRNA FASTA records from a transcriptome FASTA file.

This script extracts transcript IDs for microRNA entries from an Ensembl GTF file
and removes corresponding sequences from a transcriptome FASTA file.
"""

import argparse
import sys
import re
from Bio import SeqIO


def extract_mirna_transcript_ids(gtf_file, mirna_pattern=None):
    """
    Extract transcript IDs for microRNA entries from a GTF file.
    
    Args:
        gtf_file: Path to input GTF file
        mirna_pattern: Optional compiled regex pattern for matching miRNA in attributes.
                      If None, uses default checks (feature type and biotype only).
    
    Returns:
        Set of transcript IDs that are microRNAs
    """
    mirna_transcript_ids = set()
    
    try:
        with open(gtf_file, 'r') as fh:
            for line in fh:
                # Skip comment lines
                if line.startswith('#'):
                    continue
                
                line = line.strip()
                if not line:
                    continue
                
                fields = line.split('\t')
                if len(fields) < 9:
                    continue
                
                feature_type = fields[2]
                attributes = fields[8]
                
                # Check if this line is miRNA-related
                is_mirna = False
                
                # Check feature type
                if feature_type in ['miRNA', 'miRNA_primary_transcript']:
                    is_mirna = True
                
                # Check for gene_biotype "miRNA"
                elif re.search(r'gene_biotype\s+"miRNA"', attributes):
                    is_mirna = True
                
                # Check for transcript_biotype "miRNA"
                elif re.search(r'transcript_biotype\s+"miRNA"', attributes):
                    is_mirna = True
                
                # Check using user-provided pattern if available
                elif mirna_pattern:
                    if mirna_pattern.search(attributes):
                        is_mirna = True
                
                # If miRNA-related, extract transcript_id
                if is_mirna:
                    # Extract transcript_id from attributes
                    transcript_id_match = re.search(r'transcript_id\s+"([^"]+)"', attributes)
                    if transcript_id_match:
                        transcript_id = transcript_id_match.group(1)
                        mirna_transcript_ids.add(transcript_id)
        
        return mirna_transcript_ids
        
    except FileNotFoundError:
        print(f"Error: GTF file '{gtf_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading GTF file: {e}", file=sys.stderr)
        sys.exit(1)


def remove_mirna_from_fasta(fasta_file, output_file, mirna_transcript_ids, keep_unmatched=True):
    """
    Remove FASTA records whose transcript IDs match miRNA transcript IDs.
    
    Args:
        fasta_file: Path to input FASTA file
        output_file: Path to output FASTA file
        mirna_transcript_ids: Set of transcript IDs to remove
        keep_unmatched: If True, keep sequences that don't match any transcript ID pattern (default: True)
    """
    removed_count = 0
    kept_count = 0
    
    try:
        with open(fasta_file, 'r') as fh_in, open(output_file, 'w') as fh_out:
            for record in SeqIO.parse(fh_in, 'fasta'):
                # Extract transcript ID from FASTA header
                # FASTA headers can have various formats:
                # - ENST00000123456 (just transcript ID)
                # - ENST00000123456|other_info (transcript ID with pipe separator)
                # - >transcript_id description (transcript ID as first word)
                # - >gene_id|transcript_id|other (pipe-separated)
                
                header = record.id
                transcript_id = None
                
                # Try to extract transcript ID from header
                # Method 1: Check if header itself is in the set
                if header in mirna_transcript_ids:
                    transcript_id = header
                
                # Method 2: Check if header starts with a transcript ID (space or pipe separator)
                if not transcript_id:
                    # Try pipe-separated format: transcript_id|other
                    if '|' in header:
                        parts = header.split('|')
                        for part in parts:
                            if part in mirna_transcript_ids:
                                transcript_id = part
                                break
                    else:
                        # Try space-separated format: transcript_id description
                        first_part = header.split()[0] if ' ' in header else header
                        if first_part in mirna_transcript_ids:
                            transcript_id = first_part
                
                # Check if this record should be removed
                if transcript_id and transcript_id in mirna_transcript_ids:
                    removed_count += 1
                else:
                    # Keep the record
                    SeqIO.write(record, fh_out, 'fasta')
                    kept_count += 1
        
        print(f"Extracted {len(mirna_transcript_ids)} miRNA transcript IDs from GTF file")
        print(f"Removed {removed_count} miRNA sequences from FASTA file")
        print(f"Kept {kept_count} non-miRNA sequences")
        print(f"Output written to: {output_file}")
        
    except FileNotFoundError:
        print(f"Error: FASTA file '{fasta_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing FASTA file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Remove microRNA FASTA records from a transcriptome FASTA file using transcript IDs from a GTF file',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-g', '--gtf', action='store', dest='gtf_file', required=True,
                        metavar='', help='Input GTF file (Ensembl format)')
    parser.add_argument('-f', '--fasta', action='store', dest='fasta_file', required=True,
                        metavar='', help='Input transcriptome FASTA file')
    parser.add_argument('-o', '--output', action='store', dest='output_file', required=True,
                        metavar='', help='Output FASTA file (without microRNA sequences)')
    parser.add_argument('-p', '--pattern', action='store', dest='mirna_pattern', default=None,
                        metavar='', help='Regular expression pattern for matching miRNA in GTF attributes. '
                        'If not provided, only feature type and biotype fields are checked. '
                        'Example: \'gene_name\\s+"[^"]*(?:[Mm][Ii][Rr][_-]|[^"]*[-_][Mm][Ii][Rr][_-])\'')
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
    
    # Extract miRNA transcript IDs from GTF file
    print(f"Reading GTF file: {args.gtf_file}")
    mirna_transcript_ids = extract_mirna_transcript_ids(args.gtf_file, mirna_pattern)
    
    if not mirna_transcript_ids:
        print("Warning: No miRNA transcript IDs found in GTF file. Output will be identical to input.", file=sys.stderr)
    
    # Remove miRNA sequences from FASTA file
    print(f"Processing FASTA file: {args.fasta_file}")
    remove_mirna_from_fasta(args.fasta_file, args.output_file, mirna_transcript_ids)

