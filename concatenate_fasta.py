#!/usr/bin/env python
"""
Concatenate cDNA and miRNA-free ncRNA FASTA files and remove transcript version numbers.

This script combines cDNA and miRNA-free ncRNA FASTA files into a single output file,
removing version numbers from transcript IDs in the headers.
Supports both Ensembl style (space-separated) and GENCODE style (pipe-separated) headers.
"""

import argparse
import sys
import re
from Bio import SeqIO


def remove_version_from_header(header, split_pattern = None, version_pattern = None):
    """
    Remove version number from FASTA header transcript ID.
    
    Args:
        header: FASTA header string (without '>')
    
    Returns:
        Header with version number removed
    """
    header = split_pattern.split(header.strip())[0] if split_pattern else header.strip()
    return version_pattern.sub('', header) if version_pattern else header


def concatenate_fasta_files(cdna_file, ncrna_file, output_file):
    """
    Concatenate cDNA and ncRNA FASTA files, removing version numbers from headers.
    
    Args:
        cdna_file: Path to input cDNA FASTA file
        ncrna_file: Path to input ncRNA FASTA file
        output_file: Path to output concatenated FASTA file
    """
    total_sequences = 0
    cdna_count = 0
    ncrna_count = 0
    
    split_pattern = re.compile(r'[\s+|]')
    version_pattern = re.compile(r'\.\d+$')

    try:
        with open(output_file, 'w') as fh_out:
            # Process cDNA file
            if cdna_file:
                print(f"Processing cDNA file: {cdna_file}")
                try:
                    with open(cdna_file, 'r') as fh_in:
                        for record in SeqIO.parse(fh_in, 'fasta'):
                            # Remove version number from header
                            record.id = remove_version_from_header(record.id, split_pattern, version_pattern)
                            record.description = record.id  # Update description to match
                            
                            SeqIO.write(record, fh_out, 'fasta')
                            cdna_count += 1
                            total_sequences += 1
                except FileNotFoundError:
                    print(f"Warning: cDNA file '{cdna_file}' not found, skipping...", file=sys.stderr)
                except Exception as e:
                    print(f"Error processing cDNA file: {e}", file=sys.stderr)
                    sys.exit(1)
            
            # Process ncRNA file
            if ncrna_file:
                print(f"Processing ncRNA file: {ncrna_file}")
                try:
                    with open(ncrna_file, 'r') as fh_in:
                        for record in SeqIO.parse(fh_in, 'fasta'):
                            # Remove version number from header
                            record.id = remove_version_from_header(record.id, split_pattern, version_pattern)
                            record.description = record.id  # Update description to match
                            
                            SeqIO.write(record, fh_out, 'fasta')
                            ncrna_count += 1
                            total_sequences += 1
                except FileNotFoundError:
                    print(f"Warning: ncRNA file '{ncrna_file}' not found, skipping...", file=sys.stderr)
                except Exception as e:
                    print(f"Error processing ncRNA file: {e}", file=sys.stderr)
                    sys.exit(1)
        
        print(f"\nConcatenation complete:")
        if cdna_file:
            print(f"  cDNA sequences: {cdna_count}")
        if ncrna_file:
            print(f"  ncRNA sequences: {ncrna_count}")
        print(f"  Total sequences: {total_sequences}")
        print(f"  Output written to: {output_file}")
        
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Concatenate cDNA and ncRNA FASTA files and remove transcript version numbers',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-c', '--cdna', action='store', dest='cdna_file', default=None,
                        metavar='', help='Input cDNA FASTA file (optional)')
    parser.add_argument('-n', '--ncrna', action='store', dest='ncrna_file', default=None,
                        metavar='', help='Input ncRNA FASTA file (optional)')
    parser.add_argument('-o', '--output', action='store', dest='output_file', required=True,
                        metavar='', help='Output concatenated FASTA file')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Check that at least one input file is provided
    if not args.cdna_file and not args.ncrna_file:
        print("Error: At least one input file (--cdna or --ncrna) must be provided", file=sys.stderr)
        sys.exit(1)
    
    concatenate_fasta_files(args.cdna_file, args.ncrna_file, args.output_file)

