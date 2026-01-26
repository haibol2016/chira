#!/usr/bin/env python
"""
Download species-specific mature miRNA sequences from miRBase.

This script downloads the mature miRNA FASTA file from miRBase and extracts
sequences for a specific species based on the species code (e.g., hsa for 
human, mmu for mouse, bta for bovine, rno for rat, etc.).
"""

import argparse
import sys
import os
import gzip
import urllib.request
import urllib.error


def download_mirbase_mature(output_file, version=None, timeout=30):
    """
    Download mature miRNA FASTA file from miRBase.
    
    Args:
        output_file: Path to save the downloaded file (will be gzipped)
        version: miRBase version (e.g., "22.1"). If None, downloads CURRENT version
        timeout: Timeout in seconds for download (default: 30)
    
    Returns:
        True if successful, False otherwise
    """
    if version:
        # Use specific version URL
        url1 = f"https://www.mirbase.org/ftp/{version}/mature.fa.gz"
        url2 = f"ftp://mirbase.org/pub/mirbase/{version}/mature.fa.gz"
    else:
        # Use CURRENT version
        url1 = "https://www.mirbase.org/ftp/CURRENT/mature.fa.gz"
        url2 = "ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz"
    
    # Try primary URL first
    for url in [url1, url2]:
        try:
            print(f"Downloading from: {url}")
            urllib.request.urlretrieve(url, output_file, timeout=timeout)
            print(f"Successfully downloaded to: {output_file}")
            return True
        except urllib.error.URLError as e:
            print(f"Failed to download from {url}: {e}", file=sys.stderr)
            continue
        except Exception as e:
            print(f"Error downloading: {e}", file=sys.stderr)
            continue
    
    return False


def extract_species_mirnas(input_file, output_file, species_code):
    """
    Extract species-specific mature miRNA sequences from miRBase FASTA file.
    
    Args:
        input_file: Path to input FASTA file (can be gzipped or uncompressed)
        output_file: Path to output FASTA file (species-specific sequences)
        species_code: Species code (e.g., 'hsa' for human, 'mmu' for mouse, 'bta' for bovine)
    
    Returns:
        Number of sequences extracted
    """
    sequence_count = 0
    current_header = None
    current_sequence = []
    in_species_sequence = False
    
    # Determine if input is gzipped
    open_func = gzip.open if input_file.endswith('.gz') else open
    mode = 'rt' if input_file.endswith('.gz') else 'r'
    
    try:
        with open_func(input_file, mode) as fh_in, open(output_file, 'w') as fh_out:
            for line in fh_in:
                line = line.rstrip('\n\r')
                
                if line.startswith('>'):
                    # Write previous sequence if it was species-specific
                    if in_species_sequence and current_header and current_sequence:
                        fh_out.write(current_header + '\n')
                        fh_out.write(''.join(current_sequence) + '\n')
                        sequence_count += 1
                    
                    # Check if this header is for the target species
                    # miRBase format: >MIMAT0000061 hsa-let-7a-5p Homo sapiens let-7a-5p
                    # Pattern: -species_code- (e.g., -hsa-, -mmu-, -bta-)
                    current_header = line
                    if f'-{species_code}-' in line:
                        in_species_sequence = True
                        current_sequence = []
                    else:
                        in_species_sequence = False
                        current_sequence = []
                else:
                    # Sequence line
                    if in_species_sequence:
                        current_sequence.append(line)
            
            # Write last sequence if it was species-specific
            if in_species_sequence and current_header and current_sequence:
                fh_out.write(current_header + '\n')
                fh_out.write(''.join(current_sequence) + '\n')
                sequence_count += 1
        
        return sequence_count
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing FASTA file: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Download species-specific mature miRNA sequences from miRBase',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-s', '--species', action='store', dest='species_code', required=True,
                        metavar='', help='Species code (e.g., hsa=human, mmu=mouse, bta=bovine, rno=rat)')
    parser.add_argument('-o', '--output', action='store', dest='output_file', required=True,
                        metavar='', help='Output FASTA file for species-specific mature miRNAs')
    parser.add_argument('--mirbase-version', action='store', dest='mirbase_version', default=None,
                        metavar='', help='miRBase version (e.g., "22.1"). If not specified, downloads CURRENT version')
    parser.add_argument('--keep-full', action='store_true', dest='keep_full',
                        help='Keep the full mature.fa.gz file after extraction (default: remove)')
    parser.add_argument('--timeout', action='store', type=int, dest='timeout', default=30,
                        metavar='', help='Download timeout in seconds (default: 30)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Temporary file for downloaded mature.fa.gz
    temp_file = 'mature.fa.gz'
    
    # Download mature miRNA FASTA from miRBase
    print(f"Downloading mature miRNA sequences from miRBase...")
    if args.mirbase_version:
        print(f"Version: {args.mirbase_version}")
    else:
        print("Version: CURRENT")
    
    if not download_mirbase_mature(temp_file, args.mirbase_version, args.timeout):
        print("Error: Failed to download mature miRNA sequences from miRBase", file=sys.stderr)
        sys.exit(1)
    
    # Extract species-specific sequences
    print(f"\nExtracting {args.species_code} (species code) mature miRNAs...")
    sequence_count = extract_species_mirnas(temp_file, args.output_file, args.species_code)
    
    print(f"Extracted {sequence_count} mature miRNA sequences for {args.species_code}")
    print(f"Output written to: {args.output_file}")
    
    # Clean up temporary file unless --keep-full is specified
    if not args.keep_full:
        try:
            os.remove(temp_file)
            print(f"Removed temporary file: {temp_file}")
        except Exception as e:
            print(f"Warning: Could not remove temporary file {temp_file}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()

