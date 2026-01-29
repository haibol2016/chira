#!/usr/bin/env python
"""
Download species-specific GFF3 files from miRBase.

This script downloads the GFF3 file from miRBase for a specific species
based on the species code (e.g., hsa for human, mmu for mouse, bta for 
bovine, rno for rat, etc.).

The GFF3 file contains chromosomal coordinates of microRNAs including:
- miRNA_primary_transcript (hairpin precursor sequences)
- miRNA (mature sequences): ChiRA can make use the mature sequences for the analysis
- You don't need to convert the GFF3 file to GTF file for the analysis.
"""

import argparse
import sys
import os
import requests


def read_chromosome_mapping(mapping_file):
    """
    Read chromosome mapping from a tab-separated file.
    
    Args:
        mapping_file: Path to tab-separated file with two columns:
                      gff3_chromosome_name<tab>target_chromosome_name
    
    Returns:
        Dictionary mapping gff3 chromosome names to target chromosome names
    """
    chr_mapping = {}
    try:
        with open(mapping_file, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 2:
                    print(f"Warning: Skipping line {line_num} in mapping file (expected 2 columns, got {len(fields)})", 
                          file=sys.stderr)
                    continue
                
                gff3_chr = fields[0].strip()
                target_chr = fields[1].strip()
                
                if gff3_chr and target_chr:
                    chr_mapping[gff3_chr] = target_chr
                else:
                    print(f"Warning: Skipping line {line_num} in mapping file (empty chromosome name)", 
                          file=sys.stderr)
        
        print(f"Loaded {len(chr_mapping)} chromosome mappings from {mapping_file}")
        return chr_mapping
        
    except FileNotFoundError:
        print(f"Error: Mapping file '{mapping_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading mapping file: {e}", file=sys.stderr)
        sys.exit(1)


def convert_chromosome_names(input_file, output_file, chr_mapping):
    """
    Convert chromosome names in a GFF3 file using a mapping dictionary.
    
    Args:
        input_file: Path to input GFF3 file
        output_file: Path to output GFF3 file with converted chromosome names
        chr_mapping: Dictionary mapping source chromosome names to target names
    
    Returns:
        Number of lines converted
    """
    converted_count = 0
    total_lines = 0
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f_in, \
             open(output_file, 'w', encoding='utf-8') as f_out:
            
            for line in f_in:
                total_lines += 1
                
                # Skip comment lines and empty lines
                if line.startswith('#') or not line.strip():
                    f_out.write(line)
                    continue
                
                # GFF3 format: tab-separated, first column is chromosome/sequence name
                fields = line.rstrip('\n\r').split('\t')
                
                if len(fields) >= 1:
                    original_chr = fields[0]
                    
                    # Apply mapping if available
                    if original_chr in chr_mapping:
                        fields[0] = chr_mapping[original_chr]
                        converted_count += 1
                    
                    # Write the line (with or without conversion)
                    f_out.write('\t'.join(fields) + '\n')
                else:
                    # Write line as-is if it doesn't have expected format
                    f_out.write(line)
        
        print(f"Converted {converted_count} chromosome names out of {total_lines} total lines")
        return converted_count
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error converting chromosome names: {e}", file=sys.stderr)
        sys.exit(1)


def download_mirbase_gff3(output_file, species_code, version=None, timeout=60, chr_mapping=None):
    """
    Download species-specific GFF3 file from miRBase.
    
    Args:
        output_file: Path to save the downloaded file
        species_code: Species code (e.g., 'hsa' for human, 'mmu' for mouse, 'bta' for bovine)
        version: miRBase version number (e.g., "21"). If None, downloads CURRENT version
        timeout: Timeout in seconds for download (default: 60)
        chr_mapping: Optional dictionary mapping source chromosome names to target names
    
    Returns:
        True if successful, False otherwise
    """
    if version:
        # Use specific version URL
        # Pattern: https://www.mirbase.org/download_version_genome_files/{version}/{species_code}.gff3
        url = f"https://www.mirbase.org/download_version_genome_files/{version}/{species_code}.gff3"
    else:
        # Use CURRENT version
        # Pattern: https://www.mirbase.org/download/{species_code}.gff3
        url = f"https://www.mirbase.org/download/{species_code}.gff3"

    try:
        print(f"Downloading from: {url}")
        with requests.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()
            
            # Check if we got a valid response (not 404, etc.)
            if r.status_code != 200:
                print(f"Error: HTTP {r.status_code} - {r.reason}", file=sys.stderr)
                return False
            
            total_size = int(r.headers.get('content-length', 0))
            downloaded = 0
            first_chunk = True
            
            with open(output_file, 'w', encoding='utf-8') as f:
                for chunk in r.iter_content(chunk_size=8192, decode_unicode=True):
                    if chunk:
                        if first_chunk:
                            # Validate it's GFF3 (check for GFF version header)
                            chunk_str = chunk if isinstance(chunk, str) else chunk.decode('utf-8')
                            if not (chunk_str.startswith('##gff-version') or 
                                    chunk_str.startswith('#') or
                                    'gff' in chunk_str.lower()[:100]):
                                print(f"Warning: File may not be a valid GFF3 file (starts with: {chunk_str[:100]})", 
                                      file=sys.stderr)
                            first_chunk = False
                            chunk = chunk_str
                        
                        # Ensure chunk is string for writing
                        if isinstance(chunk, bytes):
                            chunk = chunk.decode('utf-8')
                        
                        f.write(chunk)
                        downloaded += len(chunk.encode('utf-8'))
                        
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            print(f"\rProgress: {percent:.1f}% ({downloaded}/{total_size} bytes)", end='', flush=True)
            
            if total_size > 0:
                print()
            print(f"Successfully downloaded to: {output_file}")
            
            # Apply chromosome name mapping if provided
            if chr_mapping:
                temp_file = output_file + '.tmp'
                os.rename(output_file, temp_file)
                print(f"\nApplying chromosome name mapping...")
                convert_chromosome_names(temp_file, output_file, chr_mapping)
                os.remove(temp_file)
                print(f"Chromosome names converted in output file")
            
            return True
                    
    except requests.exceptions.RequestException as e:
        print(f"Failed to download: {e}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Download species-specific GFF3 file from miRBase',
        usage='%(prog)s [-h] [-v,--version]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-s', '--species', action='store', dest='species_code', required=True,
                        metavar='', help='Species code (e.g., hsa=human, mmu=mouse, bta=bovine, rno=rat)')
    parser.add_argument('-o', '--output', action='store', dest='output_file', required=True,
                        metavar='', help='Output GFF3 file for species-specific microRNA annotations')
    parser.add_argument('--mirbase-version', action='store', dest='mirbase_version', default=None,
                        metavar='', help='miRBase version number (e.g., "21"). If not specified, downloads CURRENT version')
    parser.add_argument('--timeout', action='store', type=int, dest='timeout', default=60,
                        metavar='', help='Download timeout in seconds (default: 60)')
    parser.add_argument('-m', '--chromosome-mapping', action='store', dest='chr_mapping_file', default=None,
                        metavar='', help='Tab-separated file with two columns: gff3_chromosome_name<tab>target_chromosome_name. '
                        'Chromosome names in the GFF3 file will be converted to target names in the output.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Read chromosome mapping if provided
    chr_mapping = None
    if args.chr_mapping_file:
        chr_mapping = read_chromosome_mapping(args.chr_mapping_file)
    
    print(f"Downloading GFF3 file from miRBase...")
    print(f"Species: {args.species_code}")
    print(f"Version: {args.mirbase_version or 'CURRENT'}")
    if chr_mapping:
        print(f"Chromosome mapping: {args.chr_mapping_file}")
    
    if not download_mirbase_gff3(args.output_file, args.species_code, args.mirbase_version, args.timeout, chr_mapping):
        print("Error: Failed to download GFF3 file from miRBase", file=sys.stderr)
        sys.exit(1)
    
    print(f"Output written to: {args.output_file}")


if __name__ == "__main__":
    main()

