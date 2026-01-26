#!/usr/bin/env python

"""
Convert miRBase GFF3 format to GTF format following ENSEMBL conventions and perform coordinate liftover if needed.

The input GFF3 file is expected to be in miRBase format.
The output GTF file will be in ENSEMBL format.
The coordinate liftover is performed using the pyliftover library.
The chain file is expected to be in the chain format.
The chain file can be downloaded from UCSC Genome Browser.
The chain file is required for coordinate liftover.
The chromosome mapping file is optional and can be used to map chromosome names to new names.

"""

import argparse
import sys

# Import pyliftover (required for liftover functionality)
try:
    from pyliftover import LiftOver
except ImportError:
    LiftOver = None


def _convert_coordinates_pyliftover(lo, chrom, start, end, strand):
    """
    Convert coordinates using pyliftover library.
    
    Args:
        lo: LiftOver object from pyliftover
        chrom: Chromosome name
        start: Start position (1-based)
        end: End position (1-based)
        strand: Strand (+ or -)
    
    Returns:
        Tuple of (new_chrom, new_start, new_end, new_strand) or None if conversion fails
    """
    try:
        # pyliftover uses 0-based coordinates
        start_0based = start - 1
        end_0based = end - 1
        
        # Convert start position
        start_converted = lo.convert_coordinate(chrom, start_0based)
        if not start_converted:
            return None
        
        # Convert end position
        end_converted = lo.convert_coordinate(chrom, end_0based)
        if not end_converted:
            return None
        
        # Use the first conversion result (most likely)
        new_chrom_start, new_start_0based, new_strand_start, score_start = start_converted[0]
        new_chrom_end, new_end_0based, new_strand_end, score_end = end_converted[0]
        
        # Ensure both positions map to the same chromosome
        if new_chrom_start != new_chrom_end:
            return None
        
        # Convert back to 1-based coordinates
        new_start = new_start_0based + 1
        new_end = new_end_0based + 1
        
        # Preserve original strand if liftover doesn't change it
        new_strand = strand if new_strand_start == '+' else ('-' if strand == '+' else '+')
        
        return (new_chrom_start, new_start, new_end, new_strand)
    except Exception:
        return None


def gff3_to_gtf(gff3_file, gtf_file, chromosome_mapping_file=None,
                source_genome_version=None, target_genome_version=None,
                chain_file=None):
    """
    Convert miRBase GFF3 format to GTF format following ENSEMBL conventions.
    
    This function converts GFF3 files (which use ID and Parent attributes)
    to GTF format (which uses gene_id and transcript_id attributes).
    The output follows ENSEMBL GTF format specifications.
    
    Args:
        gff3_file: Path to input GFF3 file
        gtf_file: Path to output GTF file
        chromosome_mapping_file: Path to chromosome mapping file (default: None)
            - Two-column TSV/CSV file with tab delimiter: original_name<tab>new_name
            - Example file content:
                chr1    1
                chrX    X
                chrM    MT
                chrMT   MT
            - If provided, chromosome names in the output will be converted according to the mapping
            - If None, chromosome names are kept unchanged
        source_genome_version: Source genome version (e.g., 'hg19', 'GRCh37')
        target_genome_version: Target genome version (e.g., 'hg38', 'GRCh38')
        chain_file: Path to chain file for coordinate liftover (required if liftover is enabled)
            - Chain files can be downloaded from UCSC Genome Browser
            - Requires pyliftover library: pip install pyliftover
    
    The conversion handles:
    - ID/Parent relationships -> gene_id/transcript_id
    - GFF3 attribute format (key=value) -> GTF format (key "value")
    - Hierarchical structure (gene -> transcript -> exon/CDS/UTR)
    - Feature type mapping (e.g., miRNA -> transcript with miRNA exon)
    - Chromosome name conversion (if chromosome_mapping_file provided)
    - Genome version liftover (if source/target versions and chain file provided)
    
    Note: This function expects standard GFF3 format with ID and Parent attributes.
    """
    # Setup liftover if needed
    lo = None
    use_liftover = False
    
    if source_genome_version and target_genome_version and chain_file:
        use_liftover = True
        if LiftOver is None:
            raise ValueError("pyliftover is required for coordinate liftover. Install with: pip install pyliftover")
        try:
            lo = LiftOver(chain_file)
            print(f"Using pyliftover for coordinate conversion: {source_genome_version} -> {target_genome_version}")
        except Exception as e:
            raise ValueError(f"Failed to initialize pyliftover with chain file '{chain_file}': {e}")
    
    # Read chromosome mapping file if provided
    chromosome_mapping = {}
    if chromosome_mapping_file:
        try:
            with open(chromosome_mapping_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            chromosome_mapping[parts[0].strip()] = parts[1].strip()
        except Exception as e:
            raise ValueError(f"Error reading chromosome mapping file '{chromosome_mapping_file}': {e}")
    
    # Dictionaries to track ID mappings
    id_to_gene_id = {}  # Maps any ID to its gene_id
    id_to_transcript_id = {}  # Maps transcript ID to transcript_id
    id_to_gene_name = {}  # Maps gene ID to gene name
    id_to_biotype = {}  # Maps gene ID to biotype
    id_to_alias = {}  # Maps ID to alias
    id_to_name = {}  # Maps ID to name
    
    # First pass: build ID relationships
    with open(gff3_file, 'r') as fh_in:
        for line in fh_in:
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            attributes_str = fields[8]
            
            # Parse attributes (GFF3 format: key=value;key2=value2)
            attrs = {}
            for attr in attributes_str.split(';'):
                attr = attr.strip()
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attrs[key.strip()] = value.strip()
            
            # Track gene-level information
            if feature_type == 'gene' and 'ID' in attrs:
                gene_id = attrs['ID']
                id_to_gene_id[gene_id] = gene_id
                if 'Name' in attrs:
                    id_to_gene_name[gene_id] = attrs['Name']
                if 'biotype' in attrs:
                    id_to_biotype[gene_id] = attrs['biotype']
                elif 'gene_biotype' in attrs:
                    id_to_biotype[gene_id] = attrs['gene_biotype']
            
            # Track transcript-level information
            elif feature_type in ['transcript', 'mRNA', 'miRNA', 'tRNA', 'miRNA_primary_transcript']:
                if 'ID' in attrs:
                    transcript_id = attrs['ID']
                    id_to_transcript_id[transcript_id] = transcript_id
                    # Store alias and name
                    if 'Alias' in attrs:
                        id_to_alias[transcript_id] = attrs['Alias']
                    if 'Name' in attrs:
                        id_to_name[transcript_id] = attrs['Name']
                    # Link transcript to gene
                    if 'Parent' in attrs:
                        parent_id = attrs['Parent']
                        if parent_id in id_to_gene_id:
                            id_to_gene_id[transcript_id] = id_to_gene_id[parent_id]
                        else:
                            # If parent is not a gene, try to find gene
                            id_to_gene_id[transcript_id] = parent_id
                    else:
                        # For miRNA_primary_transcript, treat as both gene and transcript
                        if feature_type == 'miRNA_primary_transcript':
                            id_to_gene_id[transcript_id] = transcript_id
                            # Get biotype from feature type
                            if 'biotype' in attrs:
                                id_to_biotype[transcript_id] = attrs['biotype']
                            else:
                                id_to_biotype[transcript_id] = 'miRNA'
            
            # Handle Derives_from relationship (miRBase format)
            if 'Derives_from' in attrs:
                parent_id = attrs['Derives_from']
                feature_id = attrs.get('ID', '')
                if feature_id and parent_id in id_to_gene_id:
                    # Link derived feature to parent
                    id_to_gene_id[feature_id] = id_to_gene_id[parent_id]
                    id_to_transcript_id[feature_id] = parent_id
                    # Store alias and name for derived features
                    if 'Alias' in attrs:
                        id_to_alias[feature_id] = attrs['Alias']
                    if 'Name' in attrs:
                        id_to_name[feature_id] = attrs['Name']
            
            # Track child features (exon, CDS, UTR) and link to parents
            elif feature_type in ['exon', 'CDS', 'UTR', 'five_prime_UTR', 'three_prime_UTR']:
                if 'Parent' in attrs:
                    parent_id = attrs['Parent']
                    # If parent is a transcript, link to its gene
                    if parent_id in id_to_transcript_id:
                        if parent_id in id_to_gene_id:
                            id_to_gene_id[attrs.get('ID', parent_id)] = id_to_gene_id[parent_id]
                    # If parent is a gene, we'll handle it in second pass
                    elif parent_id in id_to_gene_id:
                        # Create a transcript ID from the feature ID or use parent
                        feature_id = attrs.get('ID', '')
                        if feature_id:
                            id_to_transcript_id[feature_id] = feature_id
                            id_to_gene_id[feature_id] = id_to_gene_id[parent_id]
    
    # Second pass: write GTF format
    with open(gff3_file, 'r') as fh_in, open(gtf_file, 'w') as fh_out:
        exon_number = {}  # Track exon numbers per transcript
        
        for line in fh_in:
            line = line.strip()
            # Skip comments (except keep version info as comment)
            if line.startswith('##'):
                if line.startswith('##gff-version'):
                    # Convert to GTF comment
                    fh_out.write('# GTF converted from GFF3\n')
                continue
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            
            seqname = fields[0]
            source = fields[1]
            feature_type = fields[2]
            start_str = fields[3]
            end_str = fields[4]
            score = fields[5] if fields[5] != '.' else '.'
            strand = fields[6] if fields[6] != '.' else '.'
            frame = fields[7] if fields[7] != '.' else '.'
            attributes_str = fields[8]
            
            # Parse coordinates
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                # Skip lines with invalid coordinates
                continue
            
            # Perform coordinate liftover if needed
            if use_liftover:
                converted = _convert_coordinates_pyliftover(lo, seqname, start, end, strand)
                if converted:
                    seqname, start, end, strand = converted
                else:
                    # Skip features that can't be lifted over
                    print(f"Warning: Could not lift over {seqname}:{start}-{end}, skipping", file=sys.stderr)
                    continue
            
            # Convert chromosome name if mapping provided (after liftover)
            if chromosome_mapping:
                if seqname in chromosome_mapping:
                    seqname = chromosome_mapping[seqname]
                elif chromosome_mapping_file:
                    # Only raise error if chromosome_mapping_file was explicitly provided
                    raise ValueError(f"Chromosome name {seqname} not found in chromosome mapping file")
            
            # Parse attributes
            attrs = {}
            for attr in attributes_str.split(';'):
                attr = attr.strip()
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attrs[key.strip()] = value.strip()
            
            # Skip gene features (GTF doesn't have gene lines, info is in transcript/exon)
            if feature_type == 'gene':
                continue
            
            # Determine gene_id and transcript_id
            gene_id = None
            transcript_id = None
            feature_id = attrs.get('ID', '')
            
            if feature_type in ['transcript', 'mRNA', 'miRNA', 'tRNA', 'miRNA_primary_transcript']:
                transcript_id = feature_id
                # Get gene_id from parent or direct mapping
                if 'Parent' in attrs:
                    parent_id = attrs['Parent']
                    gene_id = id_to_gene_id.get(parent_id, parent_id)
                elif 'Derives_from' in attrs:
                    # Handle Derives_from (miRBase format)
                    parent_id = attrs['Derives_from']
                    gene_id = id_to_gene_id.get(parent_id, parent_id)
                    transcript_id = parent_id  # Use parent as transcript_id
                else:
                    gene_id = id_to_gene_id.get(transcript_id, transcript_id)
            elif feature_type in ['exon', 'CDS', 'UTR', 'five_prime_UTR', 'three_prime_UTR']:
                if 'Parent' in attrs:
                    parent_id = attrs['Parent']
                    # Check if parent is a transcript
                    if parent_id in id_to_transcript_id:
                        transcript_id = parent_id
                        gene_id = id_to_gene_id.get(transcript_id, parent_id)
                    # If parent is a gene, use feature ID as transcript_id (create implicit transcript)
                    elif parent_id in id_to_gene_id:
                        gene_id = parent_id
                        # Use feature ID as transcript_id, or create one from gene_id
                        feature_id = attrs.get('ID', '')
                        if feature_id and feature_id in id_to_transcript_id:
                            transcript_id = feature_id
                        else:
                            # Create a transcript ID from gene_id
                            transcript_id = gene_id + '_transcript'
                            id_to_transcript_id[transcript_id] = transcript_id
                            id_to_gene_id[transcript_id] = gene_id
                    else:
                        # Try to find gene_id through any available mapping
                        gene_id = id_to_gene_id.get(parent_id)
                        transcript_id = parent_id if gene_id else None
                elif 'Derives_from' in attrs:
                    # Handle Derives_from (miRBase format for mature miRNAs)
                    parent_id = attrs['Derives_from']
                    transcript_id = parent_id
                    gene_id = id_to_gene_id.get(parent_id, parent_id)
            
            # Skip if we can't determine gene_id or transcript_id
            if not gene_id or not transcript_id:
                continue
            
            # Get additional attributes
            gene_name = id_to_gene_name.get(gene_id, gene_id)
            # Determine biotype - check feature type first, then stored biotype
            if feature_type == 'miRNA' or feature_type == 'miRNA_primary_transcript':
                biotype = 'miRNA'
            elif feature_type == 'tRNA':
                biotype = 'tRNA'
            else:
                biotype = id_to_biotype.get(gene_id, 'protein_coding')
            
            # Get alias and name for this feature
            feature_alias = id_to_alias.get(feature_id, '')
            feature_name = id_to_name.get(feature_id, '')
            # Also check transcript-level alias/name
            if not feature_alias:
                feature_alias = id_to_alias.get(transcript_id, '')
            if not feature_name:
                feature_name = id_to_name.get(transcript_id, '')
            
            # Handle exon numbering for GTF
            if feature_type == 'exon':
                if transcript_id not in exon_number:
                    exon_number[transcript_id] = 0
                exon_number[transcript_id] += 1
                exon_num = exon_number[transcript_id]
            else:
                exon_num = None
            
            # Build GTF attributes (ENSEMBL format: key "value";)
            gtf_attrs = []
            gtf_attrs.append('gene_id "{}"'.format(gene_id))
            gtf_attrs.append('transcript_id "{}"'.format(transcript_id))
            
            # Use feature-level name if available, otherwise use gene-level name
            if feature_name:
                gtf_attrs.append('gene_name "{}"'.format(feature_name))
            elif gene_name and gene_name != gene_id:
                gtf_attrs.append('gene_name "{}"'.format(gene_name))
            
            if biotype:
                gtf_attrs.append('gene_biotype "{}"'.format(biotype))
            
            # Add alias if available
            if feature_alias:
                gtf_attrs.append('gene_alias "{}"'.format(feature_alias))
            
            if exon_num is not None:
                gtf_attrs.append('exon_number "{}"'.format(exon_num))
            
            # Add feature-specific attributes
            if feature_type == 'CDS':
                gtf_attrs.append('feature_type "CDS"')
            elif feature_type in ['UTR', 'five_prime_UTR', 'three_prime_UTR']:
                if feature_type == 'five_prime_UTR':
                    gtf_attrs.append('feature_type "5UTR"')
                elif feature_type == 'three_prime_UTR':
                    gtf_attrs.append('feature_type "3UTR"')
                else:
                    gtf_attrs.append('feature_type "UTR"')
            elif feature_type in ['miRNA', 'tRNA']:
                gtf_attrs.append('feature_type "{}"'.format(feature_type))
            
            # Write GTF line
            gtf_line = '\t'.join([
                seqname,
                source,
                feature_type,
                str(start),
                str(end),
                score,
                strand,
                frame,
                '; '.join(gtf_attrs) + ';'
            ])
            fh_out.write(gtf_line + '\n')

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert GFF3 to GTF format following ENSEMBL conventions',
                                     usage='%(prog)s [-h] [-v,--version]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--gff3', action='store', dest='gff3', required=True, metavar='',
                        help='Input GFF3 file')
    parser.add_argument('-o', '--gtf', action='store', dest='gtf', required=True, metavar='',
                        help='Output GTF file')
    parser.add_argument('-m', '--chromosome_mapping', action='store', dest='chromosome_mapping_file', 
                        default=None, metavar='',
                        help='Path to chromosome mapping file (two-column TSV/CSV with tab delimiter: original_name<tab>new_name). '
                             'Example file content: chr1<tab>1, chrX<tab>X, chrM<tab>MT. If not specified, chromosome names are kept unchanged.')
    parser.add_argument('--source-genome', action='store', dest='source_genome_version', default=None, metavar='',
                        help='Source genome version (e.g., hg19, GRCh37). Required for liftover.')
    parser.add_argument('--target-genome', action='store', dest='target_genome_version', default=None, metavar='',
                        help='Target genome version (e.g., hg38, GRCh38). Required for liftover.')
    parser.add_argument('--chain-file', action='store', dest='chain_file', default=None, metavar='',
                        help='Path to chain file for coordinate liftover (required if source/target genome versions are specified). '
                             'Chain files can be downloaded from UCSC Genome Browser. Requires pyliftover: pip install pyliftover')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()
    
    # Validate liftover parameters
    if (args.source_genome_version or args.target_genome_version or args.chain_file):
        if not (args.source_genome_version and args.target_genome_version and args.chain_file):
            parser.error("--source-genome, --target-genome, and --chain-file must all be specified together for liftover")
    
    gff3_to_gtf(args.gff3, args.gtf, args.chromosome_mapping_file,
                args.source_genome_version, args.target_genome_version,
                args.chain_file)