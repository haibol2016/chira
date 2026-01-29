# ChiRA - Chimeric Read Analyzer

**Version**: 1.4.4 (Modified with performance optimizations)

ChiRA is a set of tools to analyze RNA-RNA interactome experimental data such as CLASH, CLEAR-CLIP, PARIS, SPLASH etc. Following are the descriptions of the each tool. Here we provide descriptions about the input and ouptput files. For the detailed description of the other parameters please look at the help texts of tools.

**Note**: This is a modified version of ChiRA (based on v1.4.3) with significant performance optimizations and new features. The original code is licensed under GPL v3, and this modified version maintains the same license. All changes are documented in the "Recent Improvements" section below and in [CHANGELOG.md](CHANGELOG.md).

## Version History
- **v1.4.4** (Current, 2026-01-26): Performance optimizations, BEDTools compatibility, sample name support, comprehensive bug fixes, and new utility scripts
- **v1.4.3** (Previous): Original version from GitHub (https://github.com/pavanvidem/chira)

See [CHANGELOG.md](CHANGELOG.md) for detailed change history.

## Recent Improvements (v1.4.4)

This version includes significant performance optimizations and new features. For detailed change history, see [CHANGELOG.md](CHANGELOG.md).

### Key Highlights

**Performance Improvements:**
- **3-10x faster** overall processing, especially in `chira_quantify.py`
- **2-5x faster** FASTQ parsing in `chira_collapse.py` (removed Biopython dependency)
- **2-5x faster** CIGAR string parsing with pre-compiled regex patterns
- **20-40% faster** BAM file processing in `chira_map.py`
- **10-50x faster** EM algorithm with optimized memory operations

**New Features:**
- Added `--sample_name` parameter to `chira_extract.py` for customizable output file names
- Output files now use format: `{sample_name}.chimeras.txt`, `{sample_name}.singletons.txt`, `{sample_name}.interactions.txt`
- Added `--gzip` option to `chira_extract.py` for compressing large output files (saves disk space, faster I/O for large files)
- Added header rows to interactions output file with column descriptions
- Added `mirna_position` column to chimeras output indicating read orientation (miRNA_first or miRNA_last)
- Automatic BEDTools version detection - works with both old and new command formats
- New utility scripts for reference file preparation (see Utility Scripts section)
- Docker support with pre-installed dependencies (see Docker Support section)
- Singularity/Apptainer support - see [SINGULARITY_SETUP.md](SINGULARITY_SETUP.md) for detailed instructions

**Compatibility:**
- BEDTools: Automatically supports both `intersectBed`/`fastaFromBed` (old) and `bedtools intersect`/`bedtools getfasta` (new)
- No code changes needed when switching between BEDTools versions

**Bug Fixes:**
- **chira_map.py**: Fixed first iteration bug that wrote empty string to unmapped FASTA file
- **chira_merge.py**: Added zero-length match checks to prevent division by zero errors
- **chira_extract.py**: Fixed multiple bugs including undefined variables, logic errors in strand/chromosome checks, index out-of-bounds in TPM cutoff, and string slicing issues
- **chira_utilities.py**: Fixed `median()` function bug where input wasn't sorted; improved `get_bedtools_command()` to verify command success
- **chira_quantify.py**: Fixed CRL iteration bug that skipped index 0
- Improved file handling with proper context managers throughout
- Fixed BEDTools command compatibility across versions with automatic detection

For complete details with line-by-line changes, please refer to [CHANGELOG.md](CHANGELOG.md).

## Installation

### Recommended: Container Installation

**The easiest way to use ChiRA is with the provided container images**, which include all dependencies pre-installed:

**Docker:**
```bash
# Run ChiRA commands
docker run --rm -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  docker.io/nemat1976/chiraplus:v0.0.1 chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

**Singularity/Apptainer (for HPC systems):**
```bash
# Pull image
singularity pull docker://docker.io/nemat1976/chiraplus:v0.0.1

# Run ChiRA commands (entrypoint script handles environment automatically)
singularity exec -B $(pwd)/data:/app/data -B $(pwd)/output:/app/output \
  chira_latest.sif chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

For detailed Singularity/Apptainer setup and usage instructions, see [SINGULARITY_SETUP.md](SINGULARITY_SETUP.md).

### Docker Image Contents

The Docker image includes:
- **Python packages**: biopython, bcbiogff, pysam, requests, pyliftover
- **Bioinformatics tools**: bwa, samtools, bedtools, intarna
- **All ChiRA scripts and utilities**: Pre-installed and executable
- **Environment setup**: Proper PATH and PYTHONPATH configuration

### Running ChiRA in Docker

**Basic usage:**
```bash
docker run --rm -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  docker.io/nemat1976/chiraplus:v0.0.1  chira_collapse.py -i data/input.fastq -o output/collapsed.fasta
```

**Interactive shell:**
```bash
docker run --rm -it -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  docker.io/nemat1976/chiraplus:v0.0.1 bash
```

**Volume mounts:**
```bash
docker run --rm \
  -v /path/to/data:/app/data \
  -v /path/to/output:/app/output \
  -v /path/to/references:/app/references \
  docker.io/nemat1976/chiraplus:v0.0.1 chira_extract.py [options]
```

### Manual Installation

If you prefer to install dependencies manually:

**Core Python packages:**
- biopython
- bcbiogff
- pysam

**Optional Python packages:**
- pyliftover (for `download_mirbase_gff3.py` coordinate liftover)
- requests (for `download_ensembl.py`)

**Command-line tools:**
- bwa (recommended for alignment)
- samtools
- bedtools
- intarna (optional, for hybridization prediction)

**Installation commands:**
```bash
pip install biopython bcbiogff pysam requests pyliftover
conda install -c bioconda bwa samtools bedtools intarna
```

For a complete list of dependencies, see [DEPENDENCIES.md](DEPENDENCIES.md).

---

## Workflow

The ChiRA pipeline consists of five main steps, from raw FASTQ files to final interaction results. The following workflow shows a typical analysis using split-reference mapping (miRNA and target transcriptomes separately).

### Step 1: Prepare Reference Files

Before starting the analysis, prepare your reference files:

**1.1 Download mature miRNA sequences:**
```bash
# Download species-specific mature miRNAs from miRBase
download_mirbase_mature.py -s hsa -o mature_mirna_hsa.fasta
```

**1.2 Download Ensembl reference files:**
```bash
# Download cDNA, ncRNA, GTF, and genome FASTA from Ensembl
download_ensembl.py -s homo_sapiens -g 110 -t 110 -o ./ensembl_files
```

**1.3 Download miRBase GFF3 file:**
```bash
# Download species-specific GFF3 file from miRBase (contains mature miRNA coordinates)
download_mirbase_gff3.py -s hsa -o hsa.gff3

# With chromosome name mapping (if needed)
download_mirbase_gff3.py -s hsa -o hsa.gff3 -m chr_mapping.txt

# Download specific version (e.g., version 21)
download_mirbase_gff3.py -s hsa -o hsa_v21.gff3 --mirbase-version 21
```

**Note:** The GFF3 file from miRBase contains mature miRNA coordinates and can be used directly with ChiRA. You don't need to convert it to GTF format. 

**1.4 Prepare target transcriptome (remove miRNAs):**
```bash
# Remove miRNA entries from Ensembl GTF
remove_mirna_hairpin_from_gtf.py -i ensembl_files/Homo_sapiens.GRCh38.110.gtf \
  -o target_transcriptome.gtf

# Remove miRNA sequences from Ensembl cDNA FASTA
remove_mirna_hairpin_from_fasta.py -g target_transcriptome.gtf \
  -f ensembl_files/Homo_sapiens.GRCh38.cdna.all.fa \
  -o target_transcriptome.fasta
```

**1.6 Combine miRNA and target GTF (for annotation):**
```bash
# Concatenate miRNA GTF with target GTF
concatenate_gtf.py -m mature_mirna.gtf -t target_transcriptome.gtf \
  -o combined_annotation.gtf
```

**1.7 Combine cDNA and ncRNA FASTA files (optional):**
```bash
# Concatenate cDNA and ncRNA FASTA files (removes version numbers from transcript IDs)
concatenate_fasta.py -c cdna.fasta -n ncrna.fasta -o combined_transcriptome.fasta
```

**Result:** You now have:
- `ref1.fasta`: Mature miRNA sequences (from miRBase)
- `ref2.fasta`: Target transcriptome sequences (without miRNAs)
- `combined_annotation.gtf`: Combined annotation file (if combining miRNA and target annotations)

---

### Step 2: Collapse Reads

Deduplicate reads from FASTQ file:

```bash
chira_collapse.py -i raw_reads.fastq -o collapsed_reads.fasta -u 12
```

**Input:** Raw FASTQ file (quality and adapter trimmed)  
**Output:** FASTA file with unique sequences and read counts

---

### Step 3: Map Reads

Map collapsed reads to reference transcriptomes:

```bash
# Build indices and map (one command)
chira_map.py -i collapsed_reads.fasta -o mapping_output \
  -f1 ref1.fasta -f2 ref2.fasta -b -a bwa -s both -p 8

# Or use pre-built indices
chira_map.py -i collapsed_reads.fasta -o mapping_output \
  -x1 index1 -x2 index2 -a bwa -s both -p 8
```

**Input:** Collapsed FASTA file  
**Output:** BAM files and `mapped.bed` file with all alignments

---

### Step 4: Merge Alignments

Merge overlapping alignments into loci:

```bash
chira_merge.py -b mapping_output/mapped.bed -o merge_output \
  -g combined_annotation.gtf -f1 ref1.fasta -f2 ref2.fasta \
  -ao 0.7 -so 0.7
```

**Input:** `mapped.bed` from Step 3  
**Output:** `loci.txt` file with merged alignments

---

### Step 5: Quantify CRLs

Build Chimeric Read Loci (CRLs) and quantify:

```bash
chira_quantify.py -b merge_output/segments.bed \
  -m merge_output/loci.txt -o quantify_output \
  -cs 0.7 -ls 10
```

**Input:** `segments.bed` and `loci.txt` from Step 4  
**Output:** `loci.txt` with TPM values and CRL assignments

---

### Step 6: Extract Interactions

Extract chimeric reads and summarize interactions:

```bash
chira_extract.py -l quantify_output/loci.txt -o extract_output \
  -f1 ref1.fasta -f2 ref2.fasta -n sample1 \
  -g combined_annotation.gtf -r -s -tc 0.1 -sc 0.5
```

**Input:** `loci.txt` from Step 5  
**Output:** 
- `sample1.chimeras.txt`: Individual chimeric reads
- `sample1.singletons.txt`: Non-chimeric reads
- `sample1.interactions.txt`: Summarized interactions (if `-s` used)

---

### Complete Workflow Example

Here's a complete workflow script for a typical analysis:

```bash
#!/bin/bash

# Set variables
SAMPLE="sample1"
SPECIES="hsa"  # human
ENSEMBL_VERSION="110"

# Step 1: Prepare references
echo "Step 1: Preparing reference files..."

# Download mature miRNA sequences (FASTA)
download_mirbase_mature.py -s $SPECIES -o ref1.fasta

# Download miRBase GFF3 (contains mature miRNA coordinates, can be used directly)
download_mirbase_gff3.py -s $SPECIES -o mirbase.gff3

download_ensembl.py -s homo_sapiens -g $ENSEMBL_VERSION -t $ENSEMBL_VERSION -o ./ensembl


# Remove miRNAs from target transcriptome
remove_mirna_hairpin_from_gtf.py -i ensembl/Homo_sapiens.GRCh38.$ENSEMBL_VERSION.gtf \
  -o target_transcriptome.gtf
remove_mirna_hairpin_from_fasta.py -g target_transcriptome.gtf \
  -f ensembl/Homo_sapiens.GRCh38.cdna.all.fa -o ref2.fasta

# Combine annotations (if you have miRNA GTF)
# concatenate_gtf.py -m mature_mirna.gtf -t target_transcriptome.gtf -o combined.gtf

# Step 2: Collapse reads
echo "Step 2: Collapsing reads..."
chira_collapse.py -i raw_reads.fastq -o collapsed.fasta -u 12

# Step 3: Map reads
echo "Step 3: Mapping reads..."
chira_map.py -i collapsed.fasta -o mapping -f1 ref1.fasta -f2 ref2.fasta \
  -b -a bwa -s both -p 8

# Step 4: Merge alignments
echo "Step 4: Merging alignments..."
chira_merge.py -b mapping/mapped.bed -o merge -g target_transcriptome.gtf \
  -f1 ref1.fasta -f2 ref2.fasta -ao 0.7 -so 0.7

# Step 5: Quantify CRLs
echo "Step 5: Quantifying CRLs..."
chira_quantify.py -b merge/segments.bed -m merge/loci.txt -o quantify \
  -cs 0.7 -ls 10

# Step 6: Extract interactions
echo "Step 6: Extracting interactions..."
chira_extract.py -l quantify/loci.txt -o extract -f1 ref1.fasta -f2 ref2.fasta \
  -n $SAMPLE -g target_transcriptome.gtf -r -s -tc 0.1 -sc 0.5

echo "Analysis complete! Results in extract/"
```

---

### Docker Workflow

If using Docker, the workflow is similar but commands are prefixed with `docker run`:

```bash
# Build image once
docker build -t chira:latest .

# Run each step with volume mounts
docker run --rm -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  chira:latest python chira_collapse.py -i data/input.fastq -o output/collapsed.fasta

docker run --rm -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  chira:latest python chira_map.py -i output/collapsed.fasta -o output/mapping \
  -f1 data/ref1.fasta -f2 data/ref2.fasta -b -a bwa

# ... continue with remaining steps
```

---

## Tool Documentation

### chira_collapse.py

**Description:** Deduplicates reads from a FASTQ file by collapsing identical sequences and their UMIs (Unique Molecular Identifiers), if available, then outputs a FASTA file with unique sequences and their read counts. This is typically the first step in the ChiRA pipeline to reduce redundancy before mapping.

**Key Features:**
- Handles UMI-tagged reads (optional UMI trimming from 5' end)
- Counts occurrences of each unique sequence
- Fast raw file parsing (no Biopython dependency)

**Required Inputs:**
- `-i, --fastq`: Quality and adapter trimmed FASTQ file

**Optional Parameters:**
- `-u, --umi_len`: Length of UMI to trim from 5' end of each read (default: 0, no UMI)

**Outputs:**
- `-o, --fasta`: FASTA file with unique sequences
  - Header format: `>sequence_id|UMI|read_count` (or `>sequence_id|read_count` if no UMI)
  - Each unique sequence appears once with its total count

**Usage Example:**
```bash
chira_collapse.py -i input.fastq -o output.fasta -u 12
```

---

### chira_map.py

**Description:** Maps reads to a reference transcriptome using either BWA-MEM (recommended) or CLAN aligner. Performs two-pass mapping (long then short segments) to capture chimeric alignments. Supports strand-specific mapping and can handle split-reference genomes (e.g., miRNA and target transcriptomes separately).

**Key Features:**
- Two alignment algorithms: BWA-MEM (default, recommended for speed) or CLAN (slower, use only if needed)
- Two-pass mapping strategy: first pass for long segments, second pass for short/unmapped segments
- Strand-specificity options: forward, reverse-complement, or both
- Handles chimeric reads with configurable overlap between segments
- Automatic index building option

**Required Inputs:**
- `-i, --query_fasta`: FASTA file containing reads (typically from `chira_collapse.py`)
- `-o, --outdir`: Output directory for BAM and BED files
- Either `-x1, --index1` (pre-built index) or `-f1, --ref_fasta1` (reference FASTA to build index)

**Optional Inputs:**
- `-x2, --index2`: Second priority index file (for split-reference)
- `-f2, --ref_fasta2`: Second priority reference FASTA file
- `-b, --build`: Build indices from reference FASTA files

**Split Reference:**
A split reference uses two separate reference FASTA files instead of one combined file. This is useful for:
- **Separating different RNA types**: For example, use `ref_fasta1` for miRNA sequences and `ref_fasta2` for target transcript sequences
- **Better chimeric detection**: Helps identify chimeric reads that span between the two reference types (e.g., miRNA-target interactions)
- **Smaller indices**: Each reference can be indexed separately, potentially reducing memory usage

**How to prepare a split reference:**
1. **Create two separate FASTA files:**
   - `ref1.fasta`: First priority reference (e.g., mature miRNA sequences from [miRBase](https://www.mirbase.org/))
   - `ref2.fasta`: Second priority reference (e.g., target transcript sequences)
   
   **Important:** When preparing `ref2.fasta` (target transcripts):
   - **Remove mature miRNA sequences** from the target transcript reference
   - **Remove miRNA hairpin sequences** from the target transcript reference
   - This prevents false-positive chimeric alignments where miRNA sequences align to both references

2. **Prepare GTF/GFF annotation file:**
   - **Include mature miRNA annotations** (e.g., `3p_mature_mir`, `5p_mature_mir`, `mature_mir`) download .gff3 from miRBase Downloads: https://www.mirbase.org/download/
   - **Exclude miRNA hairpin annotations** to avoid confusion
   - This ensures proper identification of miRNA vs target loci in downstream analysis

3. **Option A - Build indices automatically:**
   ```bash
   chira_map.py -i reads.fasta -o output_dir -f1 ref1.fasta -f2 ref2.fasta -b -a bwa
   ```
   The `-b` flag will automatically build indices for both references.

4. **Option B - Pre-build indices (recommended for repeated use):**
   ```bash
   # Build index1
   bwa index -p index1 ref1.fasta
   
   # Build index2
   bwa index -p index2 ref2.fasta
   
   # Then use pre-built indices
   chira_map.py -i reads.fasta -o output_dir -x1 index1 -x2 index2 -a bwa
   ```

**Note:** When using split reference, reads are mapped to both references separately, then the results are merged. The tool tracks which reference each alignment came from, which is important for downstream analysis in `chira_merge.py` and `chira_extract.py`.

**Key Parameters:**
- `-a, --aligner`: Alignment program (`bwa` or `clan`, default: `bwa`)
  - **Note:** BWA-MEM is significantly faster than CLAN and is recommended for most use cases. CLAN should only be used if specific alignment characteristics are required.
- `-s, --stranded`: Strand specificity (`fw`, `rc`, or `both`, default: `fw`)
- `-l1, --seed_length1`: Seed length for 1st mapping iteration (default: 12)
- `-l2, --seed_length2`: Seed length for 2nd mapping iteration (default: 16)
- `-s1, --align_score1`: Minimum alignment score for 1st iteration (default: 18)
- `-s2, --align_score2`: Minimum alignment score for 2nd iteration (default: 16)
- `-co, --chimeric_overlap`: Max bases between chimeric segments (default: 2)
- `-p, --processes`: Number of parallel processes (default: 1)

**Outputs:**
- `sorted.bam`: Sorted BAM file
- `sorted.bed`: BED file containing all alignments
- `unmapped.fasta`: FASTA file with unmapped reads (optional)

**Usage Example:**
```bash
chira_map.py -i reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both -p 8
```

---

### chira_merge.py

**Description:** Merges overlapping aligned positions to define read-concentrated loci (RCLs). If a GTF annotation file is provided, transcriptomic coordinates are converted to genomic coordinates. Segments reads into aligned portions and merges overlapping segments using configurable overlap thresholds.

**Key Features:**
- Converts transcriptomic to genomic coordinates (if GTF provided)
- Segments reads based on alignment parts
- Merges overlapping alignments into loci
- Two merging modes: segment-based or block-based (using Blockbuster)
- Filters by minimum locus size

**Required Inputs:**
- `-b, --bed`: Input BED file containing alignments (from `chira_map.py`)
- `-o, --outdir`: Output directory

**Optional Inputs:**
- `-g, --gtf`: Annotation GTF/GFF file for coordinate conversion
- `-f1, --ref_fasta1`: First priority reference FASTA file
- `-f2, --ref_fasta2`: Second priority reference FASTA file

**Key Parameters:**
- `-so, --segment_overlap`: Minimum overlap fraction to merge segments (default: 0.7)
- `-ao, --alignment_overlap`: Minimum overlap fraction to merge alignments (default: 0.7)
- `-lt, --length_threshold`: Minimum alignment length as fraction of longest (default: 0.9)
- `-co, --chimeric_overlap`: Max bases between chimeric segments (default: 2)
- `-c, --chimeric_only`: Consider only chimeric reads for merging
- `-ls, --min_locus_size`: Minimum alignments per merged locus (default: 1)
- `-bb, --block_based`: Use Blockbuster for block-based merging
- Blockbuster parameters: `-d, --distance`, `-mc, --min_cluster_height`, `-mb, --min_block_height`, `-sc, --scale`

**Outputs:**
- `segments.bed`: A BED file with reads categorized into segments
- `merged.bed`: A tabular file with merged alignments
  - Column 4: All alignments merged into that location
  - Column 5: Number of reads supporting the locus

**Usage Example:**
```bash
chira_merge.py -b mapped.bed -o output_dir -g annotation.gtf -f1 ref1.fasta -ao 0.7 -so 0.7
```

---

### chira_quantify.py

**Description:** Creates Chimeric Read Loci (CRLs) from merged BED files and quantifies them using an Expectation-Maximization (EM) algorithm to handle multi-mapping reads. Calculates TPM (Transcripts Per Million) values for each CRL.

**Key Features:**
- Builds CRLs by grouping loci that share significant read overlap
- EM algorithm for resolving multi-mapping reads
- TPM normalization for expression quantification
- Configurable CRL building thresholds

**Required Inputs:**
- `-b, --bed`: Input BED file containing alignment segments (segments.bed from `chira_merge.py`)
- `-m, --merged_bed`: Input merged BED file (merged.bed from `chira_merge.py`)
- `-o, --outdir`: Output directory

**Key Parameters:**
- `-cs, --crl_share`: Minimum fraction of locus reads that must overlap with all CRL loci to merge (default: 0.7)
- `-ls, --min_locus_size`: Minimum reads per locus to participate in CRL creation (default: 10)
- `-e, --em_threshold`: EM algorithm convergence threshold (default: 0.00001)
- `-crl, --build_crls_too`: Create CRLs in addition to quantification

**Outputs:**
- `loci.counts`: Tabular file containing reads, their CRLs, and TPM values
  - Each line represents a read-CRL association with TPM quantification
  - Used as input for `chira_extract.py`

**Usage Example:**
```bash
chira_quantify.py -b segments.bed -m loci.txt -o output_dir -cs 0.7 -ls 10
```

---

### chira_extract.py

**Description:** Extracts the best chimeric alignments for each read and optionally performs RNA-RNA hybridization prediction using IntaRNA. Summarizes interactions at the locus level and identifies miRNA-target pairs. This is the final step that produces the interaction results.

**Key Features:**
- Extracts chimeric reads with best scoring alignments
- Optional RNA-RNA hybridization prediction (IntaRNA)
- Interaction summarization at locus level
- Identifies miRNA vs target loci based on annotation
- Customizable output file names with sample name prefix

**Required Inputs:**
- `-l, --loci`: Tabular file containing CRL information ('loci.counts' from `chira_quantify.py`)
- `-o, --out`: Output directory path
- `-f1, --ref_fasta1`: First priority reference FASTA file
- `-n, --sample_name`: Sample name prefix for output files (required)

**Optional Inputs:**
- `-g, --gtf`: Annotation GTF/GFF file for locus annotation including miRBase gff3 for mature miRNA.
- `-f2, --ref_fasta2`: Second priority reference FASTA file
- `-f, --ref`: Reference genomic FASTA file 

**Key Parameters:**
- `-r, --hybridize`: Enable RNA-RNA hybridization prediction (IntaRNA)
- `-tc, --tpm_cutoff`: TPM percentile cutoff for filtering (default: 0)
- `-sc, --score_cutoff`: Hybridization score cutoff (default: 0.0)
- `-co, --chimeric_overlap`: Max bases between chimeric segments (default: 2)
- `-ns, --no_seed`: Do not enforce seed interactions
- `-m, --intarna_mode`: IntaRNA mode (`H`=heuristic, `M`=exact, `S`=seed-only, default: `H`)
- `-t, --temperature`: Temperature in Celsius for energy parameters (default: 37)
- `-sbp, --seed_bp`: Number of base pairs in seed region (default: 5)
- `-acc, --accessibility`: Compute accessibility (`C` or `N`, default: `N`)
- `-p, --processes`: Number of parallel processes (default: 1)
- `-s, --summarize`: Summarize interactions at locus level
- `-z, --gzip`: Compress output files (chimeras and singletons) with gzip (optional)

**Outputs:**

**Note:** If `--gzip` is specified, output files will have `.gz` extension (e.g., `{sample_name}.chimeras.txt.gz`, `{sample_name}.singletons.txt.gz`). Compression is applied only to final merged files, not intermediate files, for optimal performance.

**1. `{sample_name}.chimeras.txt`** (or `.txt.gz` if `--gzip` is used) - Tabular file with chimeric read information in an extended BED format

Header columns (34 total):
- `read_id`: Read identifier (from collapsed FASTQ)
- `transcript_id_1`, `transcript_id_2`: Transcript IDs for locus 1 and locus 2
- `gene_id_1`, `gene_id_2`: Gene IDs for locus 1 and locus 2
- `gene_symbol_1`, `gene_symbol_2`: Gene symbols for locus 1 and locus 2
- `annotation_region_1`, `annotation_region_2`: Annotation regions (e.g., `3p_mature_mir`, `5p_mature_mir`, `mature_mir` for miRNA; gene/exon types for targets)
- `transcript_start_1`, `transcript_end_1`, `transcript_strand_1`: Transcriptomic coordinates for locus 1
- `transcript_length_1`: Alignment length for locus 1
- `transcript_start_2`, `transcript_end_2`, `transcript_strand_2`: Transcriptomic coordinates for locus 2
- `transcript_length_2`: Alignment length for locus 2
- `read_alignment_info`: Read alignment information (format: `arm1_start,arm1_end,arm2_start,arm2_end,read_length`)
  - Contains positions of both alignment arms within the read sequence
  - **Use to determine actual read orientation**: Compare `arm1_start` vs `arm2_start`
    - If `arm1_start < arm2_start`: Read is 5' locus1 → locus2 3'
    - If `arm2_start < arm1_start`: Read is 5' locus2 → locus1 3' (reoriented in output for split reference)
- `genomic_coordinates_1`, `genomic_coordinates_2`: Genomic coordinates (if GTF provided)
- `locus_id_1`, `locus_id_2`: Locus identifiers in format `chr:start:end:strand`
- `crl_group_id_1`, `crl_group_id_2`: CRL group IDs
- `tpm_1`, `tpm_2`: TPM (Transcripts Per Million) values for each locus
- `alignment_score_1`, `alignment_score_2`: Alignment scores for each locus
- `combined_alignment_score`: Combined score (score1 × score2)
- `hybridized_sequences`: Hybridized sequences (format: `sequence1&sequence2`, or `NA` if not hybridized)
- `hybridization_structure`: Dot-bracket notation for RNA-RNA hybridization structure (or `NA`)
- `hybridization_positions`: Hybridization position information (or `NA`)
- `hybridization_mfe_kcal_mol`: Minimum free energy of hybridization in kcal/mol (or `NA`)
- `mirna_read_position`: Indicates whether miRNA is at the 5' or 3' end of the read
  - `miRNA_first`: miRNA is at the 5' end (5' miRNA → target 3')
  - `miRNA_last`: miRNA is at the 3' end (5' target → miRNA 3')
  - `NA`: Neither locus is annotated as miRNA (rare, may occur in non-split reference scenarios)

**Note on chimeric read orientation:**
The tool detects chimeric reads in both orientations:
- **5' miRNA → target 3'**: miRNA at 5' end, target at 3' end
- **5' target → miRNA 3'**: Target at 5' end, miRNA at 3' end

For **split reference** (when `-f2` is provided), the output is standardized so that:
- `locus1` always corresponds to the first reference (typically miRNA from `ref_fasta1`)
- `locus2` always corresponds to the second reference (typically target from `ref_fasta2`)

**To determine the actual read orientation**, check the `read_info` field:
- The `read_info` field contains: `arm1_start,arm1_end,arm2_start,arm2_end,read_length`
- Compare the start positions:
  - If `arm1_start < arm2_start`: The read is in **5' locus1 → locus2 3'** orientation
  - If `arm2_start < arm1_start`: The read is in **5' locus2 → locus1 3'** orientation (reoriented in output for split reference)

**Example:**
- If `read_info = "1,20,21,40,50"`: arm1 starts at position 1, arm2 starts at position 21 → Read is 5' locus1 → locus2 3'
- If `read_info = "21,40,1,20,50"`: arm2 starts at position 1, arm1 starts at position 21 → Read is 5' locus2 → locus1 3' (original orientation was target → miRNA, but output shows miRNA → target)

In the **interactions file**, both orientations are merged into a single entry to avoid duplicates, but the `read_info` field in the chimeras file preserves the actual orientation information.

**2. `{sample_name}.singletons.txt`** (or `.txt.gz` if `--gzip` is used) - Tabular file with singleton reads (non-chimeric alignments)

Header columns (14 total):
- `read_id`: Read identifier
- `transcript_id`: Transcript ID
- `gene_id`: Gene ID
- `gene_symbol`: Gene symbol
- `annotation_region`: Annotation region type
- `transcript_start`, `transcript_end`, `transcript_strand`: Transcriptomic coordinates
- `transcript_length`: Alignment length
- `read_alignment_info`: Read alignment information
- `genomic_coordinates`: Genomic coordinates (if GTF provided)
- `locus_id`: Locus identifier in format `chr:start:end:strand`
- `crl_group_id`: CRL group ID
- `tpm`: TPM value
- `alignment_score`: Alignment score

**3. `{sample_name}.interactions.txt`** - Tabular file with detected interactions (if `--summarize` used). This file is always uncompressed for compatibility with downstream analysis tools.

Header columns (24 total):
- `supporting_read_count`: Number of reads supporting this interaction
- `locus_1_chromosome`, `locus_1_start`, `locus_1_end`, `locus_1_strand`: Genomic coordinates for locus 1
- `locus_2_chromosome`, `locus_2_start`, `locus_2_end`, `locus_2_strand`: Genomic coordinates for locus 2
- `locus_1_sequence`, `locus_2_sequence`: Sequences from each locus (or `NA` if not hybridized)
- `hybridization_structure_dotbracket`: RNA-RNA hybridization structure in dot-bracket notation (format: `structure1&structure2`, or `NA`)
- `hybridization_mfe_kcal_mol`: Minimum free energy of hybridization in kcal/mol (or `NA`)
- `hybridized_sequence_segments`: Hybridized sequence segments (format: `seq1\tseq2`, or `NA\tNA`)
- `hybridization_start_positions`: Start positions of hybridization within sequences (format: `pos1&pos2`, or `NA`)
- `hybridization_genomic_coordinates`: Genomic coordinates of hybridization region (format: `chr1:start1:end1:strand1\tchr2:start2:end2:strand2`, or `NA`)
- `tpm_locus_1`, `tpm_locus_2`: TPM values for each locus
- `tpm_combined`: Combined TPM (tpm1 + tpm2)
- `alignment_score_locus_1`, `alignment_score_locus_2`: Alignment scores for each locus
- `combined_alignment_score`: Combined score (score1 × score2)
- `annotation_region_locus_1`, `annotation_region_locus_2`: Annotation regions (semicolon-separated if multiple; use to identify miRNA vs target)
- `reference_transcript_id_1`, `reference_transcript_id_2`: Reference sequence IDs (semicolon-separated if multiple)

**Note:** The interactions file includes comment lines explaining how to identify miRNA vs target loci. Check the `annotation_region_locus_1` and `annotation_region_locus_2` fields - miRNA annotations typically include: `miRNA`, `3p_mature_mir`, `5p_mature_mir`, `mature_mir`.

**Usage Example:**
```bash
# Basic usage
chira_extract.py -l loci.txt -o output_dir -f1 ref1.fasta -n sample1 \
  -g annotation.gtf -r -s -tc 0.1 -sc 0.5

# With gzip compression (recommended for large files)
chira_extract.py -l loci.txt -o output_dir -f1 ref1.fasta -n sample1 \
  -g annotation.gtf -r -s -tc 0.1 -sc 0.5 --gzip
```

**Note:** The interactions file includes comments explaining how to identify miRNA vs target loci based on the `annotation_region_locus_1` and `annotation_region_locus_2` fields. miRNA annotations typically include: `miRNA`, `3p_mature_mir`, `5p_mature_mir`, `mature_mir`.

---

## Utility Scripts

The following utility scripts are provided to help prepare reference files and annotations for ChiRA analysis:


### download_ensembl.py

**Description:** Downloads cDNA, ncRNA, GTF, and genome FASTA files from Ensembl for a given species and release versions.

**Key Features:**
- Downloads primary assembly genome FASTA (not toplevel)
- Auto-detects assembly names or accepts explicit assembly parameter
- Supports HTTP and FTP with automatic fallback
- Automatically decompresses gzipped files

**Required Inputs:**
- `-s, --species`: Species name (e.g., homo_sapiens, mus_musculus, bos_taurus)
- `-g, --genome-version`: Ensembl release version for genome/cDNA/ncRNA (e.g., 110)
- `-t, --gtf-version`: Ensembl release version for GTF annotation (e.g., 110)
- `-o, --output-dir`: Output directory for downloaded files

**Optional Parameters:**
- `-a, --assembly`: Genome assembly name (e.g., GRCh38, GRCm39). Auto-detected if not specified.
- `--keep-compressed`: Keep compressed files after decompression
- `--no-decompress`: Do not decompress gzipped files
- `--timeout`: Download timeout in seconds (default: 60)

**Usage Example:**
```bash
download_ensembl.py -s homo_sapiens -g 110 -t 110 -o ./ensembl_files
```

---

### download_mirbase_mature.py

**Description:** Downloads species-specific mature miRNA sequences from miRBase.

**Key Features:**
- Downloads from specific miRBase version or CURRENT
- Extracts sequences by species code
- Handles gzipped files automatically

**Required Inputs:**
- `-s, --species`: Species code (e.g., hsa=human, mmu=mouse, bta=bovine, rno=rat)
- `-o, --output`: Output FASTA file path

**Optional Parameters:**
- `--mirbase-version`: miRBase version (e.g., "22.1"). Default: CURRENT
- `--keep-full`: Keep the full downloaded file after extraction
- `--timeout`: Download timeout in seconds (default: 30)

**Usage Example:**
```bash
download_mirbase_mature.py -s hsa -o mature_mirna_hsa.fasta
```

---

### download_mirbase_gff3.py

**Description:** Downloads species-specific GFF3 files from miRBase containing chromosomal coordinates of microRNAs. Supports coordinate liftover between genome assemblies and chromosome name mapping.

**Key Features:**
- Downloads current or version-specific GFF3 files
- Contains both miRNA_primary_transcript (hairpin precursors) and miRNA (mature sequences)
- **Coordinate liftover**: Convert coordinates between genome assemblies (e.g., hg19 → hg38) using pyliftover
- **Chromosome name mapping**: Rename chromosomes based on a mapping file
- Processing order: Download → Liftover → Rename chromosomes
- Can be used directly with ChiRA (no GTF conversion needed)

**Required Inputs:**
- `-s, --species`: Species code (e.g., hsa=human, mmu=mouse, bta=bovine, rno=rat)
- `-o, --output`: Output GFF3 file path

**Optional Parameters:**
- `--mirbase-version`: miRBase version number (e.g., "21"). Default: CURRENT version
- `--timeout`: Download timeout in seconds (default: 60)
- `--source-genome`: Source genome assembly name for coordinate liftover (e.g., hg19, hg38, mm9, mm10)
  - Required if `--chain-file` is provided
- `--target-genome`: Target genome assembly name for coordinate liftover (e.g., hg19, hg38, mm9, mm10)
  - Required if `--chain-file` is provided
- `--chain-file`: Path to chain file for coordinate liftover
  - Chain files can be downloaded from UCSC Genome Browser (e.g., https://hgdownload.soe.ucsc.edu/downloads.html)
  - Requires `pyliftover` package (install with: `pip install pyliftover`)
  - If provided, `--source-genome` and `--target-genome` must also be provided
- `-m, --chromosome-mapping`: Tab-separated file with two columns: `gff3_chromosome_name<tab>target_chromosome_name`
  - Chromosome names in the GFF3 file will be converted to target names in the output
  - Applied after coordinate liftover (if performed)

**Usage Examples:**
```bash
# Download current version
download_mirbase_gff3.py -s hsa -o hsa.gff3

# Download specific version with chromosome mapping
download_mirbase_gff3.py -s hsa -o hsa.gff3 --mirbase-version 21 -m chr_mapping.txt

# Download with coordinate liftover (hg19 to hg38)
download_mirbase_gff3.py -s hsa -o hsa_hg38.gff3 \
  --source-genome hg19 --target-genome hg38 \
  --chain-file hg19ToHg38.over.chain

# Download with both liftover and chromosome mapping
download_mirbase_gff3.py -s hsa -o hsa_processed.gff3 \
  --source-genome hg19 --target-genome hg38 \
  --chain-file hg19ToHg38.over.chain \
  -m chr_mapping.txt
```

**Note:** The GFF3 file from miRBase contains mature miRNA coordinates and can be used directly with ChiRA. You don't need to convert it to GTF format unless you want to combine it with Ensembl annotations.

**Coordinate Liftover:**
- Liftover converts coordinates from one genome assembly to another (e.g., GRCh37/hg19 to GRCh38/hg38)
- Chain files are available from UCSC Genome Browser for common assembly conversions
- Liftover is performed first, then chromosome renaming (if provided)
- Features that cannot be lifted over will retain their original coordinates

---

### concatenate_fasta.py

**Description:** Concatenates cDNA and ncRNA FASTA files into a single output file, removing version numbers from transcript IDs in the headers.

**Key Features:**
- Combines cDNA and ncRNA FASTA files
- Removes version numbers from transcript IDs (e.g., `ENST00000123456.1` → `ENST00000123456`)
- Supports both Ensembl style (space-separated) and GENCODE style (pipe-separated) headers
- At least one input file (cDNA or ncRNA) must be provided

**Required Inputs:**
- `-o, --output`: Output concatenated FASTA file

**Optional Parameters:**
- `-c, --cdna`: Input cDNA FASTA file (optional, but at least one input must be provided)
- `-n, --ncrna`: Input ncRNA FASTA file (optional, but at least one input must be provided)

**Usage Example:**
```bash
# Concatenate both cDNA and ncRNA
concatenate_fasta.py -c cdna.fasta -n ncrna.fasta -o combined_transcriptome.fasta

# Concatenate only cDNA
concatenate_fasta.py -c cdna.fasta -o combined.fasta
```

---

### remove_mirna_hairpin_from_gtf.py

**Description:** Removes microRNA entries from an Ensembl GTF file. Used for preparing target-only transcriptome annotations.

**Key Features:**
- Identifies miRNA entries by feature type, biotype, and optional regex pattern
- Flexible pattern matching for custom miRNA identification
- Preserves comment lines (optional removal)

**Required Inputs:**
- `-i, --input`: Input GTF file
- `-o, --output`: Output GTF file (without microRNA entries)

**Optional Parameters:**
- `-p, --pattern`: Regular expression pattern for matching miRNA in GTF attributes
  - If not provided, only feature type and biotype fields are checked
  - Example: `'gene_name\s+"[^"]*(?:[Mm][Ii][Rr][_-]|[^"]*[-_][Mm][Ii][Rr][_-])'`
- `--remove-comments`: Remove comment lines from output (default: keep comments)

**Usage Example:**
```bash
remove_mirna_hairpin_from_gtf.py -i annotation.gtf -o annotation_no_mirna.gtf
```

---

### remove_mirna_hairpin_from_fasta.py

**Description:** Removes miRNA FASTA records from a transcriptome FASTA file using transcript IDs extracted from a GTF file.

**Key Features:**
- Uses transcript IDs from GTF to identify miRNA sequences
- Supports various FASTA header formats
- Reuses miRNA detection logic from `remove_mirna_hairpin_from_gtf.py`

**Required Inputs:**
- `-g, --gtf`: Input GTF file (Ensembl format)
- `-f, --fasta`: Input transcriptome FASTA file
- `-o, --output`: Output FASTA file (without microRNA sequences)

**Optional Parameters:**
- `-p, --pattern`: Regular expression pattern for matching miRNA in GTF attributes (same as `remove_mirna_hairpin_from_gtf.py`)

**Usage Example:**
```bash
remove_mirna_hairpin_from_fasta.py -g annotation.gtf -f transcriptome.fasta -o transcriptome_no_mirna.fasta
```

---

### concatenate_gtf.py

**Description:** Concatenates mature miRNA GTF file with target transcriptome GTF file, removing comment lines from the miRNA GTF.

**Key Features:**
- Combines miRNA and target GTF files for split-reference analysis
- Removes comment lines from miRNA GTF
- Optionally removes comment lines from target GTF

**Required Inputs:**
- `-m, --mirna-gtf`: Mature miRNA GTF file (can be converted from miRBase GFF3 if needed)
- `-t, --target-gtf`: Target transcriptome GTF file (output from `remove_mirna_hairpin_from_gtf.py`)
- `-o, --output`: Output combined GTF file

**Optional Parameters:**
- `--remove-target-comments`: Remove comment lines from target GTF as well

**Usage Example:**
```bash
concatenate_gtf.py -m mature_mirna.gtf -t target_no_mirna.gtf -o combined.gtf
```

---
