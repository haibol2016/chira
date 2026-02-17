# ChiRA - Chimeric Read Analyzer

**Version**: 1.4.11 (Modified with performance optimizations and parallel computing support)

ChiRA is a set of tools to analyze RNA-RNA interactome experimental data such as CLASH, CLEAR-CLIP, PARIS, SPLASH etc. Following are the descriptions of the each tool. Here we provide descriptions about the input and ouptput files. For the detailed description of the other parameters please look at the help texts of tools.

**Note**: This is a modified version of ChiRA (based on v1.4.3) with significant performance optimizations and new features. The original code is licensed under GPL v3, and this modified version maintains the same license. All changes are documented in the "Recent Improvements" section below and in [CHANGELOG.md](CHANGELOG.md).

## Version History
- **v1.4.11** (Current, 2026-02-15): MPIRE made required dependency for optimal multiprocessing performance, removed fallback code, improved memory efficiency (50-90% reduction) and startup time (2-3x faster)
- **v1.4.10** (2026-02-15): Fixed batchtools submission issues (template path handling, JSON parsing), ensured all paths are absolute for cluster jobs, and refactored scripts for better code organization
- **v1.4.9** (2026-02-15): Added `--parallel_chunks` parameter for configurable chunk parallelism in chira_map.py
- **v1.4.8** (2026-02-15): Improved chunk-based parallelization for very large transcript counts (e.g., human genome with 387K+ transcripts), variable naming consistency improvements, and bug fixes in chira_merge.py
- **v1.4.7** (2026-02-15): New utility scripts (extract_transcripts_from_genome.py), U→T conversion in download_mirbase_mature.py, gffread support, and removal of deprecated scripts
- **v1.4.6** (2026-02-15): Multiprocessing improvements, adaptive I/O buffer sizing, FASTA chunking, I/O bottleneck detection, and code refactoring
- **v1.4.5** (2026-02-02): Parallel computing support, I/O optimizations, automatic memory management, and enhanced performance
- **v1.4.4** (2026-01-26): Performance optimizations, BEDTools compatibility, sample name support, comprehensive bug fixes, and new utility scripts
- **v1.4.3** (Previous): Original version from GitHub (https://github.com/pavanvidem/chira)

See [CHANGELOG.md](CHANGELOG.md) for detailed change history.

## Recent Improvements

### v1.4.7 (2026-02-15) - Utility Script Updates & Improvements

**New Features:**
- **extract_transcripts_from_genome.py**: New utility script to extract transcript FASTA sequences from genome FASTA using gffread
  - Replaces `remove_mirna_hairpin_from_fasta.py` with a more accurate approach
  - Extracts sequences directly from genome coordinates using gffread
- **download_mirbase_mature.py**: Added automatic U→T conversion in mature miRNA sequences for ChiRA compatibility
- **concatenate_gtf.py**: Updated documentation to reflect that miRBase GFF3 format can be used directly with ChiRA

**Docker:**
- Added `gffread` package to Dockerfile for transcript extraction functionality

**Removed:**
- `remove_mirna_hairpin_from_fasta.py`: Replaced by `extract_transcripts_from_genome.py`
- `concatenate_fasta.py`: No longer needed (use standard Unix tools like `cat`)

### v1.4.6 (2026-02-15) - Multiprocessing & I/O Optimizations

**Multiprocessing Improvements:**
- **chira_quantify.py**: Changed from `ThreadPoolExecutor` to `ProcessPoolExecutor` to MPIRE for EM algorithm (2-8x faster, bypasses Python GIL, 50-90% memory reduction)
  - Parameter: `-t, --threads` (use 0 for all available cores)
  - MPIRE is required (no fallback)
- **chira_extract.py**: Changed from `multiprocessing.Process` to MPIRE for chimera extraction (50-90% memory reduction, 2-3x faster startup)
  - MPIRE is required (no fallback)
- **chira_merge.py**: Changed from `ThreadPoolExecutor` to chunk-based `multiprocessing.Pool` for transcript processing (4-8x faster)
  - Parameter changed: `-t, --threads` → `-p, --processes` (default: None, auto-detects CPU count)
  - **Chunk-based strategy**: Groups transcripts into chunks (~1000 per chunk) to reduce overhead for very large datasets (e.g., 387K+ transcripts)
- **chira_map.py**: Enhanced with FASTA chunking and intelligent process management
  - New parameter: `--chunk_fasta` for splitting large FASTA files into chunks
  - **Process-aware chunking**: Parallel chunk execution automatically limited by available processes
  - **Total processes control**: `--processes` now specifies total processes, automatically divided among jobs/chunks
  - CPU usage guidance based on system capabilities

**I/O Performance Optimizations:**
- **Adaptive buffer sizing**: Automatically calculates optimal buffer size (8-16MB) based on available RAM
  - Reduces system calls by 1000-2000x for large files
  - **10-50x faster** I/O performance for large files (150M+ reads)
  - Uses `psutil` (optional) for intelligent buffer sizing
  - Falls back to 8MB default if `psutil` unavailable
- **Progress tracking**: Reports progress every 1M reads for large BAM files

**Code Improvements:**
- Refactored `chira_extract.py` and `chira_map.py` into modular functions for better maintainability
- Better error handling with `subprocess.run()` instead of `os.system()`
- Enhanced logging and user guidance

### v1.4.5 (2026-02-02) - Parallel Computing & Performance Enhancements

**Parallel Computing Support:**
- **Multi-threading** in `chira_quantify.py`: EM algorithm parallelization (2-4x faster)
  - New parameter: `-t, --threads` (use 0 for all available cores)
  - **Multi-threading** in `chira_merge.py`: Transcript processing parallelization (2-4x faster, improved in v1.4.8 with chunk-based strategy)
  - New parameter: `-t, --threads` (use 0 for all available cores)
- **Enhanced multi-threading** in `chira_map.py`: External tools (samtools, pysam, sort)
  - Automatic memory optimization with `psutil` (optional dependency)
  - New parameter: `--sort_memory` for manual memory specification
- **Parallel sort** support: Automatic detection of GNU sort `--parallel` (2-4x faster sorting)
- **I/O optimizations**: 2MB buffer sizes for faster file I/O (20-40% faster)

**Performance Improvements:**
- **2-4x faster** EM algorithm with multi-threading for large datasets
- **2-4x faster** transcript processing with multi-threading (improved in v1.4.8 with chunk-based strategy for very large datasets)
- **2-4x faster** BAM operations with multi-threaded tools
- **2-4x faster** file sorting with GNU sort parallel support
- **20-40% faster** I/O operations with optimized buffer sizes

### v1.4.4 (2026-01-26) - Performance Optimizations

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

**Pre-built Docker Image Available:**
A functional Docker image is available at: **`docker.io/nemat1976/chiraplus:v0.0.1`**

This image includes all dependencies and is ready to use. Simply pull and run:

**Docker:**
```bash
# Pull the image (first time only)
docker pull docker.io/nemat1976/chiraplus:v0.0.1

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
- **Python packages**: biopython, bcbiogff, pysam, requests, pyliftover, psutil
- **Bioinformatics tools**: bwa, samtools, bedtools, gffread, intarna
- **All ChiRA scripts and utilities**: Pre-installed and executable
- **Environment setup**: Proper PATH and PYTHONPATH configuration

**Note:** Optional tools (blockbuster, clan) are not included in the Docker image by default but can be added if needed for specific use cases.

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
- **psutil** (highly recommended for optimal performance)
  - **chira_map.py**: Automatic memory optimization for BAM sorting and I/O bottleneck detection
  - **chira_utilities.py**: Adaptive buffer sizing (8-16MB) for 10-50x I/O performance improvement
  - Install with: `pip install psutil` or `conda install psutil`
  - Falls back to safe defaults if not available (2GB per thread for BAM sorting, 8MB buffer for I/O)
- **mpire** (REQUIRED for `chira_quantify.py` and `chira_extract.py` parallel processing)
  - Enhanced multiprocessing framework for EM algorithm parallelization and chimera extraction
  - Benefits: 50-90% memory reduction, 2-3x faster startup, better performance
  - Install with: `pip install mpire` or `conda install -c conda-forge mpire`
  - Required for parallel processing (no fallback)
- pyliftover (for `download_mirbase_gff3.py` coordinate liftover)
- requests (for `download_ensembl.py`)

**Optional R packages (for batchtools HPC cluster support):**
- **batchtools** and **jsonlite** (required for `--use_batchtools` option in `chira_map.py`)
  - Install with: `conda install -c conda-forge r-batchtools r-jsonlite`
  - Or in R: `install.packages(c("batchtools", "jsonlite"))`
  - See [BATCHTOOLS_USAGE.md](BATCHTOOLS_USAGE.md) for detailed usage instructions

**Command-line tools:**
- bwa (recommended for alignment)
- samtools
- bedtools
- gffread (for `extract_transcripts_from_genome.py`)
  - Install with: `conda install -c bioconda gffread`
  - Part of GFF utilities from Johns Hopkins University
  - Documentation: https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread
- intarna (optional, for hybridization prediction)
- blockbuster (optional, for block-based merging in `chira_merge.py`)
- clan (optional, alternative aligner to BWA in `chira_map.py`)
- **GNU coreutils** (for parallel sort support)
  - **Linux**: Usually pre-installed (GNU sort is standard)
  - **macOS**: Install via Homebrew: `brew install coreutils`
  - Required for `--parallel` option in sort (version >= 8.6)
  - Code automatically detects and uses parallel sort when available

**Installation commands:**
```bash
# Core packages (required)
pip install biopython bcbiogff pysam mpire requests

# Optional packages (highly recommended for optimal performance)
pip install psutil  # For automatic memory optimization, I/O bottleneck detection, and adaptive buffer sizing (10-50x I/O improvement)
pip install pyliftover  # For coordinate liftover in download_mirbase_gff3.py

# Command-line tools
conda install -c bioconda bwa samtools bedtools gffread intarna

# Optional tools (for specific use cases)
conda install -c bioconda blockbuster  # For block-based merging in chira_merge.py
conda install -c bioconda clan  # Alternative aligner to BWA in chira_map.py

# GNU coreutils (for parallel sort on macOS)
# Linux: Usually pre-installed
# macOS: brew install coreutils
```

For a complete list of dependencies, see [DEPENDENCIES.md](DEPENDENCIES.md).

For batchtools HPC cluster usage, see [BATCHTOOLS_USAGE.md](BATCHTOOLS_USAGE.md).

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

# Extract transcript sequences from genome FASTA using filtered GTF
# This uses gffread to extract transcripts from the genome FASTA based on the filtered GTF
extract_transcripts_from_genome.py -g target_transcriptome.gtf \
  -f ensembl_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -o target_transcriptome.fasta
```

**1.5 Combine miRNA and target GTF (for annotation):**
```bash
# Concatenate miRNA GTF with target GTF
concatenate_gtf.py -m mature_mirna.gtf -t target_transcriptome.gtf \
  -o combined_annotation.gtf
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
# Build indices and map (one command, using 8 processes for parallelization)
chira_map.py -i collapsed_reads.fasta -o mapping_output \
  -f1 ref1.fasta -f2 ref2.fasta -b -a bwa -s both -p 8

# Or use pre-built indices (using 8 processes for parallelization)
chira_map.py -i collapsed_reads.fasta -o mapping_output \
  -x1 index1 -x2 index2 -a bwa -s both -p 8

# Use all available CPU cores automatically
chira_map.py -i collapsed_reads.fasta -o mapping_output \
  -f1 ref1.fasta -f2 ref2.fasta -b -a bwa -s both -p 0
```

**Input:** Collapsed FASTA file  
**Output:** BAM files and `mapped.bed` file with all alignments

---

### Step 4: Merge Alignments

Merge overlapping alignments into loci:

```bash
# Using 8 processes for parallel transcript processing
chira_merge.py -b mapping_output/mapped.bed -o merge_output \
  -g combined_annotation.gtf -f1 ref1.fasta -f2 ref2.fasta \
  -ao 0.7 -so 0.7 -p 8

# Auto-detect CPU count (default behavior)
# Automatically creates optimal number of chunks based on transcript count
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
# Using 8 threads for parallel EM algorithm
chira_quantify.py -b merge_output/segments.bed \
  -m merge_output/loci.txt -o quantify_output \
  -cs 0.7 -ls 10 -t 8

# Use all available CPU cores automatically
chira_quantify.py -b merge_output/segments.bed \
  -m merge_output/loci.txt -o quantify_output \
  -cs 0.7 -ls 10 -t 0
```

**Input:** `segments.bed` and `loci.txt` from Step 4  
**Output:** `loci.txt` with TPM values and CRL assignments

---

### Step 6: Extract Interactions

Extract chimeric reads and summarize interactions:

```bash
# Using 8 processes for parallel extraction and hybridization
chira_extract.py -l quantify_output/loci.txt -o extract_output \
  -f1 ref1.fasta -f2 ref2.fasta -n sample1 \
  -g combined_annotation.gtf -r -s -tc 0.1 -sc 0.5 -p 8

# Use all available CPU cores automatically
chira_extract.py -l quantify_output/loci.txt -o extract_output \
  -f1 ref1.fasta -f2 ref2.fasta -n sample1 \
  -g combined_annotation.gtf -r -s -tc 0.1 -sc 0.5 -p 0
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
extract_transcripts_from_genome.py -g target_transcriptome.gtf \
  -f ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa -o ref2.fasta

# Combine annotations (if you have miRNA GTF)
# concatenate_gtf.py -m mature_mirna.gtf -t target_transcriptome.gtf -o combined.gtf

# Step 2: Collapse reads
echo "Step 2: Collapsing reads..."
chira_collapse.py -i raw_reads.fastq -o collapsed.fasta -u 12

# Step 3: Map reads (using 8 processes for parallelization)
echo "Step 3: Mapping reads..."
chira_map.py -i collapsed.fasta -o mapping -f1 ref1.fasta -f2 ref2.fasta \
  -b -a bwa -s both -p 8

# Step 4: Merge alignments (using 8 processes for parallel transcript processing)
echo "Step 4: Merging alignments..."
chira_merge.py -b mapping/mapped.bed -o merge -g target_transcriptome.gtf \
  -f1 ref1.fasta -f2 ref2.fasta -ao 0.7 -so 0.7 -p 8

# Step 5: Quantify CRLs (using 8 threads for parallelization)
echo "Step 5: Quantifying CRLs..."
chira_quantify.py -b merge/segments.bed -m merge/loci.txt -o quantify \
  -cs 0.7 -ls 10 -t 8

# Step 6: Extract interactions (using 8 processes for parallelization)
echo "Step 6: Extracting interactions..."
chira_extract.py -l quantify/loci.txt -o extract -f1 ref1.fasta -f2 ref2.fasta \
  -n $SAMPLE -g target_transcriptome.gtf -r -s -tc 0.1 -sc 0.5 -p 8

echo "Analysis complete! Results in extract/"
```

---

### Docker Workflow

If using Docker, the workflow is similar but commands are prefixed with `docker run`. Use the pre-built image:

```bash
# Pull the pre-built image (first time only)
docker pull docker.io/nemat1976/chiraplus:v0.0.1

# Run each step with volume mounts
docker run --rm -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  docker.io/nemat1976/chiraplus:v0.0.1 chira_collapse.py -i data/input.fastq -o output/collapsed.fasta

docker run --rm -v $(pwd)/data:/app/data -v $(pwd)/output:/app/output \
  docker.io/nemat1976/chiraplus:v0.0.1 chira_map.py -i output/collapsed.fasta -o output/mapping \
  -f1 data/ref1.fasta -f2 data/ref2.fasta -b -a bwa

# ... continue with remaining steps
```

**Note:** Alternatively, you can build the image from the provided Dockerfile:
```bash
docker build -t chira:latest .
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
  - **`fw`** (forward/transcript strand): Use for stranded libraries where reads align to the transcript strand (recommended for most protocols like CLASH, CLEAR-CLIP, PARIS, SPLASH)
  - **`rc`** (reverse complement): Use for stranded libraries where reads align to the reverse complement strand (rare, check your protocol)
  - **`both`**: Use only for unstranded libraries where strand information is not preserved (rare for interactome protocols)
  - **Why it matters**: Stranded mapping filters out alignments on the wrong strand, reducing false positives and improving chimeric read detection accuracy
- `-l1, --seed_length1`: Seed length for 1st mapping iteration (default: 12)
- `-l2, --seed_length2`: Seed length for 2nd mapping iteration (default: 16)
- `-s1, --align_score1`: Minimum alignment score for 1st iteration (default: 18)
- `-s2, --align_score2`: Minimum alignment score for 2nd iteration (default: 16)
- `-co, --chimeric_overlap`: Max bases between chimeric segments (default: 2)
- `-p, --processes`: Total number of CPU processes/threads to use (default: auto-detects CPU count)
  - **Automatic distribution**: The script automatically divides total processes among parallel BWA jobs
    - **Without chunking**: Total processes divided among parallel BWA jobs (2-4 jobs depending on indices)
    - **With chunking**: Total processes divided among parallel chunks, then each chunk uses its allocated processes
  - **Multi-threading**: Enables parallel processing for `samtools view`, `pysam.merge`, `pysam.sort`, and `sort` commands
  - **Performance**: 2-4x faster for large BAM files and sorting operations
  - **Memory optimization**: Use `--sort_memory` to specify memory per thread (e.g., "2G", "3G"), or install `psutil` for automatic optimization
  - **Recommendation**: Set to total number of CPU cores for optimal performance
- `--sort_memory`: Memory per thread for BAM sorting (e.g., "2G", "3G", optional)
  - If not specified, automatically calculates based on available RAM (requires `psutil`)
  - Falls back to safe default (2GB per thread) if `psutil` unavailable
  - **Note**: Total memory used = sort_memory × processes, so ensure sufficient RAM
- `--chunk_fasta`: Split input FASTA into N chunks for parallel processing (optional, recommended for very large files >1GB)
  - **How it works**: 
    - Creates N chunks from the input FASTA file
    - Processes chunks in batches, with parallel execution controlled by `--parallel_chunks` (default: 2)
    - Each chunk processes all BWA jobs sequentially
    - Remaining chunks are processed in subsequent batches automatically
  - **Example 1**: `--chunk_fasta 10 --processes 32 --parallel_chunks 2` (default)
    - Creates 10 chunks from FASTA
    - Runs 2 chunks in parallel (16 processes each, using all 32 processes)
    - Processes remaining 8 chunks in subsequent batches of 2
  - **Example 2**: `--chunk_fasta 10 --processes 32 --parallel_chunks 4`
    - Creates 10 chunks from FASTA
    - Runs 4 chunks in parallel (8 processes each, using all 32 processes)
    - Processes remaining 6 chunks in subsequent batches of 4
  - **Benefits**: Better I/O performance and memory efficiency for large datasets
  - Each chunk is processed independently through all BWA jobs, then results are merged
- `--parallel_chunks`: Number of chunks to process simultaneously when using `--chunk_fasta` (default: 2)
  - **How to set**: Based on available memory and CPU resources
    - **Default: 2** (recommended for most systems)
    - Small systems (<16 cores, <32GB RAM): `--parallel_chunks 1`
    - Medium systems (16-32 cores, 32-64GB RAM): `--parallel_chunks 2` (default)
    - Large systems (>32 cores, >64GB RAM): `--parallel_chunks 2-4`
    - Very large systems (>64 cores, >128GB RAM): `--parallel_chunks 4-8`
  - **Memory consideration**: Each chunk needs ~2-4GB RAM, so total memory ≈ `parallel_chunks × 4GB`
  - **CPU consideration**: Processes per chunk = `--processes / --parallel_chunks`
    - Each chunk should get at least 4 processes for optimal BWA performance
    - Example: `--processes 32 --parallel_chunks 4` → 8 processes per chunk (good)
    - Example: `--processes 8 --parallel_chunks 4` → 2 processes per chunk (too few, will auto-reduce)
  - **Note**: Only takes effect when `--chunk_fasta` is specified
- `--use_batchtools`: Enable batchtools for HPC cluster job submission (optional, for LSF/SLURM clusters)
  - **Requirements**: R with `batchtools` and `jsonlite` packages installed
  - **Benefits**: Submit chunk jobs to cluster scheduler for true parallel processing across cluster nodes
  - **Path handling**: All file paths are automatically converted to absolute paths for cluster job execution
  - **Additional parameters** (required when `--use_batchtools` is specified):
    - `--batchtools_queue`: LSF queue name (e.g., `long`, `short`)
    - `--batchtools_cores`: Cores per LSF job (e.g., `8`)
    - `--batchtools_memory`: Total memory per LSF job (e.g., `8GB`, automatically converted to per-core for LSF)
    - `--batchtools_walltime`: Walltime limit per job (e.g., `240:00`)
    - `--batchtools_max_parallel`: Max concurrent running jobs (optional, default: unlimited)
    - `--batchtools_conda_env`: Conda environment path (optional, auto-detected if not specified)
    - `--batchtools_template`: LSF template file path (optional, default: `lsf_custom.tmpl`, can use built-in `"lsf-simple"`)
  - **Example**: See [BATCHTOOLS_USAGE.md](BATCHTOOLS_USAGE.md) for detailed usage instructions and examples

**Outputs:**
- `sorted.bam`: Sorted BAM file
- `sorted.bed`: BED file containing all alignments
- `unmapped.fasta`: FASTA file with unmapped reads (optional)

**Usage Example:**
```bash
# Basic usage with 32 total processes (automatically divided among BWA jobs)
chira_map.py -i reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both -p 32

# Use all available CPU cores automatically (default behavior)
chira_map.py -i reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both

# With manual memory specification for BAM sorting
chira_map.py -i reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both -p 8 --sort_memory 3G

# For very large FASTA files (>1GB), use chunking for better I/O performance
# Creates 10 chunks, processes them in batches (default: 2 chunks at a time)
chira_map.py -i large_reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both -p 32 --chunk_fasta 10

# Example 1: Default behavior (2 chunks in parallel)
# - Creates 10 chunks from FASTA
# - Runs 2 chunks in parallel (16 processes each, using all 32 processes)
# - Processes remaining 8 chunks in subsequent batches of 2

# Example 2: Custom parallel chunks (4 chunks in parallel)
# - Creates 10 chunks from FASTA
# - Runs 4 chunks in parallel (8 processes each, using all 32 processes)
# - Processes remaining 6 chunks in subsequent batches of 4
chira_map.py -i large_reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both -p 32 --chunk_fasta 10 --parallel_chunks 4

# Example 3: Using batchtools for HPC cluster submission (LSF)
# - Splits FASTA into 20 chunks
# - Submits 20 independent LSF jobs (one per chunk)
# - Each job runs on different cluster node with 8 cores and 8GB memory
# - All paths automatically converted to absolute paths for cluster execution
chira_map.py -i large_reads.fasta -o output_dir \
   -f1 ref1.fasta -f2 ref2.fasta -a bwa -s both \
   --chunk_fasta 20 --use_batchtools \
   --batchtools_queue long \
   --batchtools_cores 8 \
   --batchtools_memory 8GB \
   --batchtools_walltime 240:00
```

**Understanding Chunking and Process Management:**

When using `--chunk_fasta`, the script uses a two-stage approach:

1. **FASTA Splitting**: The input FASTA is split into N chunks (as specified by `--chunk_fasta`)
   - Each chunk contains a subset of reads distributed round-robin
   - Chunks are typically 1-3GB each for optimal I/O performance

2. **Parallel Processing**: Chunks are processed in batches, with parallelism controlled by `--parallel_chunks` (default: 2)
   - Number of chunks running simultaneously is set by `--parallel_chunks` (default: 2)
   - Processes per chunk = `--processes / --parallel_chunks`
   - Each chunk should get at least 4 processes for optimal BWA performance
   - If processes per chunk < 4, the number of parallel chunks is automatically reduced
   - Example: 32 processes with `--parallel_chunks 2` → 2 chunks run in parallel (16 processes each)
   - Example: 32 processes with `--parallel_chunks 4` → 4 chunks run in parallel (8 processes each)
   - Example: 8 processes with `--parallel_chunks 2` → 1 chunk runs at a time (8 processes each, auto-reduced)

3. **Batch Processing**: If you have more chunks than `--parallel_chunks`, remaining chunks are processed in subsequent batches
   - Example: 10 chunks with `--parallel_chunks 2` → batches of 2: (2, 2, 2, 2, 2)
   - Example: 10 chunks with `--parallel_chunks 4` → batches of 4: (4, 4, 2)
   - Example: 10 chunks with `--parallel_chunks 1` → batches of 1: (1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

**Why This Approach:**
- **Prevents oversubscription**: Won't try to run more chunks than you have processes for
- **Efficient resource use**: All available processes are utilized without waste
- **Better I/O**: Chunking still provides I/O benefits even with limited parallelism
- **User control**: Adjust `--parallel_chunks` to match your system's memory and CPU resources
- **Automatic fallback**: If processes per chunk < 4, automatically reduces parallel chunks

**Recommendations:**
- For large files (>1GB): Use `--chunk_fasta 10-20` to create manageable chunks
- Set `--processes` to your total CPU core count for optimal performance
- Set `--parallel_chunks` based on your system:
  - **Small systems** (<16 cores, <32GB RAM): `--parallel_chunks 1`
  - **Medium systems** (16-32 cores, 32-64GB RAM): `--parallel_chunks 2` (default)
  - **Large systems** (>32 cores, >64GB RAM): `--parallel_chunks 2-4`
  - **Very large systems** (>64 cores, >128GB RAM): `--parallel_chunks 4-8`
- Ensure each chunk gets at least 4 processes: `--processes / --parallel_chunks >= 4`

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
- `-p, --processes`: Number of parallel processes for transcript processing (default: None, auto-detects CPU count)
  - **Chunk-based multiprocessing**: Groups transcripts into chunks (~1000 transcripts per chunk) and processes chunks in parallel
  - **Scalability**: Efficiently handles very large datasets (e.g., human genome with 387K+ transcripts)
  - **Performance**: 4-8x faster for datasets with many transcripts
  - **Note**: Changed from `-t, --threads` in v1.4.6 to reflect multiprocessing implementation
  - **Optimization**: For human genome (387K transcripts), creates ~388 manageable chunks instead of processing each transcript individually, reducing overhead significantly

**Outputs:**
- `segments.bed`: A BED file with reads categorized into segments
- `merged.bed`: A tabular file with merged alignments
  - Column 4: All alignments merged into that location
  - Column 5: Number of reads supporting the locus

**Usage Example:**
```bash
# Basic usage with 8 processes (auto-detects CPU count if not specified)
chira_merge.py -b mapped.bed -o output_dir -g annotation.gtf -f1 ref1.fasta -ao 0.7 -so 0.7 -p 8

# Auto-detect CPU count (default behavior)
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
- `-t, --threads`: Number of processes for parallel processing (default: 1, use 0 for all available cores)
  - **Multiprocessing**: Uses MPIRE WorkerPool to parallelize EM algorithm E-step (multimapped reads) and aggregation step (bypasses Python GIL)
  - **MPIRE benefits**: 50-90% memory reduction, 2-3x faster startup, better performance (required dependency)
  - **Performance**: 2-8x faster for large datasets with many multimapping reads
  - **Automatic**: Falls back to sequential processing for small datasets (<500 reads or num_processes × 50) to avoid process overhead

**Outputs:**
- `loci.counts`: Tabular file containing reads, their CRLs, and TPM values
  - Each line represents a read-CRL association with TPM quantification
  - Used as input for `chira_extract.py`

**Usage Example:**
```bash
# Basic usage with 8 threads
chira_quantify.py -b segments.bed -m loci.txt -o output_dir -cs 0.7 -ls 10 -t 8

# Use all available CPU cores automatically
chira_quantify.py -b segments.bed -m loci.txt -o output_dir -cs 0.7 -ls 10 -t 0
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
- **Automatically converts U (uracil) to T (thymine)** in sequences for ChiRA compatibility (ChiRA expects DNA sequences)

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

### extract_transcripts_from_genome.py

**Description:** Extracts transcript FASTA sequences from a genome FASTA file using gffread based on a filtered GTF file.

This script uses gffread (from [GFF utilities](https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)) to extract transcript sequences directly from the genome FASTA file based on transcript features in a filtered GTF file (e.g., output from `remove_mirna_hairpin_from_gtf.py`). This approach is more accurate than filtering pre-extracted transcript sequences, as it extracts sequences directly from the genome coordinates.

**Key Features:**
- Uses gffread to extract transcript sequences from genome FASTA
- Works with filtered GTF files (e.g., miRNA-removed GTF from `remove_mirna_hairpin_from_gtf.py`)
- Produces transcriptome FASTA without miRNA sequences (if GTF was filtered)
- Automatically validates that gffread is available

**Required Inputs:**
- `-g, --gtf`: Filtered GTF file (e.g., output from `remove_mirna_hairpin_from_gtf.py`)
- `-f, --genome-fasta`: Genome FASTA file (e.g., primary assembly from Ensembl)
- `-o, --output`: Output transcript FASTA file

**Dependencies:**
- Requires `gffread` (available via conda: `conda install -c bioconda gffread`)
- gffread is part of the GFF utilities package from Johns Hopkins University

**Usage Example:**
```bash
# Extract transcripts from genome using filtered GTF (without miRNAs)
extract_transcripts_from_genome.py -g target_transcriptome.gtf \
  -f Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -o target_transcriptome.fasta
```

**Note:** This script replaces the previous `remove_mirna_hairpin_from_fasta.py` approach. Instead of filtering pre-extracted transcript sequences, it extracts sequences directly from the genome FASTA based on the filtered GTF coordinates, ensuring accuracy and completeness.

---

### concatenate_gtf.py

**Description:** Concatenates mature miRNA GTF/GFF3 file with target transcriptome GTF file, removing comment lines from the miRNA GTF.

**Key Features:**
- Combines miRNA and target GTF files for split-reference analysis
- Removes comment lines from miRNA GTF
- Optionally removes comment lines from target GTF
- Note: miRBase provides GFF3 format, which ChiRA can handle directly. This script accepts GTF format.

**Required Inputs:**
- `-m, --mirna-gtf`: Mature miRNA GTF/GFF3 file (e.g., from miRBase via `download_mirbase_gff3.py`)
  - Note: miRBase GFF3 format can be used directly with ChiRA. This script accepts GTF format.
- `-t, --target-gtf`: Target transcriptome GTF file (output from `remove_mirna_hairpin_from_gtf.py`)
- `-o, --output`: Output combined GTF file

**Optional Parameters:**
- `--remove-target-comments`: Remove comment lines from target GTF as well

**Usage Example:**
```bash
concatenate_gtf.py -m mature_mirna.gtf -t target_no_mirna.gtf -o combined.gtf
```

---
