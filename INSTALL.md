# ChiRA Installation Guide

## Installing ChiRA in a Conda Environment

### Step 1: Create and Activate Conda Environment

```bash
# Create a new conda environment with Python 3.9 (or 3.10, 3.11)
conda create -n chira python=3.9
conda activate chira
```

### Step 2: Install Required External Tools

ChiRA requires several command-line tools that must be installed separately:

```bash
# Install required tools from bioconda
conda install -c bioconda bwa samtools bedtools

# Optional: Install additional tools if needed
conda install -c bioconda intarna gffread  # Only if using specific features
```

### Step 3: Install ChiRA Package

#### Option A: Install from Source (Recommended for Development)

```bash
# Navigate to ChiRA directory
cd /path/to/chira

# Install in editable mode (allows code modifications)
pip install -e .

# Or install normally
pip install .
```

#### Option B: Install with Optional Dependencies

```bash
# Install with all optional dependencies (recommended for best performance)
pip install -e .[optional]

# This installs:
# - psutil (for memory optimization and I/O performance)
# - requests (for downloading Ensembl files)
# - pyliftover (for coordinate liftover)
```

### Step 4: Verify Installation

```bash
# Check that ChiRA scripts are available
chira_map.py --version
chira_merge.py --version
chira_extract.py --version

# Check that external tools are available
bwa
samtools --version
bedtools --version
```

## Complete Installation Example

```bash
# 1. Create environment
conda create -n chira python=3.9
conda activate chira

# 2. Install external tools
conda install -c bioconda bwa samtools bedtools

# 3. Navigate to ChiRA directory
cd /path/to/chira

# 4. Install ChiRA with optional dependencies
pip install -e .[optional]

# 5. Verify installation
chira_map.py --version
```

## Installing for Batchtools Support

If you plan to use batchtools for LSF cluster processing:

```bash
# Activate your conda environment
conda activate chira

# Install R and batchtools
conda install -c conda-forge r-base
# Then in R:
# install.packages("batchtools")
# Or via conda:
conda install -c conda-forge r-batchtools
```

## Troubleshooting

### Python Dependencies

If you encounter import errors:

```bash
# Reinstall dependencies
pip install --upgrade biopython bcbio-gff pysam

# Install optional dependencies for better performance
pip install psutil
```

### External Tools Not Found

Ensure tools are in your PATH:

```bash
# Check if tools are accessible
which bwa
which samtools
which bedtools

# If not found, ensure conda environment is activated
conda activate chira
```

### Permission Errors

If you get permission errors during installation:

```bash
# Install to user directory
pip install --user -e .

# Or use conda environment (recommended)
conda activate chira
pip install -e .
```

## Uninstalling

```bash
# Deactivate environment
conda deactivate

# Remove environment
conda env remove -n chira

# Or uninstall package only
pip uninstall chira
```

## Dependencies Summary

### Required Python Packages (auto-installed)
- `biopython` - FASTA/sequence parsing
- `bcbio-gff` - GFF/GTF annotation parsing
- `pysam` - BAM file manipulation

### Optional Python Packages (recommended)
- `psutil` - Memory optimization and I/O performance
- `requests` - Downloading Ensembl files
- `pyliftover` - Coordinate liftover

### Required External Tools (install via conda)
- `bwa` - Sequence alignment
- `samtools` - BAM file processing
- `bedtools` - Genomic interval operations

### Optional External Tools
- `intarna` - RNA-RNA interaction prediction
- `gffread` - GFF utilities

## Notes

- **Python Version**: Requires Python >= 3.6 (Python 3.9+ recommended)
- **Conda Channel**: Use `bioconda` channel for bioinformatics tools
- **Editable Install**: Use `pip install -e .` during development to see code changes immediately
- **Batchtools**: Required only if using `--use_batchtools` option in `chira_map.py`

