# ChiRA External Dependencies

This document lists all external dependencies required by the ChiRA package.

## Python Packages

The following Python packages need to be installed via pip or conda:

### Core Dependencies

1. **Biopython** (`Bio`)
   - Used in: `chira_extract.py`, `chira_utilities.py`, `remove_mirna_hairpin_from_fasta.py`
   - Import: `from Bio import SeqIO`
   - Purpose: Reading FASTA files (used in `extract_reflengths()` and for parsing FASTA sequences in `chira_extract.py` and `remove_mirna_hairpin_from_fasta.py`)
   - Note: `chira_collapse.py` no longer uses Biopython - uses raw file parsing for better performance

2. **bcbiogff** (`BCBio`)
   - Used in: `chira_extract.py`, `chira_merge.py`
   - Import: `from BCBio import GFF`
   - Purpose: Parsing GFF/GTF annotation files

3. **pysam**
   - Used in: `chira_map.py`
   - Import: `import pysam`
   - Purpose: Reading and manipulating BAM files, merging BAM files, sorting BAM files

### Optional Dependencies

4. **psutil**
   - Used in: `chira_map.py` (optional, for automatic memory detection)
   - Import: `import psutil`
   - Purpose: Automatically calculating optimal memory allocation for BAM sorting based on available system RAM
   - Note: If not available, falls back to safe default (2GB per thread). Install with `pip install psutil` or `conda install psutil` for automatic memory optimization.
   - Benefit: Prevents memory exhaustion and optimizes performance by using available RAM efficiently

5. **pyliftover**
   - Used in: `download_mirbase_gff3.py` (optional, only when using coordinate liftover)
   - Import: `from pyliftover import LiftOver`
   - Purpose: Converting coordinates between different genome versions
   - Note: Only required if using `--source-genome`, `--target-genome`, and `--chain-file` options in `download_mirbase_gff3.py`

6. **requests**
   - Used in: `download_ensembl.py`
   - Import: `import requests`
   - Purpose: Downloading files from Ensembl via HTTP/HTTPS

## External Command-Line Tools

The following command-line tools must be installed and available in the system PATH:

### Alignment Tools

1. **BWA (Burrows-Wheeler Aligner)**
   - Used in: `chira_map.py`
   - Commands: `bwa mem`, `bwa index`
   - Purpose: Mapping reads to reference transcriptome

2. **CLAN (Chimeric Long-read Aligner)**
   - Used in: `chira_map.py` (optional)
   - Commands: `clan_search`, `clan_output`, `clan_index`
   - Purpose: Alternative alignment tool for chimeric reads
   - Note: Only needed if using CLAN aligner instead of BWA (not recommended due to performance)

### Bioinformatics Utilities

3. **samtools**
   - Used in: `chira_map.py`
   - Commands: `samtools view`
   - Purpose: Converting SAM to BAM format

4. **BEDTools**
   - Used in: `chira_extract.py`, `chira_merge.py`
   - Commands: `bedtools intersect` (or `intersectBed` in older versions), `bedtools getfasta` (or `fastaFromBed` in older versions)
   - Purpose: Extracting sequences from BED files and calculating overlaps
   - Note: Code automatically detects and supports both old and new BEDTools command formats

5. **blockbuster.x**
   - Used in: `chira_merge.py` (optional, for block-based merging)
   - Purpose: Clustering and merging alignments

### RNA Interaction Tools

6. **IntaRNA**
   - Used in: `chira_extract.py` (optional, for hybridization)
   - Purpose: Predicting RNA-RNA interactions and hybridization

### System Utilities

7. **sort** (standard Unix/Linux tool)
   - Used in: `chira_extract.py`, `chira_map.py`, `chira_merge.py`, `chira_quantify.py`
   - Purpose: Sorting files (used via `os.system()`)
   - Note: GNU sort (part of GNU coreutils) version >= 8.6 supports `--parallel` option for faster sorting of large files. The code automatically detects and uses parallel sort when available, falling back to standard sort otherwise.
   - Requirement: GNU coreutils (standard on Linux, available via Homebrew on macOS: `brew install coreutils`)

8. **cat** (standard Unix/Linux tool)
   - Used in: `chira_extract.py`
   - Purpose: Concatenating files (used via `os.system()`)

9. **mv** (standard Unix/Linux tool)
   - Used in: `chira_merge.py`
   - Purpose: Moving/renaming files (used via `os.system()`)

## Installation Recommendations

### Python Packages

**Core packages:**
```bash
pip install biopython bcbiogff pysam requests
# or
conda install -c bioconda biopython bcbiogff pysam requests
```

**Optional packages:**
```bash
# For automatic memory optimization in chira_map.py
pip install psutil
# or
conda install -c bioconda psutil

# For coordinate liftover in gff3_to_gtf.py
pip install pyliftover
# or
conda install -c bioconda pyliftover
```

### Command-Line Tools

Most tools can be installed via conda:
```bash
conda install -c bioconda bwa samtools bedtools intarna
```

**GNU coreutils** (for parallel sort support):
- **Linux**: Usually pre-installed (GNU sort is standard)
- **macOS**: Install via Homebrew: `brew install coreutils` (provides `gsort` command, or use `gcoreutils` to get GNU versions of standard commands)
- **Note**: GNU sort version >= 8.6 is required for `--parallel` option. The code automatically detects GNU sort and uses parallel sorting when available.

For CLAN and blockbuster, please refer to their respective documentation for installation instructions.

## Parallel Computing Support

ChiRA scripts support parallel processing to improve performance on multi-core systems:

### Multi-Threading Support

1. **chira_quantify.py**
   - Uses Python's `ThreadPoolExecutor` for parallelizing the EM algorithm
   - Command-line option: `-t, --threads` (default: 1, use 0 for all available cores)
   - Parallelizes: E-step (multimapped reads) and aggregation step
   - Benefit: 2-4x faster for large datasets with many multimapping reads

2. **chira_merge.py**
   - Uses Python's `ThreadPoolExecutor` for parallelizing chromosome processing
   - Command-line option: `-t, --threads` (default: 1, use 0 for all available cores)
   - Parallelizes: Chromosome processing in overlap-based merging, parallel sort for blockbuster-based merging
   - Benefit: 2-4x faster for datasets with many chromosomes

3. **chira_map.py**
   - Uses multi-threading for external tools (`samtools`, `pysam`, `sort`)
   - Command-line option: `-p, --processes` (default: 1)
   - Parallelizes: `samtools view`, `pysam.merge`, `pysam.sort`, and `sort` commands
   - Automatic memory optimization: Uses `psutil` (optional) to calculate optimal memory per thread for BAM sorting
   - Benefit: 2-4x faster for large BAM files and sorting operations

### Multi-Processing Support

4. **chira_extract.py**
   - Uses Python's `multiprocessing.Process` for parallelizing read processing
   - Command-line option: `-p, --processes` (default: 1)
   - Parallelizes: Chimeras extraction and hybridization steps
   - Also uses parallel sort (GNU sort `--parallel`) for merging and interaction summary
   - Benefit: Linear speedup with number of processes (up to available cores)

### Performance Recommendations

- **For large datasets**: Use `-t 0` or `-p 0` to automatically use all available CPU cores
- **For memory-constrained systems**: Specify thread/process count explicitly to control memory usage
- **For BAM sorting**: Install `psutil` for automatic memory optimization, or use `--sort_memory` to specify memory per thread manually

## Script-Specific Dependencies

### Core ChiRA Scripts

- **chira_collapse.py**: No external dependencies (uses standard library only)
- **chira_map.py**: Requires `pysam`, `bwa`, `samtools` (CLAN optional, `psutil` optional for memory optimization)
- **chira_merge.py**: Requires `bcbiogff`, `BEDTools` (blockbuster optional)
- **chira_extract.py**: Requires `biopython`, `bcbiogff`, `BEDTools` (IntaRNA optional)
- **chira_quantify.py**: No external dependencies (uses standard library only)
- **chira_utilities.py**: Requires `biopython`

### Utility Scripts

- **download_ensembl.py**: Requires `requests`
- **download_mirbase_gff3.py**: Requires `pyliftover` (optional, only for coordinate liftover)
- **download_mirbase_mature.py**: No external dependencies (uses standard library only)
- **remove_mirna_hairpin_from_gtf.py**: No external dependencies (uses standard library only)
- **remove_mirna_hairpin_from_fasta.py**: Requires `biopython`
- **concatenate_gtf.py**: No external dependencies (uses standard library only)

## Notes

- Some tools are optional depending on usage:
  - **CLAN**: Only needed if using CLAN aligner instead of BWA (not recommended)
  - **blockbuster.x**: Only needed if using block-based merging method
  - **IntaRNA**: Only needed if using the `--hybridize` option in `chira_extract.py`
  - **pyliftover**: Only needed if using coordinate liftover in `download_mirbase_gff3.py`
  - **psutil**: Optional but recommended for automatic memory optimization in `chira_map.py`

- Standard library modules used (no installation needed):
  - `argparse`, `os`, `sys`, `collections`, `multiprocessing`, `itertools`, `datetime`, `subprocess`, `math`, `re`, `copy`, `gzip`, `shutil`, `ftplib`, `urllib`, `tempfile`, `time`, `concurrent.futures` (for ThreadPoolExecutor)

- Parallel computing features:
  - All parallel computing features are backward compatible (default to single-threaded/process)
  - Threading is used for CPU-bound tasks that can be safely parallelized (EM algorithm, chromosome processing)
  - Multiprocessing is used for I/O-bound tasks and external tool parallelization
  - GNU sort (GNU coreutils >= 8.6) parallel support is automatically detected and used when available
  - On macOS, GNU coreutils can be installed via Homebrew (`brew install coreutils`) to enable parallel sort

- All utility scripts (`download_ensembl.py`, `download_mirbase_mature.py`, `gff3_to_gtf.py`, etc.) are standalone and can be used independently of the main ChiRA pipeline.
