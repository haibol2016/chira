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

4. **pyliftover**
   - Used in: `download_mirbase_gff3.py` (optional, only when using coordinate liftover)
   - Import: `from pyliftover import LiftOver`
   - Purpose: Converting coordinates between different genome versions
   - Note: Only required if using `--source-genome`, `--target-genome`, and `--chain-file` options in `download_mirbase_gff3.py`

5. **requests**
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

For CLAN and blockbuster, please refer to their respective documentation for installation instructions.

## Script-Specific Dependencies

### Core ChiRA Scripts

- **chira_collapse.py**: No external dependencies (uses standard library only)
- **chira_map.py**: Requires `pysam`, `bwa`, `samtools` (CLAN optional)
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

- Standard library modules used (no installation needed):
  - `argparse`, `os`, `sys`, `collections`, `multiprocessing`, `itertools`, `datetime`, `subprocess`, `math`, `re`, `copy`, `gzip`, `shutil`, `ftplib`, `urllib`, `tempfile`, `time`

- All utility scripts (`download_ensembl.py`, `download_mirbase_mature.py`, `gff3_to_gtf.py`, etc.) are standalone and can be used independently of the main ChiRA pipeline.
