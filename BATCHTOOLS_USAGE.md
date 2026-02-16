# Using batchtools for LSF Cluster Processing

## Overview

ChiRA now supports R batchtools for submitting chunk processing jobs to LSF clusters. This is simpler and more reliable than Dask for LSF environments.

## Prerequisites

1. **R installed** with `batchtools` package:
   ```r
   install.packages("batchtools")
   ```

2. **LSF scheduler** available on your cluster

3. **Shared filesystem** accessible to all compute nodes

4. **Conda environment** (optional but recommended) with ChiRA dependencies

## Usage

### Basic Usage

```bash
python chira_map.py \
  -i input.fasta \
  -o output_dir \
  --chunk_fasta 20 \
  --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 8 \
  --batchtools_memory 8GB
```

### All batchtools Options

- `--use_batchtools`: Enable batchtools mode (instead of Dask or single-node)
- `--batchtools_queue`: LSF queue name (default: `long`)
- `--batchtools_cores`: Cores per job (default: auto-calculated)
- `--batchtools_memory`: Memory per job, e.g., `8GB` (default: auto-calculated)
- `--batchtools_walltime`: Walltime limit, e.g., `240:00` (default: `240:00`)
- `--batchtools_max_parallel`: Max concurrent chunks (default: unlimited, all submitted at once)
- `--batchtools_conda_env`: Conda environment name (default: auto-detect from `CONDA_DEFAULT_ENV`)

### Example: Full Command

```bash
python chira_map.py \
  -i reads.fasta \
  -o results \
  -x1 index1.fa \
  -x2 index2.fa \
  --chunk_fasta 20 \
  --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 8 \
  --batchtools_memory 8GB \
  --batchtools_walltime 240:00 \
  --batchtools_conda_env chira_env
```

## How It Works

1. **Chunk Preparation**: Python splits FASTA into chunks
2. **Job Submission**: Python calls R script that uses batchtools to submit LSF jobs
3. **Job Execution**: Each LSF job runs `process_chunk_batchtools.py` to process one chunk
4. **Result Collection**: Python polls for completion markers and collects BAM files
5. **Merging**: Python merges chunk BAMs into final BAM files

## Monitoring Jobs

### View Submitted Jobs

```bash
# View all your jobs
bjobs -u $USER

# View jobs with specific prefix
bjobs -J chira_bt_*
```

### Check Batchtools Registry

The registry directory is created in your output directory:
```
output_dir/batchtools_registry_<timestamp>/
```

You can check job status in R:
```r
library(batchtools)
reg <- loadRegistry("output_dir/batchtools_registry_<timestamp>")
getStatus(reg)
```

## Advantages Over Dask

1. **Simpler**: Direct LSF job submission, no worker management
2. **More Reliable**: Each chunk is an independent LSF job
3. **Better Visibility**: All jobs visible in `bjobs`
4. **Easier Debugging**: Each job has its own log file
5. **No Network Overhead**: No task serialization/deserialization

## Troubleshooting

### R/batchtools Not Found

Ensure R is in your PATH:
```bash
which Rscript
```

Install batchtools if needed:
```r
install.packages("batchtools")
```

### Jobs Not Starting

Check LSF queue status:
```bash
bqueues
```

Check job requirements match your LSF configuration:
- `--batchtools_cores` should match `#BSUB -n` in your job script
- `--batchtools_memory` should match `#BSUB -R rusage[mem=...]`

### Conda Environment Not Found

Specify conda environment explicitly:
```bash
--batchtools_conda_env chira_env
```

Or ensure `CONDA_DEFAULT_ENV` is set:
```bash
export CONDA_DEFAULT_ENV=chira_env
```

### nodename Parameter

The `nodename = "localhost"` parameter in `submit_chunks_batchtools.R` (line 51) specifies the hostname where the LSF scheduler runs.

**For most clusters: Keep as `"localhost"`**
- You submit jobs from the login node where LSF commands are available
- LSF scheduler runs on the same node
- This is the default and works for 99% of LSF clusters

**Only change if:**
- Your cluster documentation specifies a different submission host
- LSF scheduler runs on a different machine
- You get connection errors when submitting jobs

**To auto-detect your hostname:**
```r
nodename = Sys.info()["nodename"]
```

**To manually specify:**
```r
nodename = "your-lsf-hostname"
```

## Template Files

**You typically DON'T need a custom template file.** The code uses batchtools' built-in `"lsf-simple"` template, which works for most LSF clusters.

**Only create a custom template if:**
- The built-in template doesn't work with your LSF setup
- You need to load specific modules (e.g., `module load python/3.9`)
- You need custom LSF directives
- Your cluster has non-standard LSF configuration

**To use a custom template:**
1. Create a template file (see `lsf_custom.tmpl` for an example)
2. Modify `submit_chunks_batchtools.R` line 47:
   ```r
   template = "/path/to/your/custom_lsf.tmpl"  # Instead of "lsf-simple"
   ```

The built-in template handles:
- Queue specification (`#BSUB -q`)
- Resource requirements (`#BSUB -n`, `#BSUB -R rusage[mem=...]`)
- Walltime (`#BSUB -W`)
- Job name (`#BSUB -J`)
- Log files (`#BSUB -o`, `#BSUB -e`)

## Files Created

- `process_chunk_batchtools.py`: Standalone script for processing one chunk
- `submit_chunks_batchtools.R`: R script that submits jobs via batchtools
- `lsf_custom.tmpl`: Example custom template (only needed if built-in doesn't work)
- `output_dir/batchtools_registry_<timestamp>/`: Batchtools registry directory

