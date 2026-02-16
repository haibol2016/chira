# Using batchtools for LSF Cluster Parallel Processing

## Overview

ChiRA supports R batchtools as a backend for submitting chunk processing jobs to LSF clusters. This enables **true parallel computing** by distributing chunks across multiple cluster nodes, each running independently. This is simpler and more reliable than Dask for LSF environments.

**Key Benefits:**
- **Parallel Processing**: Each chunk runs as an independent LSF job on different cluster nodes
- **Scalability**: Process hundreds of chunks simultaneously across the cluster
- **Resource Control**: Specify cores, memory, and walltime per job
- **Visibility**: All jobs visible in `bjobs` for easy monitoring
- **Reliability**: Each job is independent - failures don't affect other jobs

## Prerequisites

1. **R installed** with `batchtools` package:
   ```r
   install.packages("batchtools")
   ```

2. **LSF scheduler** available on your cluster

3. **Shared filesystem** accessible to all compute nodes

4. **Conda environment** (optional but recommended) with ChiRA dependencies

## Quick Start: Running with Batchtools

### Basic Example

```bash
python chira_map.py \
  -i input.fasta \
  -o output_dir \
  --index1 index1.fa \
  --index2 index2.fa \
  --chunk_fasta 20 \
  --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 8 \
  --batchtools_memory 8GB
```

**What happens:**
1. FASTA file is split into 20 chunks
2. 20 independent LSF jobs are submitted (one per chunk)
3. Each job runs on a different cluster node with 8 cores and 8GB memory
4. Jobs run in parallel across the cluster
5. Results are automatically collected and merged

### All batchtools Options

- `--use_batchtools`: **REQUIRED** - Enable batchtools mode for LSF cluster parallel processing
- `--batchtools_queue`: LSF queue name (default: `long`)
- `--batchtools_cores`: **Cores per LSF job** (default: 8)
  - Each chunk job runs with this many cores on a cluster node
  - Example: `--batchtools_cores 16` → each job uses 16 cores
  - **Independent of `--processes`** (which applies to main job only)
- `--batchtools_memory`: **Total memory per LSF job** (automatically converted to per-core for LSF), e.g., `16GB` (default: auto-calculated)
  - **IMPORTANT**: You specify TOTAL memory, but LSF `rusage[mem=...]` is PER CORE
  - Code automatically converts: `total_memory ÷ cores = per_core_memory`
  - Example: `--batchtools_cores 8 --batchtools_memory 16GB` → LSF gets `rusage[mem=2GB]` (16GB ÷ 8 = 2GB per core)
  - Auto-calculation: `cores × 0.5GB` total (e.g., 8 cores → 4GB total → 0.5GB per core)
- `--batchtools_walltime`: Walltime limit per job, e.g., `240:00` (default: `240:00`)
- `--batchtools_max_parallel`: **Max concurrent running jobs** (default: unlimited, all submitted at once)
  - Limits how many jobs run simultaneously
  - Example: `--batchtools_max_parallel 2` → only 2 jobs run at a time, others wait
  - Useful if cluster has job limits per user
- `--batchtools_conda_env`: Conda environment path (default: auto-detect from `CONDA_DEFAULT_ENV`)
- `--batchtools_template`: LSF template file (default: `lsf_custom.tmpl` in ChiRA directory)
  - Uses proven template based on InPAS implementation
  - Only change if you need custom LSF directives

### Example: Full Command with Parallel Processing

```bash
python chira_map.py \
  -i reads.fasta \
  -o results \
  --index1 index1.fa \
  --index2 index2.fa \
  --chunk_fasta 20 \
  --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 8 \
  --batchtools_memory 8GB \
  --batchtools_walltime 240:00 \
  --batchtools_max_parallel 10 \
  --batchtools_conda_env ~/miniconda3/envs/chira
```

**This command:**
- Splits `reads.fasta` into 20 chunks
- Submits 20 LSF jobs to the `long` queue
- Each job uses 8 cores and 8GB total memory (1GB per core)
- Limits to 10 concurrent jobs (if cluster has limits)
- All 20 chunks process in parallel across cluster nodes

### Example: Limiting Concurrent Jobs

If your cluster limits concurrent jobs per user, use `--batchtools_max_parallel`:

```bash
python chira_map.py \
  -i large.fasta \
  -o results \
  --index1 index1.fa \
  --index2 index2.fa \
  --chunk_fasta 50 \
  --use_batchtools \
  --batchtools_queue long \
  --batchtools_cores 16 \
  --batchtools_memory 16GB \
  --batchtools_max_parallel 2
```

**This ensures:**
- Only 2 jobs run at a time (respects cluster limits)
- All 50 chunks will be processed sequentially in batches of 2
- Each batch waits for completion before submitting the next

## How Parallel Processing Works

1. **Chunk Preparation**: Python splits FASTA file into N chunks (specified by `--chunk_fasta`)
2. **Job Submission**: Python calls R script (`submit_chunks_batchtools.R`) that uses batchtools to submit N independent LSF jobs
3. **Parallel Execution**: Each LSF job runs on a **different cluster node**:
   - Job 1 processes chunk 1 on node A
   - Job 2 processes chunk 2 on node B
   - Job 3 processes chunk 3 on node C
   - ... and so on
4. **Resource Allocation**: Each job gets:
   - `--batchtools_cores` cores
   - `--batchtools_memory` total memory
   - Runs independently with its own resources
5. **Result Collection**: Python monitors completion and collects BAM files from each chunk
6. **Merging**: Python merges all chunk BAMs into final BAM files

**Key Point**: This is **true parallel processing** - chunks run simultaneously on different nodes, not sequentially on one node.

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

**Default Template**: ChiRA uses `lsf_custom.tmpl` by default, which is based on the proven InPAS implementation. This template is included with ChiRA and works for most LSF clusters.

**The default template (`lsf_custom.tmpl`) handles:**
- Queue specification (`#BSUB -q`)
- Resource requirements (`#BSUB -n`, `#BSUB -R rusage[mem=...]`)
- Walltime (`#BSUB -W`)
- Job name (`#BSUB -J`)
- Log files (`#BSUB -o`, `#BSUB -e`)
- Host spanning (`#BSUB -R "span[hosts=1]"`)
- Threading environment variables (OMP_NUM_THREADS, etc.)

**To use a different template:**
```bash
--batchtools_template /path/to/your/custom.tmpl
```

**To use built-in "lsf-simple" (not recommended):**
```bash
--batchtools_template lsf-simple
```

**Only create a custom template if:**
- You need to load specific modules (e.g., `module load python/3.9`)
- You need additional custom LSF directives
- Your cluster has non-standard LSF configuration

See `lsf_custom.tmpl` for the default template structure.

## Files Created

- `process_chunk_batchtools.py`: Standalone script for processing one chunk
- `submit_chunks_batchtools.R`: R script that submits jobs via batchtools
- `lsf_custom.tmpl`: Example custom template (only needed if built-in doesn't work)
- `output_dir/batchtools_registry_<timestamp>/`: Batchtools registry directory

