# Testing Batchtools Submission Interactively

This guide shows how to test batchtools LSF job submission interactively to diagnose why jobs aren't visible in `bjobs`.

**Note**: All file paths in configuration files are automatically converted to absolute paths. This ensures cluster jobs can correctly resolve files regardless of working directory.

## Quick Test: Standalone R Script

### Step 1: Run the Interactive Test Script

```bash
# Make sure you're in the ChiRA directory
cd /path/to/chira

# Run the interactive test script
Rscript test_batchtools_interactive.R
```

This script will:
- ✓ Check if batchtools is installed
- ✓ Test LSF commands (`bqueues`, `bsub`)
- ✓ Verify template file exists
- ✓ Create a test registry
- ✓ Configure LSF cluster functions
- ✓ Submit a single test job
- ✓ Check if job appears in `bjobs`

### Step 2: Check the Output

Look for these key indicators:

**Success:**
```
✓ LSF batch ID found: 12345
Verify with: bjobs 12345
✓ Job visible in LSF:
```

**Failure:**
```
✗ No LSF batch ID found! Job was not submitted to LSF.
```

## Manual Testing: Test with Minimal Configuration

### Step 1: Create Test Config Files

Create a minimal test configuration:

```bash
# Create test directory
mkdir -p /tmp/batchtools_test
cd /tmp/batchtools_test

# Create minimal config.json
# Note: All paths should be absolute paths (handled automatically by chira_map.py)
cat > config.json << 'EOF'
{
  "reg_dir": "/tmp/batchtools_test/registry",
  "queue": "normal",
  "cores_per_job": 1,
  "memory_per_job": "1GB",
  "walltime": "00:05",
  "conda_env": "",
  "python_script": "/path/to/chira/process_chunk_batchtools.py",
  "chunk_dir": "/tmp/batchtools_test/chunks",
  "alignment_job_types_json": "[]",
  "per_chunk_processes": 1,
  "job_name_prefix": "test_job",
  "max_parallel": 1,
  "template_file": "/path/to/chira/lsf_custom.tmpl"
}
EOF

# Create minimal chunks.json
cat > chunks.json << 'EOF'
[
  {
    "chunk_file": "/tmp/test.fasta",
    "chunk_idx": 0
  }
]
EOF
```

### Step 2: Test R Script Directly

```bash
# Run the R script with test configs
Rscript /path/to/chira/submit_chunks_batchtools.R config.json chunks.json
```

### Step 3: Check Output

Look for:
- LSF batch IDs in the output
- Any error messages
- Registry directory contents

### Step 4: Verify in LSF

```bash
# Check if jobs are visible
bjobs -u $USER

# Or check specific job IDs if printed
bjobs <job_id>
```

## Testing from Python

### Step 1: Create a Minimal Test Script

Create `test_batchtools_python.py`:

```python
#!/usr/bin/env python
import sys
import os
sys.path.insert(0, '/path/to/chira')

from chira_map import submit_chunks_with_batchtools
import argparse

# Create minimal args object
class Args:
    def __init__(self):
        self.batchtools_queue = "normal"  # Change to your queue
        self.batchtools_cores = 1
        self.batchtools_memory = "1GB"
        self.batchtools_walltime = "00:05"
        self.batchtools_conda_env = ""
        self.batchtools_max_parallel = 1
        self.batchtools_template = None
        self.outdir = "/tmp/batchtools_test"

args = Args()

# Create a dummy chunk file
chunk_dir = "/tmp/batchtools_test/chunks"
os.makedirs(chunk_dir, exist_ok=True)
chunk_file = os.path.join(chunk_dir, "chunk_000.fasta")
with open(chunk_file, 'w') as f:
    f.write(">test\nATCG\n")

chunk_files = [chunk_file]
alignment_job_types = []  # Empty for testing
per_chunk_processes = 1

# Submit
reg_dir, job_ids = submit_chunks_with_batchtools(
    args, chunk_files, chunk_dir, alignment_job_types, per_chunk_processes
)

if reg_dir and job_ids:
    print(f"✓ Submitted {len(job_ids)} jobs")
    print(f"Registry: {reg_dir}")
    print(f"Job IDs: {job_ids}")
    
    # Check LSF
    import subprocess
    result = subprocess.run(['bjobs', '-u', os.environ.get('USER', '')], 
                          capture_output=True, text=True)
    print("\nLSF jobs:")
    print(result.stdout)
else:
    print("✗ Submission failed")
```

### Step 2: Run the Test

```bash
python test_batchtools_python.py
```

## Debugging Checklist

If jobs aren't visible in `bjobs`, check:

### 1. LSF Commands Available?
```bash
which bsub
which bqueues
bjobs -u $USER
```

### 2. Template File Valid?
```bash
# Check template exists
ls -l /path/to/chira/lsf_custom.tmpl

# Check template syntax (should have #BSUB directives)
head -20 /path/to/chira/lsf_custom.tmpl
```

### 3. Queue Name Correct?
```bash
# List available queues
bqueues

# Check queue details
bqueues -l normal  # Replace 'normal' with your queue
```

### 4. Batchtools Registry Logs
```bash
# Check registry directory
ls -la /tmp/batchtools_test/registry/

# Check logs
ls -la /tmp/batchtools_test/registry/logs/

# View log files
cat /tmp/batchtools_test/registry/logs/*.log
```

### 5. Test Direct bsub Submission
```bash
# Try submitting a job directly with bsub
bsub -q normal -n 1 -R "rusage[mem=1GB]" -J test_job "echo 'test'"
bjobs -J test_job
```

If direct `bsub` works but batchtools doesn't, the issue is in batchtools configuration.

## Common Issues and Solutions

### Issue: "No LSF batch ID found"

**Possible causes:**
1. Template file not found or invalid
2. Queue name incorrect
3. LSF commands not in PATH
4. Permission issues

**Solution:**
- Check template file path
- Verify queue name with `bqueues`
- Ensure `bsub` is in PATH
- Check batchtools logs

### Issue: "Cluster functions are NULL"

**Possible causes:**
1. Template file syntax error
2. batchtools version incompatibility

**Solution:**
- Check template file syntax
- Update batchtools: `install.packages("batchtools")`

### Issue: Jobs submitted but not visible

**Possible causes:**
1. Jobs immediately failed
2. Jobs in different queue
3. Jobs expired/finished

**Solution:**
- Check `bjobs -a` (all jobs including finished)
- Check `bhist` for job history
- Check batchtools registry status

## Advanced: Test in R Console

You can also test interactively in R:

```r
library(batchtools)

# Create registry
reg <- makeRegistry(file.dir = "/tmp/test_reg", conf.file = NA)

# Configure LSF
reg$cluster.functions <- makeClusterFunctionsLSF(
  template = "/path/to/chira/lsf_custom.tmpl",
  scheduler.latency = 1,
  fs.latency = 65
)

# Create test job
test_fun <- function(x) x * 2
ids <- batchMap(fun = test_fun, x = 1, reg = reg)

# Submit
resources <- list(
  ncpus = 1,
  mpp = "1GB",
  queue = "normal",
  walltime = 300
)
submitJobs(ids, resources = resources, reg = reg)

# Check status
Sys.sleep(2)
getJobTable(ids = ids, reg = reg)

# Check LSF
system("bjobs -u $USER")
```

## Next Steps

Once you identify the issue:
1. Fix the configuration (queue name, template path, etc.)
2. Re-run the test
3. Verify jobs appear in `bjobs`
4. Test with actual ChiRA workflow

