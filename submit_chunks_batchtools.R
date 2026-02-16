#!/usr/bin/env Rscript
# R script to submit chunk processing jobs using batchtools for LSF cluster
# Usage: Rscript submit_chunks_batchtools.R <config_json> <chunks_json>
#
# This script uses batchtools to submit LSF jobs for each chunk.
# Each job runs process_chunk_batchtools.py to process one chunk.

suppressPackageStartupMessages({
  library(batchtools)
  library(jsonlite)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript submit_chunks_batchtools.R <config_json> <chunks_json>")
}

config_file <- args[1]
chunks_file <- args[2]

# Read configuration
config <- fromJSON(config_file)
chunks <- fromJSON(chunks_file)

# Extract configuration
reg_dir <- config$reg_dir
queue <- config$queue
cores_per_job <- config$cores_per_job
memory_per_job <- config$memory_per_job
walltime <- config$walltime
conda_env <- config$conda_env
python_script <- config$python_script
chunk_dir <- config$chunk_dir
alignment_job_types_json <- config$alignment_job_types_json
per_chunk_processes <- config$per_chunk_processes
job_name_prefix <- config$job_name_prefix

# Create batchtools registry
if (dir.exists(reg_dir)) {
  unlink(reg_dir, recursive = TRUE)
}
reg <- makeRegistry(file.dir = reg_dir, work.dir = getwd())

# Configure LSF cluster functions
# Uses built-in "lsf-simple" template by default
# If you need a custom template, create a file and specify its path here
# Example: template = "/path/to/custom_lsf.tmpl"
reg$cluster.functions <- makeClusterFunctionsLSF(
  template = "lsf-simple",  # Built-in template - works for most LSF setups
  nodename = "localhost"    # Hostname where LSF scheduler runs
  # For most clusters: keep as "localhost" (you submit from login node)
  # Only change if: LSF runs on a different host, or cluster docs specify otherwise
  # To auto-detect: nodename = Sys.info()["nodename"]
)

# Define job function
process_chunk_job <- function(chunk_file, chunk_idx, chunk_dir, 
                              alignment_job_types_json, per_chunk_processes, 
                              python_script, conda_env) {
  # Build command
  cmd <- paste(
    "python3", python_script,
    shQuote(chunk_file),
    chunk_idx,
    shQuote(chunk_dir),
    shQuote(alignment_job_types_json),
    per_chunk_processes
  )
  
  # Add conda activation if specified
  if (!is.null(conda_env) && conda_env != "") {
    conda_init <- 'eval "$(conda shell.bash hook)" 2>/dev/null || source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true'
    conda_activate <- paste('conda activate', conda_env)
    cmd <- paste(conda_init, "&&", conda_activate, "&&", cmd)
  }
  
  # Execute command
  system(cmd)
}

# Submit jobs
cat(sprintf("Submitting %d chunk jobs to LSF cluster...\n", nrow(chunks)))

# Prepare job resources
resources <- list(
  queue = queue,
  cores = cores_per_job,
  memory = memory_per_job,
  walltime = walltime,
  job.name = paste0(job_name_prefix, "_chunk")
)

# Submit jobs
ids <- batchMap(
  fun = process_chunk_job,
  chunk_file = chunks$chunk_file,
  chunk_idx = chunks$chunk_idx,
  MoreArgs = list(
    chunk_dir = chunk_dir,
    alignment_job_types_json = alignment_job_types_json,
    per_chunk_processes = per_chunk_processes,
    python_script = python_script,
    conda_env = conda_env
  ),
  reg = reg
)

# Submit jobs with resources
submitJobs(ids, resources = resources, reg = reg)

cat(sprintf("Submitted %d jobs. Job IDs: %s\n", length(ids), paste(ids, collapse = ", ")))
cat(sprintf("Registry directory: %s\n", reg_dir))
cat(sprintf("Monitor jobs with: bjobs -J %s*\n", job_name_prefix))
cat(sprintf("Check status in R: loadRegistry('%s'); getStatus(reg)\n", reg_dir))

# Write job IDs to file for Python to read
writeLines(as.character(ids), file.path(reg_dir, "job_ids.txt"))

# Write registry path for Python
writeLines(reg_dir, file.path(reg_dir, "registry_path.txt"))

