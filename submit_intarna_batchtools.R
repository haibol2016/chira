#!/usr/bin/env Rscript
# R script to submit IntaRNA hybridization jobs using batchtools (LSF cluster)
# Usage: Rscript submit_intarna_batchtools.R <config_json> <jobs_json>
#
# Follows the same pattern as submit_chunks_batchtools.R (chira_map.py).
# Each job runs IntaRNA once on a pairwise FASTA (query.fa + target.fa) for one chunk.

suppressPackageStartupMessages({
  library(batchtools)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript submit_intarna_batchtools.R <config_json> <jobs_json>")
}

config_file <- args[1]
jobs_file <- args[2]

if (!file.exists(config_file)) {
  stop(sprintf("ERROR: Config file does not exist: %s", config_file))
}
if (!file.exists(jobs_file)) {
  stop(sprintf("ERROR: Jobs file does not exist: %s", jobs_file))
}

tryCatch({
  config <- fromJSON(config_file)
}, error = function(e) {
  cat(sprintf("ERROR: Failed to parse config JSON: %s\n", config_file))
  cat(sprintf("Error: %s\n", as.character(e)))
  stop(e)
})

tryCatch({
  jobs <- fromJSON(jobs_file)
}, error = function(e) {
  cat(sprintf("ERROR: Failed to parse jobs JSON: %s\n", jobs_file))
  cat(sprintf("Error: %s\n", as.character(e)))
  stop(e)
})

reg_dir <- config$reg_dir
queue <- config$queue
cores_per_job <- config$cores_per_job
memory_per_job <- config$memory_per_job
walltime <- config$walltime
conda_env <- config$conda_env
intarna_params <- config$intarna_params
job_name_prefix <- config$job_name_prefix
max_parallel <- ifelse(is.null(config$max_parallel), nrow(jobs), config$max_parallel)
template_file <- if (is.null(config$template_file) || config$template_file == "") {
  if (file.exists("lsf_custom.tmpl")) {
    normalizePath("lsf_custom.tmpl")
  } else {
    "lsf-simple"
  }
} else {
  config$template_file
}

if (dir.exists(reg_dir)) {
  unlink(reg_dir, recursive = TRUE, force = TRUE)
}

reg <- makeRegistry(
  file.dir = reg_dir,
  conf.file = NA,
  work.dir = getwd(),
  seed = 1
)

if (template_file == "lsf-simple" || file.exists(template_file)) {
  cat(sprintf("Configuring LSF cluster functions with template: %s\n", template_file))
  reg$cluster.functions <- makeClusterFunctionsLSF(
    template = template_file,
    scheduler.latency = 1,
    fs.latency = 65
  )
  if (is.null(reg$cluster.functions)) {
    stop("ERROR: Cluster functions are NULL after configuration!")
  }
} else {
  stop(sprintf("Template file '%s' does not exist. Use --batchtools_template to specify a valid template.", template_file))
}

run_intarna_job <- function(n, query_fa, target_fa, output_csv, params) {
  args <- c(strsplit(params, " ")[[1]], "-q", query_fa, "-t", target_fa, "--outPairwise", "--out", output_csv)
  args <- args[args != ""]
  ret <- system2("IntaRNA", args = args, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(ret, "status")) && attr(ret, "status") != 0) {
    stop(paste(ret, collapse = "\n"))
  }
  invisible(NULL)
}

cat(sprintf("Submitting %d IntaRNA jobs to LSF cluster...\n", nrow(jobs)))

walltime_seconds <- if (grepl(":", walltime)) {
  parts <- strsplit(walltime, ":")[[1]]
  as.numeric(parts[1]) * 3600 + as.numeric(parts[2]) * 60
} else {
  as.numeric(walltime) * 60
}

resources <- list(
  ncpus = cores_per_job,
  mpp = memory_per_job,
  queue = queue,
  walltime = walltime_seconds,
  job.name = paste0(job_name_prefix, "_intarna")
)

ids <- batchMap(
  fun = run_intarna_job,
  n = jobs$n,
  query_fa = jobs$query_fa,
  target_fa = jobs$target_fa,
  output_csv = jobs$output_csv,
  more.args = list(params = intarna_params),
  reg = reg
)

if (max_parallel < length(ids)) {
  num_batches <- ceiling(length(ids) / max_parallel)
  remaining_ids <- ids
  all_submitted <- c()

  for (batch_idx in 1:num_batches) {
    batch_size <- min(max_parallel, length(remaining_ids))
    batch_ids <- remaining_ids[1:batch_size]
    remaining_ids <- remaining_ids[-(1:batch_size)]

    cat(sprintf("Submitting batch %d/%d: %d jobs...\n", batch_idx, num_batches, length(batch_ids)))
    tryCatch({
      submitJobs(batch_ids, resources = resources, reg = reg)
      all_submitted <- c(all_submitted, batch_ids)
    }, error = function(e) {
      cat(sprintf("ERROR submitting batch %d: %s\n", batch_idx, as.character(e)))
      stop(e)
    })

    if (length(remaining_ids) > 0) {
      cat(sprintf("Waiting for batch %d to complete...\n", batch_idx))
      while (!waitForJobs(ids = batch_ids, sleep = 30, timeout = Inf, stop.on.error = FALSE, reg = reg)) {
        job_table <- getJobTable(ids = batch_ids, reg = reg)
        cat(sprintf("  %d done, %d running, %d error. Waiting...\n",
          sum(job_table$done, na.rm = TRUE),
          sum(job_table$running, na.rm = TRUE),
          sum(job_table$error, na.rm = TRUE)))
        Sys.sleep(30)
      }
    }
  }
  submitted_ids <- all_submitted
} else {
  cat(sprintf("Submitting all %d jobs at once...\n", length(ids)))
  tryCatch({
    submitJobs(ids, resources = resources, reg = reg)
    submitted_ids <- ids
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", as.character(e)))
    stop(e)
  })

  cat("Waiting for all jobs to complete...\n")
  waitForJobs(ids = ids, sleep = 30, timeout = Inf, stop.on.error = FALSE, reg = reg)
}

job_table <- getJobTable(ids = submitted_ids, reg = reg)
lsf_ids <- job_table$batch.id[!is.na(job_table$batch.id)]
cat(sprintf("\nRegistry: %s\n", reg_dir))
cat(sprintf("LSF job IDs: %s\n", ifelse(length(lsf_ids) > 0, paste(lsf_ids, collapse = ", "), "NONE")))

writeLines(as.character(ids), file.path(reg_dir, "job_ids.txt"))
writeLines(reg_dir, file.path(reg_dir, "registry_path.txt"))
