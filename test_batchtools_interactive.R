#!/usr/bin/env Rscript
# Interactive test script for batchtools LSF submission
# Usage: Rscript test_batchtools_interactive.R

suppressPackageStartupMessages({
  library(batchtools)
  library(jsonlite)
})

cat("=== Testing batchtools LSF Configuration ===\n\n")

# Test 1: Check if batchtools is available
cat("1. Testing batchtools package...\n")
if (!requireNamespace("batchtools", quietly = TRUE)) {
  stop("ERROR: batchtools package not installed. Install with: install.packages('batchtools')")
}
cat("   ✓ batchtools package loaded\n\n")

# Test 2: Check LSF commands
cat("2. Testing LSF commands...\n")
test_bqueues <- system("bqueues > /dev/null 2>&1", intern = FALSE)
test_bsub <- system("which bsub > /dev/null 2>&1", intern = FALSE)

if (test_bqueues == 0) {
  cat("   ✓ bqueues command available\n")
  # Show available queues
  queues <- system("bqueues | tail -n +2 | awk '{print $1}'", intern = TRUE)
  cat("   Available queues:", paste(queues, collapse = ", "), "\n")
} else {
  cat("   ✗ bqueues command NOT available (LSF may not be configured)\n")
}

if (test_bsub == 0) {
  cat("   ✓ bsub command available\n")
} else {
  cat("   ✗ bsub command NOT available\n")
}
cat("\n")

# Test 3: Check template file
cat("3. Testing template file...\n")
script_dir <- dirname(normalizePath(commandArgs()[4]))
template_file <- file.path(script_dir, "lsf_custom.tmpl")

if (file.exists(template_file)) {
  cat("   ✓ Template file found:", template_file, "\n")
  template_content <- readLines(template_file, n = 5)
  cat("   Template preview (first 5 lines):\n")
  for (line in template_content) {
    cat("     ", line, "\n")
  }
} else {
  cat("   ✗ Template file NOT found:", template_file, "\n")
  cat("   Will use built-in 'lsf-simple' template\n")
}
cat("\n")

# Test 4: Create a test registry
cat("4. Testing registry creation...\n")
test_reg_dir <- file.path(tempdir(), "batchtools_test_registry")
if (dir.exists(test_reg_dir)) {
  unlink(test_reg_dir, recursive = TRUE)
}

tryCatch({
  reg <- makeRegistry(
    file.dir = test_reg_dir,
    conf.file = NA,
    work.dir = getwd(),
    seed = 1
  )
  cat("   ✓ Registry created successfully\n")
  
  # Test 5: Configure cluster functions
  cat("\n5. Testing LSF cluster functions configuration...\n")
  if (file.exists(template_file)) {
    reg$cluster.functions <- makeClusterFunctionsLSF(
      template = template_file,
      scheduler.latency = 1,
      fs.latency = 65
    )
    cat("   ✓ LSF cluster functions configured with custom template\n")
  } else {
    reg$cluster.functions <- makeClusterFunctionsLSF(
      template = "lsf-simple",
      scheduler.latency = 1,
      fs.latency = 65
    )
    cat("   ✓ LSF cluster functions configured with built-in template\n")
  }
  
  if (is.null(reg$cluster.functions)) {
    stop("ERROR: Cluster functions are NULL!")
  }
  cat("   Cluster functions type:", class(reg$cluster.functions)[1], "\n")
  
  # Test 6: Create a simple test job
  cat("\n6. Testing job creation and submission...\n")
  test_fun <- function(x) {
    Sys.sleep(1)
    return(x * 2)
  }
  
  # Submit a single test job
  ids <- batchMap(fun = test_fun, x = 1:1, reg = reg)
  cat("   ✓ Test job created (ID:", ids, ")\n")
  
  # Prepare resources
  resources <- list(
    ncpus = 1,
    mpp = "1GB",
    queue = "normal",  # Use a common queue name - adjust if needed
    walltime = 300,   # 5 minutes
    job.name = "test_job"
  )
  
  cat("   Submitting test job with resources:\n")
  cat("     Queue:", resources$queue, "\n")
  cat("     Cores:", resources$ncpus, "\n")
  cat("     Memory:", resources$mpp, "\n")
  cat("     Walltime:", resources$walltime, "seconds\n")
  
  # Try to submit
  tryCatch({
    submitJobs(ids, resources = resources, reg = reg)
    cat("   ✓ Job submitted successfully\n")
    
    # Wait a moment for LSF to register
    Sys.sleep(3)
    
    # Check job status
    job_table <- getJobTable(ids = ids, reg = reg)
    cat("\n   Job status:\n")
    print(job_table[, c("job.id", "batch.id", "submitted", "started", "done", "error", "running", "queued")])
    
    if (!is.na(job_table$batch.id[1])) {
      lsf_job_id <- job_table$batch.id[1]
      cat("\n   ✓ LSF batch ID found:", lsf_job_id, "\n")
      cat("   Verify with: bjobs", lsf_job_id, "\n")
      
      # Try to query LSF directly
      bjobs_cmd <- paste("bjobs", lsf_job_id)
      bjobs_result <- system(bjobs_cmd, intern = TRUE)
      if (length(bjobs_result) > 1) {
        cat("   ✓ Job visible in LSF:\n")
        for (line in bjobs_result) {
          cat("     ", line, "\n")
        }
      } else {
        cat("   ⚠ Job not found in bjobs output (may have finished already)\n")
      }
    } else {
      cat("\n   ✗ No LSF batch ID found! Job was not submitted to LSF.\n")
      cat("   Check batchtools logs in:", test_reg_dir, "/logs/\n")
    }
    
  }, error = function(e) {
    cat("   ✗ ERROR submitting job:", as.character(e), "\n")
    cat("   Check:\n")
    cat("     1. LSF queue name is correct (currently:", resources$queue, ")\n")
    cat("     2. You have permission to submit jobs\n")
    cat("     3. Template file is valid\n")
  })
  
  # Cleanup
  cat("\n7. Cleaning up test registry...\n")
  unlink(test_reg_dir, recursive = TRUE)
  cat("   ✓ Test complete\n")
  
}, error = function(e) {
  cat("   ✗ ERROR:", as.character(e), "\n")
  if (dir.exists(test_reg_dir)) {
    cat("   Registry directory:", test_reg_dir, "\n")
    cat("   Check logs in:", file.path(test_reg_dir, "logs"), "\n")
  }
})

cat("\n=== Test Summary ===\n")
cat("If all tests passed, batchtools should be able to submit jobs to LSF.\n")
cat("If jobs don't have LSF batch IDs, check:\n")
cat("  1. LSF queue name (use: bqueues to see available queues)\n")
cat("  2. Template file syntax\n")
cat("  3. Batchtools logs in registry directory\n")
cat("  4. LSF permissions and configuration\n")

