#!/usr/bin/env Rscript
# Script to print the default LSF template used by batchtools
# Usage: Rscript print_lsf_template.R

suppressPackageStartupMessages({
  library(batchtools)
})

cat("=== Default LSF Template (lsf-simple) ===\n\n")

# Try to find and read the template file
template_path <- system.file("templates", "lsf-simple.tmpl", package = "batchtools")

if (template_path != "") {
  cat("Template file location:\n")
  cat(template_path, "\n\n")
  
  cat("Template content:\n")
  cat("---\n")
  template_lines <- readLines(template_path)
  cat(paste(template_lines, collapse = "\n"))
  cat("\n---\n")
} else {
  cat("ERROR: Could not find lsf-simple.tmpl template file.\n")
  cat("Trying alternative method...\n\n")
  
  # Alternative: Check what templates are available
  templates_dir <- system.file("templates", package = "batchtools")
  if (templates_dir != "") {
    cat("Available templates in batchtools package:\n")
    cat(list.files(templates_dir), sep = "\n")
    cat("\n")
    
    # Try to read any LSF template
    lsf_templates <- list.files(templates_dir, pattern = "lsf.*\\.tmpl", full.names = TRUE)
    if (length(lsf_templates) > 0) {
      cat("Found LSF template(s):\n")
      for (tmpl in lsf_templates) {
        cat("\n---", basename(tmpl), "---\n")
        cat(paste(readLines(tmpl), collapse = "\n"))
        cat("\n")
      }
    }
  } else {
    cat("ERROR: Could not find batchtools templates directory.\n")
    cat("batchtools may not be installed or accessible.\n")
  }
}

cat("\n=== Template Variables ===\n")
cat("Batchtools replaces these variables in the template:\n")
cat("  <%= resources %> - Resource requirements (queue, cores, memory, walltime)\n")
cat("  <%= job.name %> - Job name\n")
cat("  <%= log.file %> - Log file path\n")
cat("  <%= rscript %> - R script to execute\n")
cat("  <%= array %> - Array job specification (if applicable)\n")

