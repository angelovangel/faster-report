#!/usr/bin/env Rscript
#============
#
# this just renders the faster-report.Rmd file,
# note that this uses system pandoc for rendering, not the Rstudio one
#  # nolint: trailing_whitespace_linter.
#
#============
require(funr)
require(optparse)
require(rmarkdown)
require(stringr)
require(bslib)
require(bsicons)
require(knitr)
require(shiny)
require(scales)
require(reactable)
require(dplyr)
require(sparkline)
require(htmlwidgets)
require(jsonlite)
require(parallel) # may be ships with R, so not in the environment.yml
require(parallelMap)
#require(renv)


# 1. Parse command-line arguments
option_list <- list(
  make_option(c('--path', '-p'), help = 'path to folder with fastq files', type = 'character', default = NULL),
  make_option(c('--regex', '-r'), help = 'regex pattern to match fastq files', type = 'character', default = 'fast(q|q.gz)$'),
  make_option(c('--type', '-t'), help = "seq platform ('illumina', 'ont', 'pacbio')", default = 'ont'),
  make_option(c('--rundate', '-d'), help = 'Run date', type = 'character', default = NULL),
  make_option(c('--flowcell', '-f'), help = 'Flow cell ID', type = 'character', default = NULL),
  make_option(c('--basecall', '-b'), help = 'Basecaller model', type = 'character', default = NULL),
  make_option(c('--user', '-u'), help = 'User', type = 'character', default = NULL),
  make_option(c('--save_raw', '-s'), help = 'save raw csv data used for plotting', type = 'logical', default = FALSE),
  make_option(c('--subsample', '-x'), help = 'subsample reads for kmers calculation', type = 'double', default = 1.0),
  make_option(c('--outfile','-o'), help = 'name of output report file', type = 'character', default = 'faster-report.html'),
  make_option(c('--git_commit', '-g'), help = 'git commit hash to display in report', type = 'character', default = 'NA')
)

opts <- parse_args(OptionParser(option_list = option_list))

if (is.null(opts$path)) {
  stop("At least a path to a folder with fastq files is required (use option '-p path/to/folder')", call. = FALSE)
}

# 2. Harmonize sequencer platform strings
if (opts$type == 'illumina') {
  opts$type <- 'Illumina'
  opts$basecall <- 'NA'
} else if (opts$type == 'ont') {
  opts$type <- 'Nanopore'
} else if (opts$type == 'pacbio') {
  opts$type <- 'PacBio'
}

# =========================================================================
# 3. UNIVERSAL LOGIC: Track down where the execution binaries live
# =========================================================================
initial_args <- commandArgs(trailingOnly = FALSE)
file_arg <- initial_args[grep("--file=", initial_args)]

if (length(file_arg) > 0) {
  # Resolved script home directory (works on both Docker paths and Singularity mappings)
  script_dir <- dirname(sub("--file=", "", file_arg))
} else {
  script_dir <- getwd()
}

template_filename <- "faster-report.Rmd"
target_source_path <- file.path(script_dir, template_filename)

# Fallback: if the combined directory search fails, assume it's right in the working folder
if (!file.exists(target_source_path)) {
  target_source_path <- template_filename
}

# Read text lines and generate an isolated workspace copy to drop Singularity read-only bugs
local_writable_rmd <- "local-execution-report.Rmd"
writeLines(readLines(target_source_path), local_writable_rmd)

# =========================================================================
# 4. EXPOSE THE CODE DIRECTORY: Pass script_dir directly to the Rmd scope
# =========================================================================
# By setting this globally, you can adjust the execution target paths safely
.GlobalEnv$pipeline_bin_dir <- script_dir

# Render the report safely inside the current writable workspace folder
rmarkdown::render(
  input              = local_writable_rmd, 
  output_file        = opts$outfile,
  output_dir         = getwd(),
  intermediates_dir  = getwd(),
  knit_root_dir      = getwd(),
  params = list(
    fastq_dir     = normalizePath(opts$path, mustWork = FALSE),
    fastq_pattern = opts$regex,
    sequencer     = opts$type,
    rundate       = opts$rundate,
    flowcell      = opts$flowcell,
    basecall      = opts$basecall,
    user          = opts$user,
    rawdata       = opts$save_raw,
    subsample     = opts$subsample,
    git_commit    = opts$git_commit
  )
)

# Clean up our generated local file
if (file.exists(local_writable_rmd)) file.remove(local_writable_rmd)
