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

calldir <- getwd()
scriptdir  <-  dirname(funr::sys.script())
#setwd(scriptdir)
#renv::load()

option_list <- list(
  make_option(c('--path', '-p'), help = 'path to folder with fastq files [%default]', type = 'character', default = NULL),
  make_option(c('--regex', '-r'), help = 'regex pattern to match fastq files [%default]', type = 'character', default = 'fast(q|q.gz)$'),
  make_option(c('--type', '-t'), help = "seq platform used, can be one of 'illumina', 'ont' or 'pacbio' [%default]", default = 'ont'),
  make_option(c('--rundate', '-d'), help = 'Run date', type = 'character', default = NULL),
  make_option(c('--flowcell', '-f'), help = 'Flow cell ID', type = 'character', default = NULL),
  make_option(c('--basecall', '-b'), help = 'Basecaller model', type = 'character', default = NULL),
  make_option(c('--user', '-u'), help = 'User', type = 'character', default = NULL),
  make_option(c('--save_raw', '-s'), help = 'save raw csv data used for plotting [%default]', type = 'logical', default = FALSE),
  make_option(c('--subsample', '-x'), help = 'subsample reads for kmers calculation [%default]', type = 'double', default = 1.0),
  make_option(c('--outfile','-o'), help = 'name of output report file [%default]', type = 'character', default = 'faster-report.html'),
  make_option(c('--git_commit', '-g'), help = 'git commit hash to display in report [%default]', type = 'character', default = 'NA')
  )

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

if (is.null(opts$path)){
  print_help(opt_parser)
  stop("At least a path to a folder with fastq files is required (use option '-p path/to/folder')", call.=FALSE)
}

# complicated case to parse correct fastq path when calling and script directories are not the same
# check if abs or relative path was provided
if (R.utils::isAbsolutePath(opts$path)) {
  fastqpath <- opts$path
} else {
  fastqpath <- normalizePath(file.path(calldir, opts$path))
}

print(paste0("call dir: ", calldir))
print(paste0("fastq path: ", fastqpath))

# change to match parameter used in Rmd
if (opts$type == 'illumina') {
  opts$type <- 'Illumina'
  opts$basecall <- 'NA'
} else if (opts$type == 'ont') {
  opts$type <- 'Nanopore'
} else if (opts$type == 'pacbio') {
  opts$type <- 'PacBio'
}
# Set your input template and exact target output paths
rmd_template_path <- file.path(scriptdir, "faster-report.Rmd")
tmp_md             <- file.path(calldir, "faster-report.knit.md")
final_output       <- file.path(calldir, opts$outfile)

# 2. Re-create your parameters environment
# This allows your Rmd file to use 'params$variable' exactly as it did before
knit_env <- new.env(parent = globalenv())
knit_env$params <- list(
  fastq_dir     = fastqpath,
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

# ==========================================
# STAGE 1: Process code chunks & knit to MD
# ==========================================
# We explicitly define 'output = tmp_md' to force it into the writable calldir
knitr::knit(
  input = rmd_template_path, 
  output = tmp_md, 
  envir = knit_env
)

# ==========================================
# STAGE 2: Convert MD to final report format
# ==========================================
# This reads the newly created markdown file and exports the finished document
rmarkdown::pandoc_convert(
  input = tmp_md, 
  to = "html", # 
  output = final_output, 
  options = c("--self-contained") # Keeps images/styles embedded inside one file
)

# Clean up the intermediate file in your directory
if (file.exists(tmp_md)) {
  file.remove(tmp_md)
}
#  Copy the Rmd file from the read-only script dir to your writable work dir
# local_rmd <- file.path(calldir, "faster-report.Rmd")
# file.copy(from = file.path(scriptdir, "faster-report.Rmd"), to = local_rmd, overwrite = TRUE)

# # render the rmarkdown, using fastq-report.Rmd as template
# rmarkdown::render(input = local_rmd,
#                   output_file = opts$outfile,
#                   #output_dir = calldir, # important when knitting in docker
#                   #intermediates_dir = calldir, # important when knitting in docker
#                   #knit_root_dir = calldir, # important when knitting in docker
#                   #envir = new.env(),
#                   params = list(
#                     fastq_dir = fastqpath,
#                     fastq_pattern = opts$regex,
#                     sequencer = opts$type,
#                     rundate = opts$rundate,
#                     flowcell = opts$flowcell,
#                     basecall = opts$basecall,
#                     user = opts$user,
#                     rawdata = opts$save_raw,
#                     subsample = opts$subsample,
#                     git_commit = opts$git_commit
#                   )
# )
