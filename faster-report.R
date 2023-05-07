#!/usr/bin/env Rscript
#============
#
# this just renders the faster-report.Rmd file,
# note that this uses system pandoc for rendering, not the Rstudio one
#
#============

library(optparse)
require(rmarkdown)

option_list <- list(
  make_option(c('--path', '-p'), help = 'path to folder with fastq files [%default]', type = 'character', default = NULL),
  make_option(c('--regex', '-r'), help = 'regex pattern to match fastq files [%default]', type = 'character', default = '*.fastq'),
  make_option(c('--type', '-t'), help = "seq platform used, can be one of 'illumina', 'ont' or 'pacbio' [%default]", default = 'ont'),
  make_option(c('--save_raw', '-s'), help = 'save raw csv data used for plotting [%default]', type = 'logical', default = FALSE)
  )

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

if (is.null(opts$path)){
  print_help(opt_parser)
  stop("At least a path to a folder with fastq files is required (use option '-p path/to/folder')", call.=FALSE)
}

# change to match parameter. used in Rmd
if (opts$type == 'illumina') {
  opts$type <- 'Illumina'
}

# render the rmarkdown, using fastq-report.Rmd as template
rmarkdown::render(input = "faster-report.Rmd",
                  output_file = "faster-report.html",
                  output_dir = getwd(), # important when knitting in docker
                  knit_root_dir = getwd(), # important when knitting in docker
                  params = list(
                    fastq_dir = opts$path,
                    fastq_pattern = opts$regex,
                    sequencer = opts$type,
                    rawdata = opts$save_raw
                  )
)
