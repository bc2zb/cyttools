#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
cyttools.R [--number=<number>] [--min=<min>] [--max=<max>]
cyttools.R (-h | --help | --version)
cyttools.R DIR
cyttools.R --cluster=<algorithm>

Description:   This program performs automated high parameter cytometry data analysis.
Options:
--version       Show the current version.
--number=<num>  [default: 1] The number of random numbers to generate.
--min=<min>     [default: 0] The lowest value a random number can have.
--max=<max>     [default: 1] The highest value a random number can have.
--cluster=<algorithm>    [default: FlowSOM] The algorithm to use for clustering.

Arguements:

DIR    Provide directory for files to be analyzed.
" -> doc


args    <- docopt(doc)
algList <- c("FlowSOM", "FlowType")

# returns version if version is requested
if(args$`--version` == T){
  cat("Version is pre alpha\n")
}else if(is.null(args$DIR)){ # checks for directory provided
  cat("\nERROR: Please provide a directory or file\n")
  
}else if(args$`--cluster` %in% algList){
  cat(paste("\nPreparing to run analysis using ", args$`--cluster`, "\n", sep = ""))
}else if(args$`--cluster` %in% algList == F){
  cat(paste(c("\nCurrent algorithms for clustering are:", algList), collapse = "\n"), "\n")  
}else{
  cat(paste("\nRunning analysis on ", args$DIR, "\n", sep = ""))
}