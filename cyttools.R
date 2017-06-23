#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
cyttools.R [--number=<number>] [--min=<min>] [--max=<max>]
cyttools.R (-h | --help | --version)
cyttools.R DIR

Description:   This program performs automated high parameter cytometry data analysis.
Options:
--version       Show the current version.
--number=<num>  [default: 1] The number of random numbers to generate.
--min=<min>     [default: 0] The lowest value a random number can have.
--max=<max>     [default: 1] The highest value a random number can have.

Arguements:

DIR    Provide directory for files to be analyzed.
" -> doc


args    <- docopt(doc)

# returns version if version is requested
if(args$`--version` == T){
  cat("Version is pre alpha\n")
}else if(is.null(args$DIR)){ # checks for directory provided
  cat("\nERROR: Please provide a directory or file\n")
  
}else{
  cat(paste("\nRunning analysis on ", args$DIR, "\n", sep = ""))
}