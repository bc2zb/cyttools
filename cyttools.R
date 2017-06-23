#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
cyttools.R (-h | --help | --version)
cyttools.R [--cluster=<algorithm>] DIR

Description:   This program performs automated high parameter cytometry data analysis.
Options:
--version       Show the current version.
--cluster=<algorithm>    [default: FlowSOM] The algorithm to use for clustering.

Arguements:

DIR    Provide directory for files to be analyzed.
" -> doc


args    <- docopt(doc)
algList <- c("FlowSOM", "FlowType")

# returns version if version is requested
if(args$`--version` == T){
  cat("\nVersion is pre alpha\n")

}else if(is.null(args$DIR)){ # checks for directory provided
  cat("\nERROR:\nPlease provide a directory or file\n")

}else if(args$`--cluster` %in% algList == F){ # checks that algorithm is available
  cat(paste(c("\nERROR:","\nCurrent algorithms for clustering are:", algList), collapse = "\n"), "\n")  

# ADD CHECK FOR EXISTANCE OF FCS FILES PRINT NUMBER OF FCS FILES IN DIRECTORY
}else{ # Analysis begins
  
  cat(paste("\nPreparing to run analysis using ", args$`--cluster`, " on ", args$DIR, "\n", sep = ""))
}