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
--cluster=<algorithm>    [default: FlowType] The algorithm to use for clustering.

Arguements:

DIR    Provide directory for files to be analyzed.
" -> doc


args <- docopt(doc)
algList <- c("FlowSOM", "FlowType")


if(args$`--version` == T){ # returns version if version is requested
  cat("\nVersion is pre alpha\n")

}else if(is.null(args$DIR)){ # checks for directory provided TBD: Change to check for FCS files in the directory
  cat("\nERROR:\nPlease provide a directory or file\n")

}else if(args$`--cluster` %in% algList == F){ # checks that algorithm is available
  cat(paste(c("\nERROR:","\nCurrent algorithms for clustering are:", algList), collapse = "\n"), "\n")  

}else{ # Analysis begins
  
  # create new sub directory for each instance of cyttools
  RESULTS_DIR <- file.path(getwd(), "cyttoolsResults", gsub("\ ", "_", Sys.time()))
  RESULTS_DIR <- paste(RESULTS_DIR, "/", sep = "")
  dir.create(RESULTS_DIR, showWarnings = F, recursive = T)
  argsFileName <- paste(RESULTS_DIR, "cyttools.args", ".Rdata", sep = "")
  save(args, file = argsFileName)
  
  cat(paste("\nPreparing to run analysis using ", args$`--cluster`, " on ", args$DIR, "\n", sep = ""))
  
  if(args$`--cluster` == "FlowType"){
<<<<<<< HEAD
    COMMAND <- paste("Rscript Flowtype.R", RESULTS_DIR)
    #COMMAND <- "Rscript Flowtype.R -h"
    system(command = COMMAND)
  }else{
    cat(paste(c("\nWARNING:","\nAlgorithm not found:", args$`--cluster`), collapse = "\n"), "\n")  
    }
}
=======
    source("FlowType.R")
  }else{
    cat(paste(c("\nWARNING:","\nAlgorithm not found:", args$`--cluster`), collapse = "\n"), "\n")  
    }
}

# create new sub directory for each instance of cyttools

dir.create(file.path(getwd(), "cyttoolsResults"), showWarnings = FALSE)
argsFileName <- paste(file.path(getwd(), "cyttoolsResults"), "/", Sys.time(), "_", "cyttools.args", ".Rdata", sep = "")
save(args, file = argsFileName)
>>>>>>> 0b8fab20884663e8ee8d43a68887adba1c038b7c
