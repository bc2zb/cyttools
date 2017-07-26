#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
cyttools.R (-h | --help | --version)
cyttools.R --makePanelBlank DIR
cyttools.R --makeMetaDataBlank DIR
cyttools.R --computeNRS DIR PANEL
cyttools.R [--cluster=<algorithm>] DIR PANEL
cyttools.R --compDiffAbndnc PANEL FEATURETABLE METADATA
cyttools.R --compDiffExpr PANEL FEATURETABLE METADATA

Description:   This program performs automated high parameter cytometry data analysis.
Options:
--version       Show the current version.
--makePanelBlank         Produce a panel design file based on FCS files in DIR
--makeMetaDataBlank      Produce a meta data file based on FCS files in DIR
--computeNRS             Compute non redundancy score for parameters
--cluster=<algorithm>    [default: FlowType] The algorithm to use for clustering.
--compDiffAbndnc         Test for differential abundance
--compDiffExpr         Test for differential expression

Arguments:

DIR           Provide directory for files to be analyzed
PANEL         Provide a panel design file, use --makePanelBlank to generate and edit as needed
METADATA      Provide a meta data file, use --makeMetaDataBlank to generate and edit as needed
FEATURETABLE  Provide a expression or abundance feature table, these are outputs of --cluster command

" -> doc


args <- docopt(doc)
algList <- c("FlowSOM", "FlowType")


if(args$`--version` == T){ # returns version if version is requested
  cat("\nVersion is pre alpha\n")

}else if(args$`--cluster` %in% algList == F){ # checks that algorithm is available
  cat(paste(c("\nERROR:","\nCurrent algorithms for clustering are:", algList), collapse = "\n"), "\n")  

}else{ # Analysis begins
  
  # create new sub directory for each instance of cyttools
  RESULTS_DIR <- file.path(getwd(), "cyttoolsResults", gsub("\ ", "_", Sys.time()))
  RESULTS_DIR <- paste(RESULTS_DIR, "/", sep = "")
  dir.create(RESULTS_DIR, showWarnings = F, recursive = T)
  argsFileName <- paste(RESULTS_DIR, "cyttools.args", ".Rdata", sep = "")
  save(args, file = argsFileName)
  
  if(args$`--makePanelBlank` == T){
    
    COMMAND <- paste("Rscript MakePanelBlank.R", RESULTS_DIR)
    system(command = COMMAND)
  }else if(args$`--makeMetaDataBlank` == T){
    
    COMMAND <- paste("Rscript makeMetaDataBlank.R", RESULTS_DIR)
    system(command = COMMAND)
  }else if(args$`--computeNRS` == T){
    
    COMMAND <- paste("Rscript NonRedundancyScoreComputation.R", RESULTS_DIR, paste("'", args$PANEL, "'", sep = ""))

    system(command = COMMAND)
  }else if(args$`--compDiffAbndnc` == T){
    
    COMMAND <- paste("Rscript compDiffAbndnc.R",
                     RESULTS_DIR,
                     paste("'", args$PANEL, "'", sep = ""),
                     paste("'", args$FEATURETABLE, "'", sep = ""),
                     paste("'", args$METADATA, "'", sep = ""))
    
    system(command = COMMAND)
  }else if(args$`--compDiffExpr` == T){
    
    COMMAND <- paste("Rscript compDiffExpr.R",
                     RESULTS_DIR,
                     paste("'", args$PANEL, "'", sep = ""),
                     paste("'", args$FEATURETABLE, "'", sep = ""),
                     paste("'", args$METADATA, "'", sep = ""))
    
    system(command = COMMAND)
  }else if(args$`--cluster` == "FlowSOM"){
    
    cat(paste("\nPreparing to run analysis using ", args$`--cluster`, " on ", args$DIR, "\n\n", sep = ""))
    COMMAND <- paste("Rscript FlowSOM.R", RESULTS_DIR, paste("'", args$PANEL, "'", sep = ""))
    system(command = COMMAND)    
  }else if(args$`--cluster` == "FlowType"){
    
    cat(paste("\nPreparing to run analysis using ", args$`--cluster`, " on ", args$DIR, "\n\n", sep = ""))
    COMMAND <- paste("Rscript Flowtype.R", RESULTS_DIR, paste("'", args$PANEL, "'", sep = ""))
    system(command = COMMAND)
  }else{
    cat(paste(c("\nWARNING:","\nAlgorithm not found:", args$`--cluster`), collapse = "\n"), "\n")  
    }
}
