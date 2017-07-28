#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
cyttools.R (-h | --help | --version)
cyttools.R --makePanelBlank DIR OUT
cyttools.R --makeMetaDataBlank DIR OUT
cyttools.R [--transform=<bool>] --computeNRS DIR PANEL OUT
cyttools.R [--transform=<bool>] --cluster=<algorithm> DIR PANEL OUT
cyttools.R --compDiffAbndnc PANEL FEATURETABLE METADATA OUT
cyttools.R --compDiffExpr PANEL FEATURETABLE METADATA OUT

Description:   This program performs automated high parameter cytometry data analysis.
Options:
--version                Show the current version.
--makePanelBlank         Produce a panel design file based on FCS files in DIR
--makeMetaDataBlank      Produce a meta data file based on FCS files in DIR
--transform=<true>       [default: TRUE] Transform data using arcsinH and cofactors specified in PANEL
--computeNRS             Compute non redundancy score for parameters
--cluster=<algorithm>    [default: FlowSOM] The algorithm to use for clustering.
--compDiffAbndnc         Test for differential abundance
--compDiffExpr           Test for differential expression

Arguments:

DIR           Provide directory for files to be analyzed
OUT           Provide directory for results to be saved to
PANEL         Provide a panel design file, use --makePanelBlank to generate and edit as needed
METADATA      Provide a meta data file, use --makeMetaDataBlank to generate and edit as needed
FEATURETABLE  Provide a expression or abundance feature table, these are outputs of --cluster command

" -> doc


args <- docopt(doc)
algList <- c("FlowSOM",
             "FlowType", "ParFlowType",
             "BatchFlowTypeDataPrep", "BatchFlowType", "BatchFlowTypeDataMerge")


if(args$`--version` == T){ # returns version if version is requested
  cat("\nVersion is 0.2\n")

}else if(args$`--cluster` %in% algList == F){ # checks that algorithm is available
  cat(paste(c("\nERROR:","\nCurrent algorithms for clustering are:", algList), collapse = "\n"), "\n")  

}else{ # Analysis begins
  
  # create new sub directory for each instance of cyttools
  RESULTS_DIR <- args$OUT
  dir.create(RESULTS_DIR, showWarnings = F, recursive = T)
  argsFileName <- paste(RESULTS_DIR, "cyttools.args", ".Rdata", sep = "")
  save(args, file = argsFileName)
  
  if(args$`--makePanelBlank` == T){
    
    COMMAND <- paste("Rscript MakePanelBlank.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--makeMetaDataBlank` == T){
    
    COMMAND <- paste("Rscript MakeMetaDataBlank.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--computeNRS` == T){
    
    COMMAND <- paste("Rscript NonRedundancyScoreComputation.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--compDiffAbndnc` == T){
    
    COMMAND <- paste("Rscript compDiffAbndnc.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--compDiffExpr` == T){
    
    COMMAND <- paste("Rscript compDiffExpr.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--cluster` == "FlowSOM"){
    
    cat(paste("\nPreparing to run analysis using ", args$`--cluster`, " on ", args$DIR, "\n\n", sep = ""))
    COMMAND <- paste("Rscript FlowSOM.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)    
  }else if(args$`--cluster` == "FlowType"){
    
    cat(paste("\nPreparing to run analysis using ", args$`--cluster`, " on ", args$DIR, "\n\n", sep = ""))
    COMMAND <- paste("Rscript Flowtype.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--cluster` == "ParFlowType"){
    
    cat(paste("\nPreparing to run analysis using ", args$`--cluster`, " on ", args$DIR, "\n\n", sep = ""))
    COMMAND <- paste("Rscript ParFlowtype.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--cluster` == "BatchFlowType"){
    
    cat(paste("\nPreparing to run analysis using ", args$`--cluster`, " on ", args$DIR, "\n\n", sep = ""))
    COMMAND <- paste("Rscript BatchFlowtype.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--cluster` == "BatchFlowTypeDataPrep"){
    
    cat(paste("\nPreparing to run analysis using ", args$`--cluster`, " on ", args$DIR, "\n\n", sep = ""))
    COMMAND <- paste("Rscript BatchFlowTypeDataPrep.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else if(args$`--cluster` == "BatchFlowTypeDataMerge"){
    
    cat(paste("\nPreparing to run analysis using ", args$`--cluster`, " on ", args$DIR, "\n\n", sep = ""))
    COMMAND <- paste("Rscript BatchFlowTypeDataMerge.R",
                     paste("'", RESULTS_DIR, "'", sep = ""))
    system(command = COMMAND)
  }else{
    cat(paste(c("\nWARNING:","\nAlgorithm not found:", args$`--cluster`), collapse = "\n"), "\n")  
    }
}
