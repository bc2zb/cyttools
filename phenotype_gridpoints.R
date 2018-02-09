#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
phenotype_gridpoints.R (-h | --help | --version)
phenotype_gridpoints.R DIR

Description:   This script phenotypes grid points from FlowSOM with FlowType and performs consensus clustering
Options:
--version       Show the current version.

Arguments:

DIR    Provide directory where cyttools.args.Rdata file is located
" -> doc


args <- docopt(doc)

ARGS_DIR <- args$DIR

cat("\nLoading arguments from", ARGS_DIR, "\n")

load(paste(ARGS_DIR, "cyttools.args.Rdata", sep = ""))

RESULTS_DIR <- args$OUT

source("cyttoolsFunctions.R")

##########################################################################
############################ R code goes here ############################
##########################################################################

# read in panel design
targets <- read.delim(args$PANEL)
colsToCheck <- c("Ignore", "TransformCofactor", "Lineage", "Functional", "NRS")
if(checkDesignCols(targets, colsToCheck)){
  missingCols <- colsToCheck[which(colsToCheck %in% colnames(targets) == F)]
  cat("\n\nERROR: PANEL file does not include required columns.
      \n\nMissing Columns:", missingCols,
      "\n\nPlease run cyttools.R --makePanelBlank and cyttools.R --computeNRS to generate compatible panel file.\n\nStopping cyttools.R\n\n")
  q()
}

# read in clusterd FCS files
cluster_dir <- args$CLUSTERDIR # grabs directory from initial cyttools call
file <- list.files(cluster_dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory

lineage_markers <- targets$name[targets$Lineage == 1]
functional_markers <- targets$name[targets$Functional == 1]

if(args$transform == T){
  cluster.flowSet.trans <- read.flowSet.transVS(targets, file, ncdf = F)
}else{
  cluster.flowSet.trans <- read.flowSet(file)
}

##########################################################################
############################     End code     ############################
##########################################################################

workspaceFile <- paste(RESULTS_DIR, "phenotype_gridpoints.Workspace.Rdata", sep = "")

save.image(file = workspaceFile)
