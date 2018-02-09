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

# your R code!

##########################################################################
############################     End code     ############################
##########################################################################

workspaceFile <- paste(RESULTS_DIR, "phenotype_gridpoints.Workspace.Rdata", sep = "")

save.image(file = workspaceFile)
