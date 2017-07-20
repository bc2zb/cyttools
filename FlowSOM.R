#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
FLowSOM.R (-h | --help | --version)
FLowSOM.R DIR PANEL

Description:   This script is a template for making docopts compatible Rscripts
Options:
--version       Show the current version.

Arguments:

DIR    Provide directory
PANEL  Provide a panel design file, use --makePanelBlank to generate and edit as needed

" -> doc


args <- docopt(doc)

RESULTS_DIR <- args$DIR

cat("\nLoading arguments from", args$DIR, "\n")

load(paste(RESULTS_DIR, "cyttools.args.Rdata", sep = ""))

library(flowCore)
library(cytofCore)
library(flowVS)
library(flowType)
library(tidyverse)
library(reshape2)
library(FlowSOM)

dir <- args$DIR # grabs directory from initial cyttools call
file <- list.files(dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory
flowSet <- read.flowSet(file, transformation = F) # reads in files as flowSet, required for flowType

targets <- read.delim(args$PANEL)

# if NRS column is missing from the panel, calculate NRS
if(length(grep("NRS", colnames(targets))) == 0){
  
  COMMAND <- paste("Rscript NonRedundancyScoreComputation.R", RESULTS_DIR, paste("'", args$PANEL, "'", sep = ""))
  system(command = COMMAND)
  targets <- read.delim(paste(RESULTS_DIR, "nrsPanelFile.txt", sep = ""))
}

lineage_markers <- targets$name[targets$Lineage == 1]
functional_markers <- targets$name[targets$Functional == 1]

flowSet.trans <- transFlowVS(flowSet,
                             as.character(targets$name[which(targets$Lineage == 1 |
                                                               targets$Functional == 1)]),
                             rep(5, length(targets$name[which(targets$Lineage == 1 |
                                                                targets$Functional == 1)])))


# order the markers using NRS, dropping markers set to "1" in Ignore column of panel design
lineage_markers_ord <- targets$name
lineage_markers_ord <- lineage_markers_ord[lineage_markers_ord %in% targets$name[which(targets$Ignore == 0)]]
lineage_markers_ord <- lineage_markers_ord[order(targets$NRS[which(targets$name %in% lineage_markers_ord)], decreasing = T)]

fsom <- ReadInput(flowSet.trans, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom,
                colsToUse = lineage_markers[lineage_markers %in% targets$name[targets$Ignore == 0]],
                xdim = length(lineage_markers[lineage_markers %in% targets$name[targets$Ignore == 0]]),
                ydim = length(lineage_markers[lineage_markers %in% targets$name[targets$Ignore == 0]]))

flowSOM.res <- BuildMST(som)



workspaceFile <- paste(RESULTS_DIR, "FlowSOMWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)
