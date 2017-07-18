#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
NonRedundancyScoreComputation.R (-h | --help | --version)
NonRedundancyScoreComputation.R DIR PANEL

Description:   This script is a template for making docopts compatible Rscripts
Options:
--version       Show the current version.

Arguments:

DIR    Provide directory
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

dir <- args$DIR # grabs directory from initial cyttools call
file <- list.files(dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory
flowSet <- read.flowSet(file, transformation = F) # reads in files as flowSet, required for flowType

targets <- read.delim(args$PANEL)

lineage_markers <- targets$name[targets$Lineage == 1]
functional_markers <- targets$name[targets$Functional == 1]

flowSet.trans <- transFlowVS(flowSet,
                             as.character(targets$name[which(targets$Lineage == 1 |
                                                               targets$Functional == 1)]),
                             rep(5, length(targets$name[which(targets$Lineage == 1 |
                                                                targets$Functional == 1)])))

## Define a function that calculates the NRS per sample
NRS <- function(x, ncomp = 3){
  pr <- prcomp(x, center = TRUE, scale. = FALSE)
  score <- rowSums(outer(rep(1, ncol(x)),
                         pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
  return(score)
}

## Calculate the score
nrs_sample <- fsApply(flowSet.trans[,colnames(flowSet.trans) %in% lineage_markers], NRS, use.exprs = TRUE)

nrs <- colMeans(nrs_sample, na.rm = TRUE)

## Plot the NRS for ordered markers
lineage_markers_ord <- names(sort(nrs, decreasing = TRUE))
nrs_sample <- data.frame(nrs_sample)
nrs_sample$sample_id <- rownames(nrs_sample)

nrs_all <- rep(NA, nrow(targets))
for(i in 1:length(nrs_all)){
  if(length(which(names(nrs) == targets$name[i])) == 0)
    nrs_all[i] <- NA
  else{
    nrs_all[i] <- nrs[which(names(nrs) == targets$name[i])]
  }
}

targets$NRS <- nrs_all

nrsFile <- paste(RESULTS_DIR, "nrsPanelFile.txt", sep = "")

write.table(targets, file = nrsFile, sep = "\t", quote = F, row.names = F)

workspaceFile <- paste(RESULTS_DIR, "NonRedundancyScoreComputationWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)