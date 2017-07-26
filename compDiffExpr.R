#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
compDiffExpr.R (-h | --help | --version)
compDiffExpr.R DIR PANEL FEATURETABLE METADATA

Description:   This script calculates differential abundance using a feature table
Options:
--version       Show the current version.

Arguments:

DIR             Provide directory
PANEL           Provide a panel design file
FEATURETABLE    Provide expression feature table file
METADATA        Provide meta data file
" -> doc


args <- docopt(doc)

RESULTS_DIR <- args$DIR

cat("\nLoading arguments from", args$DIR, "\n")

load(paste(RESULTS_DIR, "cyttools.args.Rdata", sep = ""))

library(limma)
library(tidyverse)
library(reshape2)

logitTransform <- function(p) { log(p/(1-p)) }
panelDesign <- read.delim(args$PANEL)
md <- read.delim(args$METADATA)
nodeExprTable <- read.delim(args$FEATURETABLE, row.names = 1)

targets <- md

targets$TimePoint <- factor(targets$TimePoint, levels = unique(targets$TimePoint))
targets$Condition <- factor(targets$Condition, levels = unique(targets$Condition))
targets$SampleID <- factor(targets$SampleID, levels = unique(targets$SampleID))

exprDesign <- paste(targets$Condition, targets$TimePoint, sep = "_")

design <- model.matrix(~ 0 + exprDesign + targets$SampleID + targets$Group)

colnames(design) <- gsub("exprDesign", "Cnd", colnames(design))
colnames(design) <- gsub("targets\\$SampleID|targets\\$Group", "BatchEffect", colnames(design))

# Automate generation of contrast matrix
cont.matrix <- matrix(nrow = length(colnames(design)), ncol = (length(colnames(design)) - 1)*(length(colnames(design)))/2) # levels X Contrasts
prevNumContrasts <- 1
for ( i in 1:nrow(cont.matrix)){
  numContrasts <- nrow(cont.matrix) - i
  endContrasts <- prevNumContrasts + numContrasts - 1
  if (numContrasts > 0){
    cont.matrix[i,c(prevNumContrasts:endContrasts)] <- 1
  }
  k <- i
  l <- i + numContrasts - 1
  if ( k <= l){
    for ( j in k:l){
      cont.matrix[j+1,c(prevNumContrasts:endContrasts)[j - i + 1]] <- -1      
    }
  }
  prevNumContrasts <- prevNumContrasts + numContrasts
  
}  
cont.matrix[is.na(cont.matrix)] <- 0
row.names(cont.matrix) <- colnames(design)
colnames(cont.matrix) <- paste("Contrast", c(1:ncol(cont.matrix)), sep = "")

# remove comparisons that don't make sense from an experimental design perspective
contMatrixRowIDs <- colsplit(row.names(cont.matrix), "\\_", c("Condition", "TimePoint"))

validContrastIndex <- vector(length = ncol(cont.matrix))
humanReadableColNames <- vector(length = ncol(cont.matrix))
for ( i in 1:ncol(cont.matrix)){
  contrastColumn <- cont.matrix[,i]
  firstCondition <- contMatrixRowIDs$Condition[cont.matrix[,i] == 1]
  secondCondition <- contMatrixRowIDs$Condition[cont.matrix[,i] == -1]
  firstTimePoint <- contMatrixRowIDs$TimePoint[cont.matrix[,i] == 1]
  secondTimePoint <- contMatrixRowIDs$TimePoint[cont.matrix[,i] == -1]
  
  if (firstCondition == secondCondition & firstTimePoint != secondTimePoint){
    validContrastIndex[i] <- T
    humanReadableColNames[i] <- paste(firstCondition,
                                      paste(firstTimePoint,
                                            secondTimePoint,
                                            sep = "_vs_"),
                                      sep = ".")
  }else if(firstCondition != secondCondition & firstTimePoint == secondTimePoint){
    validContrastIndex[i] <- T
    humanReadableColNames[i] <- paste(firstTimePoint,
                                      paste(firstCondition,
                                            secondCondition,
                                            sep = "_vs_"),
                                      sep = ".")
  }else{validContrastIndex[i] <- F}
}

cont.matrix <- cont.matrix[,validContrastIndex]
colnames(cont.matrix) <- humanReadableColNames[validContrastIndex]

cont.matrix <- cont.matrix[,grep("BatchEffect", colnames(cont.matrix), invert = T)]

fit <- lmFit(nodeExprTable, design = design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

diffExprStatsTable <- data.frame()

for ( i in 1:length(colnames(fit2))){
  nextPart <- topTable(fit2, coef = i, number = Inf)
  nextPart$Source <- rep(colnames(fit2)[i], nrow(nextPart))
  nextPart$RowNames <- row.names(nextPart)
  nextPart <- separate(nextPart,
                       RowNames,
                       c("Mapping", "Metal"),
                       "_")
  diffExprStatsTable <- rbind(diffExprStatsTable, nextPart)
}

nodeExprStatsFile <- paste(RESULTS_DIR, "nodeDifferentialExpressionTable.txt", sep = "")

write.table(diffExprStatsTable, nodeExprStatsFile, sep = "\t", quote = F, row.names = F)

workspaceFile <- paste(RESULTS_DIR, "compDiffAbndncWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)