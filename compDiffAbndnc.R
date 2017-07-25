#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
compDiffAbndnc.R (-h | --help | --version)
compDiffAbndnc.R DIR FEATURETABLE METADATA

Description:   This script calculates differential abundance using a feature table
Options:
--version       Show the current version.

Arguments:

DIR             Provide directory
FEATURETABLE    Provide feature table file
METADATE        Provide meta data file
" -> doc


args <- docopt(doc)

RESULTS_DIR <- args$DIR

cat("\nLoading arguments from", args$DIR, "\n")

load(paste(RESULTS_DIR, "cyttools.args.Rdata", sep = ""))

library(limma)
library(tidyverse)

logitTransform <- function(p) { log(p/(1-p)) }
panelDesign <- read.delim(paste(RESULTS_DIR, "nrsPanelFile.txt", sep = ""))
resultsTable <- read.delim(args$FEATURETABLE)
md <- read.delim(args$METADATA)

# ALL THIS SHOULD BE DONE BY THE CLUSTERING CALL OR THE CLUSTERING CALL SHOULD PASS ITS DATA ONTO ANOTHER SCRIPT TO DO SO
exprTable <- melt(resultsTable,
                  measure.vars = colnames(resultsTable)[colnames(resultsTable) %in% panelDesign$name[panelDesign$Ignore == 0] == T],
                  id.vars = colnames(resultsTable)[colnames(resultsTable) %in% panelDesign$name == F],
                  variable.name = "Metal",
                  value.name = "Intensity"
)
nodeExprTable <- acast(exprTable,
                       Mapping + Metal ~ FileNames,
                       fun.aggregate = median,
                       value.var = "Intensity" 
)

countTable <- resultsTable[,grep("Mapping|FileNames", colnames(resultsTable))]
countTable <- table(countTable$Mapping, countTable$FileNames)

colnames(countTable) <- md$FileName

props_table <- t(t(countTable) / colSums(countTable))

propData <- logitTransform(props_table[,c(1:18)])

targets <- md[c(1:18),]

targets$TimePoint <- factor(targets$TimePoint, levels = c("DBC_4HR", "DBC_24HR"))
targets$Condition <- factor(targets$Condition, levels = c("UNSTIM", "10NG", "25NG"))
targets$SampleID <- factor(targets$SampleID, levels = unique(targets$SampleID))

exprDesign <- paste(targets$Condition, targets$TimePoint, sep = "_")

design <- model.matrix(~ 0 + exprDesign + targets$SampleID)

colnames(design) <- gsub("exprDesign", "Cnd", colnames(design))
colnames(design) <- gsub("targets\\$SampleID", "", colnames(design))

fit <- lmFit(propData, design = design)
cont.matrix <- makeContrasts(Low_4HR = Cnd10NG_DBC_4HR - CndUNSTIM_DBC_4HR,
                             High_4HR = Cnd25NG_DBC_4HR - CndUNSTIM_DBC_4HR,
                             Low_24HR = Cnd10NG_DBC_24HR - CndUNSTIM_DBC_24HR,
                             High_24HR = Cnd25NG_DBC_24HR - CndUNSTIM_DBC_24HR,
                             levels = design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

diffAbndncStatsTable <- data.frame()

for ( i in 1:length(colnames(fit2))){
  nextPart <- topTable(fit2, coef = i, number = Inf)
  nextPart$Source <- rep(colnames(fit2)[i], nrow(nextPart))
  nextPart$Mapping <- row.names(nextPart)
  diffAbndncStatsTable <- rbind(diffAbndncStatsTable, nextPart)
}

fit <- lmFit(nodeExprTable[,c(1:18)], design = design)
cont.matrix <- makeContrasts(Low_4HR = Cnd10NG_DBC_4HR - CndUNSTIM_DBC_4HR,
                             High_4HR = Cnd25NG_DBC_4HR - CndUNSTIM_DBC_4HR,
                             Low_24HR = Cnd10NG_DBC_24HR - CndUNSTIM_DBC_24HR,
                             High_24HR = Cnd25NG_DBC_24HR - CndUNSTIM_DBC_24HR,
                             levels = design)

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


workspaceFile <- paste(RESULTS_DIR, "compDiffAbndncWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)
