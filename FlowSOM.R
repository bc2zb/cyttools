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

colsToUse <- which(targets$name %in% lineage_markers[lineage_markers %in% targets$name[targets$Ignore == 0]] == T)

som <- BuildSOM(fsom,
                colsToUse = colsToUse,
                xdim = length(colsToUse),
                ydim = length(colsToUse))

som$prettyColnames <- targets$desc

flowSOM.res <- BuildMST(som)

ResultsTable <- as.data.frame(flowSOM.res$data)
fileLabels <- flowSOM.res$metaData %>% unlist() %>% matrix(ncol = 2, byrow = T)
fileNames <- vector(length = nrow(ResultsTable))
for(i in 1:nrow(fileLabels)){
  fileNames[fileLabels[i,1]:fileLabels[i,2]] <- file[i]
}

ResultsTable$FileNames <- fileNames
ResultsTable$Mapping <- flowSOM.res$map$mapping[,1]
ResultsTable$DistToNode <- flowSOM.res$map$mapping[,2]

ResultsTableFile <- paste(RESULTS_DIR, "FlowSOMResultsTable.txt", sep = "")

write.table(ResultsTable, ResultsTableFile, sep = "\t", quote = F, row.names = F)

panelDesign <- targets

nodeExprTable <- melt(ResultsTable,
                  measure.vars = colnames(ResultsTable)[colnames(ResultsTable) %in% panelDesign$name[panelDesign$Ignore == 0] == T],
                  id.vars = colnames(ResultsTable)[colnames(ResultsTable) %in% panelDesign$name == F],
                  variable.name = "Metal",
                  value.name = "Intensity"
) %>% acast(.,
            Mapping + Metal ~ FileNames,
            fun.aggregate = median,
            value.var = "Intensity"
)

countTable <- ResultsTable[,grep("Mapping|FileNames", colnames(ResultsTable))]
countTable <- table(countTable$Mapping, countTable$FileNames)

props_table <- t(t(countTable) / colSums(countTable))

nodeExprTableFile <- paste(RESULTS_DIR, "nodeExpressionFeatureTable.txt", sep = "")

write.table(nodeExprTable, nodeExprTableFile, sep = "\t", quote = F, row.names = T)

nodeAbndncFeatureTableFile <- paste(RESULTS_DIR, "nodeAbundanceFeatureTable.txt", sep = "")

write.table(props_table, nodeAbndncFeatureTableFile, sep = "\t", quote = F, row.names = T)

workspaceFile <- paste(RESULTS_DIR, "FlowSOMWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)
