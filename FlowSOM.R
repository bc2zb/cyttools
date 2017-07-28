#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
FLowSOM.R (-h | --help | --version)
FLowSOM.R DIR

Description:   This script is a template for making docopts compatible Rscripts
Options:
--version       Show the current version.

Arguments:

DIR    Provide directory

" -> doc


args <- docopt(doc)

ARGS_DIR <- args$DIR

cat("\nLoading arguments from", ARGS_DIR, "\n")

load(paste(ARGS_DIR, "cyttools.args.Rdata", sep = ""))

RESULTS_DIR <- args$OUT

source("cyttoolsFunctions.R")

file <- list.files(dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory

targets <- read.delim(args$PANEL)
colsToCheck <- c("Ignore", "TransformCofactor", "Lineage", "Functional", "NRS")
if(checkDesignCols(targets, colsToCheck)){
  missingCols <- colsToCheck[which(colsToCheck %in% colnames(targets) == F)]
  cat("\n\nERROR: PANEL file does not include required columns.
      \n\nMissing Columns:", missingCols,
      "\n\nPlease run cyttools.R --makePanelBlank and cyttools.R --computeNRS to generate compatible panel file.\n\nStopping cyttools.R\n\n")
  q()
}

dir <- args$DIR # grabs directory from initial cyttools call
file <- list.files(dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory

lineage_markers <- targets$name[targets$Lineage == 1]
functional_markers <- targets$name[targets$Functional == 1]

if(args$transform == T){
  flowSet.trans <- read.flowSet.transVS(targets, file)
}else{
  flowSet.trans <- read.flowSet(file)
}

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
