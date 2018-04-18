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

DIR    Provide directory for cytools.args.Rdata to be found

" -> doc


args <- docopt(doc)

ARGS_DIR <- args$DIR

cat("\nLoading arguments from", ARGS_DIR, "\n")

load(paste(ARGS_DIR, "cyttools.args.Rdata", sep = ""))

RESULTS_DIR <- args$OUT

source("cyttoolsFunctions.R")

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
  flowSet.trans <- read.flowSet(file, transformation = F, truncate_max_range = F)
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

panelDesign <- targets

medExprData <- flowSOM.res$map$medianValues[,panelDesign$Ignore == 0]
medExprData <- t(medExprData)
colnames(medExprData) <- c(1:ncol(medExprData))

nodeExprTable <- ResultsTable %>%
  group_by_at(vars(FileNames, Mapping)) %>%
  summarise_at(colnames(ResultsTable)[colnames(ResultsTable) %in% targets$name[targets$Functional == 1 | targets$Lineage == 1]],
               median) %>%
  gather(Metal,
         Intensity,
         -Mapping, 
         -FileNames) %>%
  spread(FileNames,
         Intensity)

countTable <- ResultsTable[,grep("Mapping|FileNames", colnames(ResultsTable))]
countTable <- table(countTable$Mapping, countTable$FileNames)
props_table <- t(t(countTable) / colSums(countTable))

ResultsTable <- ResultsTable %>% left_join(flowSOM.res$MST$l %>%
                                             as.data.frame() %>%
                                             setNames(c("cyttools_dim_x", "cyttools_dim_y")) %>%
                                             rownames_to_column("Mapping") %>%
                                             mutate(Mapping = as.numeric(Mapping)))
dir.create(paste0(RESULTS_DIR, "CLUSTERED_FCS/"),
           showWarnings = F)
for( files in file){
  rawFCS <- read.FCS(files, transformation = F)
  clusterData <- ResultsTable %>%
    dplyr::filter(FileNames == files) %>%
    select(Mapping, DistToNode, cyttools_dim_x, cyttools_dim_y)
  clusterFCS <- flowCore::cbind2(rawFCS, as.matrix(clusterData))
  row.names(pData(parameters(clusterFCS))) <- paste0("$P", c(1:nrow(pData(parameters(clusterFCS)))))
  out.fcs.file <- paste0(RESULTS_DIR, "CLUSTERED_FCS/clustered_", basename(files))
  write.FCS(clusterFCS, out.fcs.file)
}

# # write out results
# ResultsTableFile <- paste(RESULTS_DIR, "FlowSOMResultsTable.txt", sep = "")
# nodeExprTableFile <- paste(RESULTS_DIR, "nodeExpressionFeatureTable.txt", sep = "")
# nodeAbndncFeatureTableFile <- paste(RESULTS_DIR, "nodeAbundanceFeatureTable.txt", sep = "")
# nodeCountFeatureTableFile <- paste(RESULTS_DIR, "nodeCountFeatureTable.txt", sep = "")
# nodeMedianFeatureTableFile <- paste(RESULTS_DIR, "nodeMedianFeatureTable.txt", sep = "")
# 
# write.table(ResultsTable, ResultsTableFile, sep = "\t", quote = F, row.names = F)
# write.table(nodeExprTable, nodeExprTableFile, sep = "\t", quote = F, row.names = T)
# write.table(props_table, nodeAbndncFeatureTableFile, sep = "\t", quote = F, row.names = T)
# write.table(countTable, nodeCountFeatureTableFile, sep = "\t", quote = F, row.names = T)
# write.table(flowSOM.res$map$medianValues, nodeMedianFeatureTableFile, sep = "\t", quote = F, row.names = T)

# workspaceFile <- paste(RESULTS_DIR, "FlowSOMWorkspace.Rdata", sep = "")
# 
# save.image(file = workspaceFile)
