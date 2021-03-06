#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
ParFlowType.R (-h | --help | --version)
ParFlowType.R DIR

Description:   This script applies the flowType algorthim in parallel using mclapply() to high parameter cytometry data
Options:
--version       Show the current version.

Arguments:

DIR    Provide directory for cytools.args.Rdata to be found, automatically generated by invoking cyttools.R
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


lineage_markers <- targets$name[targets$Lineage == 1]
functional_markers <- targets$name[targets$Functional == 1]

dir <- args$DIR # grabs directory from initial cyttools call
file <- list.files(dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory

if(args$transform == T){
  flowSet.trans <- read.flowSet.transVS(targets, file)
}else{
  flowSet.trans <- read.flowSet(file)
}

# order the markers using NRS, dropping markers set to "1" in Ignore column of panel design
lineage_markers_ord <- targets$name
lineage_markers_ord <- lineage_markers_ord[lineage_markers_ord %in% targets$name[which(targets$Ignore == 0)]]
lineage_markers_ord <- lineage_markers_ord[order(targets$NRS[which(targets$name %in% lineage_markers_ord)], decreasing = T)]
if(length(lineage_markers_ord) > 12){
  lineage_markers_ord <- lineage_markers_ord[1:12]
}

colsToUse <- which(targets$name %in% lineage_markers_ord == T)

ResList <- mclapply(fs.as.list(flowSet.trans),
                    flowType,
                    PropMarkers = colsToUse,
                    MFIMarkers = colsToUse,
                    MarkerNames = targets$desc[colsToUse],
                    MemLimit = 16,
                    mc.cores = 8)

panelDesign <- targets

phenotype.names=unlist(lapply(ResList[[1]]@PhenoCodes,function(x){
  return(decodePhenotype(x,
                         as.character(targets$desc[colsToUse]),
                         ResList[[1]]@PartitionsPerMarker))}))
names(ResList[[1]]@PhenoCodes)=phenotype.names

nodeExprTable <- lapply(ResList, function(x){return(apply(x@Partitions, 1, paste, collapse =""))})

allExprData <- fsApply(flowSet.trans, function(x){return(x)}, use.exprs = T)

nodeMappings <- unlist(nodeExprTable) %>% data.frame(CellId = names(.), Mapping = .)

nodeExprTable <- cbind(allExprData, nodeMappings)

nodeExprTable$NodeNames <- recoderFunc(nodeExprTable$Mapping,
                                       ResList[[1]]@PhenoCodes,
                                       names(ResList[[1]]@PhenoCodes))

nodeExprTable$FileNames <- gsub(".fcs[0-9]*", ".fcs", nodeExprTable$CellId)

ResultsTableFile <- paste(RESULTS_DIR, "FlowTypeResultsTable.txt", sep = "")

write.table(nodeExprTable, ResultsTableFile, sep = "\t", quote = F, row.names = F)

subPopsExprTable <- matrix(ncol = length(ResList))
for ( i in 1:length(ResList[[1]]@PhenoCodes)){
  subPopCode <- gsub("0", ".", ResList[[1]]@PhenoCodes[i])
  subPopData <- nodeExprTable[grep(subPopCode, nodeExprTable$Mapping),]
  subPopExprTable <- melt(subPopData,
                          measure.vars = colnames(subPopData)[colnames(subPopData) %in% panelDesign$name[panelDesign$Ignore == 0] == T],
                          id.vars = colnames(subPopData)[colnames(subPopData) %in% panelDesign$name == F],
                          variable.name = "Metal",
                          value.name = "Intensity"
  )
  subPopExprTable$Mapping <- rep(gsub("\\.", "0", subPopCode), nrow(subPopExprTable))
  subPopNodeExprTable <- acast(subPopExprTable,
                               Mapping + Metal ~ FileNames,
                               fun.aggregate = median,
                               value.var = "Intensity" 
  )
  subPopsExprTable <- rbind(subPopsExprTable, subPopNodeExprTable)
}

nodeExprTableFile <- paste(RESULTS_DIR, "nodeExpressionFeatureTable.txt", sep = "")

write.table(subPopsExprTable, nodeExprTableFile, sep = "\t", quote = F, row.names = T)

all.proportions <- matrix(0,length(ResList[[1]]@CellFreqs),length(ResList))
for (i in 1:length(ResList))
  all.proportions[,i] = ResList[[i]]@CellFreqs / ResList[[i]]@CellFreqs[1]

colnames(all.proportions) <- names(ResList)
row.names(all.proportions) <- phenotype.names

nodeAbndncFeatureTableFile <- paste(RESULTS_DIR, "nodeAbundanceFeatureTable.txt", sep = "")

write.table(all.proportions, nodeAbndncFeatureTableFile, sep = "\t", quote = F, row.names = T)

PhenoCodes <- data.frame(PhenoCodes = ResList[[1]]@PhenoCodes,
                         Names = names(ResList[[1]]@PhenoCodes))

PhenoCodesFile <- paste(RESULTS_DIR, "PhenoCodes.txt", sep = "")

write.table(PhenoCodes, nodeAbndncFeatureTableFile, sep = "\t", quote = F, row.names = F)

workspaceFile <- paste(RESULTS_DIR, "FlowTypeWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)