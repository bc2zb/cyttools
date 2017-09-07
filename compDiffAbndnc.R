#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
compDiffAbndnc.R (-h | --help | --version)
compDiffAbndnc.R DIR

Description:   This script calculates differential abundance using a feature table
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

panelDesign <- read.delim(args$PANEL)
targets <- read.delim(args$METADATA)

colsToCheck <- c("TimePoint", "Condition", "SampleID", "FileName", "Group")
if(checkDesignCols(targets, colsToCheck)){
  missingCols <- colsToCheck[which(colsToCheck %in% colnames(targets) == F)]
  cat("\n\nERROR: PANEL file does not include required columns.
      \n\nMissing Columns:", missingCols,
      "\n\nPlease run cyttools.R --makeMetaDataBlank to generate compatible meta data file.\n\nStopping cyttools.R\n\n")
  q()
}

targets$TimePoint <- factor(targets$TimePoint, levels = unique(targets$TimePoint))
targets$Condition <- factor(targets$Condition, levels = unique(targets$Condition))
targets$SampleID <- factor(targets$SampleID, levels = unique(targets$SampleID))

if(length(unique(targets$TimePoint)) > 1){
  exprDesign <- paste(targets$Condition, targets$TimePoint, sep = "_")
}else{
  exprDesign <- targets$Condition
}

if(all(table(targets$SampleID) == 1)){
  if(length(unique(targets$Group)) <= 1){
      design <- model.matrix(~ 0 + exprDesign)
  }else{
    design <- model.matrix(~targets$Group + exprDesign)
  }
}else{
  design <- model.matrix(~targets$SampleID + exprDesign)
}
colnames(design) <- gsub("targets\\$SampleID|targets\\$Group", "BatchEffect", colnames(design))

orderList <- gsub("\\s", ".", targets$FileName)

props_table <- read.delim(args$FEATURETABLE, row.names = 1)
props_table <- props_table[,orderVectorByListOfTerms(colnames(props_table), orderList)]

propData <- logitTransform(props_table)

fit <- lmFit(propData, design = design)
if(checkLmFit(fit, propData)){
  design <- model.matrix(~ 0 + exprDesign + targets$SampleID)

  colnames(design) <- gsub("exprDesign", "Cnd", colnames(design))
  colnames(design) <- gsub("targets\\$SampleID|targets\\$Group", "BatchEffect", colnames(design))

  fit <- lmFit(propData, design = design)

  if(checkLmFit(fit, propData)){
    design <- model.matrix(~ 0 + exprDesign)

    colnames(design) <- gsub("exprDesign", "Cnd", colnames(design))
    colnames(design) <- gsub("targets\\$SampleID|targets\\$Group", "BatchEffect", colnames(design))

    fit <- lmFit(propData, design = design)

    cat("\n\nWARNING: Limma cannot estimate random effects due to sample size, using simplified model for analysis\n\n")
  }else{
    cat("\n\nWARNING: Limma cannot estimate group effects due to sample size, using simplified model for analysis\n\n")}
}
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

  if(all(unlist(lapply(c(firstTimePoint, secondTimePoint), is.na))) &
     firstCondition != secondCondition){
    validContrastIndex[i] <- T
    humanReadableColNames[i] <- paste(firstCondition,
                                      secondCondition,
                                      sep = "_vs_")
  }else if(any(unlist(lapply(c(firstCondition, secondCondition, firstTimePoint, secondTimePoint), is.na)))){
    validContrastIndex[i] <- F
  }else if(firstCondition == secondCondition & firstTimePoint != secondTimePoint){
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

cont.matrix.2 <- cont.matrix[,validContrastIndex] %>% as.matrix()
colnames(cont.matrix) <- humanReadableColNames[validContrastIndex]

if(any(grepl("BatchEffect", humanReadableColNames))){
  cont.matrix <- cont.matrix[,grep("BatchEffect", colnames(cont.matrix), invert = T)]
}

if(is.null(dim(cont.matrix))){
  fit2 <- eBayes(fit)
}else{
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
}

diffAbndncStatsTable <- data.frame()

for ( i in 1:length(colnames(fit2))){
  nextPart <- topTable(fit2, coef = i, number = Inf)
  nextPart$Source <- rep(colnames(fit2)[i], nrow(nextPart))
  nextPart$Mapping <- row.names(nextPart)
  diffAbndncStatsTable <- rbind(diffAbndncStatsTable, nextPart)
}

nodeAbndncStatsFile <- paste(RESULTS_DIR, "nodeDifferentialAbundanceTable.txt", sep = "")

write.table(diffAbndncStatsTable, nodeAbndncStatsFile, sep = "\t", quote = F, row.names = F)

workspaceFile <- paste(RESULTS_DIR, "compDiffAbndncWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)
