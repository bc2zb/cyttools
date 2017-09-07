#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
compDiffCount.R (-h | --help | --version)
compDiffCount.R DIR

Description:   This script calculates differential counts using a feature table
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

exprDesign <- targets$Condition

orderList <- gsub("\\s", ".", targets$FileName)

props_table <- read.delim(args$FEATURETABLE, row.names = 1)
props_table <- props_table[,orderVectorByListOfTerms(colnames(props_table), orderList)]

propData <- DGEList(props_table[(props_table %>% apply(1, max)) > 100,],
                    lib.size = colSums(props_table))


diffAbndncStatsTable <- tibble()

for (baseline in levels(exprDesign)){
  
  tmpExprDesign <- relevel(exprDesign, baseline)
  design <- model.matrix(~targets$Group + tmpExprDesign)
  colnames(design) <- gsub("tmpExprDesign", "Cnd.", colnames(design))
  colnames(design) <- gsub("targets\\$SampleID|targets\\$Group", "BatchEffect", colnames(design))

  fit <- estimateDisp(propData, design)
  fit <- glmQLFit(fit, design, robust=TRUE)
  
  for(experimental in (levels(tmpExprDesign)[-1])){
    res <- glmQLFTest(fit, coef = paste0("Cnd.", experimental))
    topTable <- topTags(res, n = Inf)$table %>%
      rownames_to_column("ClusterID") %>%
      mutate(Condition = rep(experimental, nrow(topTable)),
                                    Baseline = rep(baseline, nrow(topTable)))
    diffAbndncStatsTable <- diffAbndncStatsTable %>% bind_rows(topTable)
                
  }

}
nodeAbndncStatsFile <- paste(RESULTS_DIR, "nodeDifferentialCountTable.txt", sep = "")

write.table(diffAbndncStatsTable, nodeAbndncStatsFile, sep = "\t", quote = F, row.names = F)

workspaceFile <- paste(RESULTS_DIR, "compDiffCountWorkspace.Rdata", sep = "")

save.image(file = workspaceFile)
