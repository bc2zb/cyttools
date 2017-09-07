#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
compDiffExpr.R (-h | --help | --version)
compDiffExpr.R DIR

Description:   This script calculates differential abundance using a feature table
Options:
--version       Show the current version.

Arguments:

DIR             Provide directory where cyttools.args.Rdata file is located
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

nodeExprTable <- read.delim(args$FEATURETABLE, row.names = 1)
nodeExprTable <- nodeExprTable[,orderVectorByListOfTerms(colnames(nodeExprTable), c(colnames(nodeExprTable)[1:2],
                                                                                    orderList))]

nodeExprTable <- nodeExprTable %>%
  mutate(RowID = paste(ConsensusCluster, Metal, sep = "_")) %>%
  select(-ConsensusCluster,
         -Metal) %>%
  column_to_rownames("RowID")

diffExprStatsTable <- tibble()

for (baseline in levels(exprDesign)){
  
  tmpExprDesign <- relevel(exprDesign, baseline)
  design <- model.matrix(~targets$Group + tmpExprDesign)
  colnames(design) <- gsub("tmpExprDesign", "Cnd.", colnames(design))
  colnames(design) <- gsub("targets\\$SampleID|targets\\$Group", "BatchEffect", colnames(design))
  
  fit <- lmFit(nodeExprTable, design = design)
  res <- eBayes(fit)
  
  for(experimental in (levels(tmpExprDesign)[-1])){
    topTable <- topTable(res, coef = paste0("Cnd.", experimental), number = Inf) %>%
      rownames_to_column("ClusterID") %>% 
      separate(ClusterID, c("ClusterID", "Observation"), "_") %>%
      mutate(Condition = rep(experimental, nrow(.)),
             Baseline = rep(baseline, nrow(.)))
    diffExprStatsTable <- diffExprStatsTable %>% bind_rows(topTable)
    
  }
  
}

diffExprStatsTable <- diffExprStatsTable %>%
  select(ClusterID, Observation, P.Value, adj.P.Val, Condition, Baseline, logFC) %>%
  mutate(PValue = P.Value,
         FDR = adj.P.Val,
         adj.P.Val = NULL,
         P.Value = NULL)

nodeExprStatsFile <- paste(RESULTS_DIR, "nodeDifferentialExpressionTable.txt", sep = "")

write.table(diffExprStatsTable, nodeExprStatsFile, sep = "\t", quote = F, row.names = F)

# workspaceFile <- paste(RESULTS_DIR, "compDiffExprWorkspace.Rdata", sep = "")
# 
# save.image(file = workspaceFile)
