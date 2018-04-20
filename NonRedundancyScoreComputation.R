#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
NonRedundancyScoreComputation.R (-h | --help | --version)
NonRedundancyScoreComputation.R DIR

Description:   This script is a template for making docopts compatible Rscripts
Options:
--version       Show the current version.

Arguments:

DIR    Provide directory where cyttools.args.Rdata file is located
" -> doc


args <- docopt(doc)

ARGS_DIR <- args$DIR

cat("\nLoading arguments from", ARGS_DIR, "\n")

load(paste(ARGS_DIR, "cyttools.args.Rdata", sep = ""))

source("cyttoolsFunctions.R")

dir <- args$DIR # grabs directory from initial cyttools call
file <- list.files(dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory
targets <- read.delim(args$PANEL)
colsToCheck <- c("Ignore", "TransformCofactor", "Lineage", "Functional")
if(checkDesignCols(targets, colsToCheck)){
  cat("\n\nERROR: PANEL file does not include required columns.
      \n\nMissing Columns:", colsToCheck[which(colsToCheck %in% colnames(targets) == F)],
      "\n\nPlease run cyttools.R --makePanelBlank to generate compatible panel file.\n\nStopping cyttools.R\n\n")
  q()
}

lineage_markers <- targets$name[targets$Lineage == 1]
functional_markers <- targets$name[targets$Functional == 1]

if(args$transform == "logicle"){
  # read in fcs files
  ncfs <- read.ncdfFlowSet(file)
  
  chnls <- colnames(ncfs)[grep("SSC|FSC|Time|\\-H", colnames(ncfs), invert = T)]
  safe_estimate_logicle <- safely(estimateLogicle)
  transFuncts <- fsApply(ncfs, safe_estimate_logicle, channels = paste0("^", chnls)) %>%
    modify_depth(1, 1) %>%
    discard(is_null)
  
  safe_transform <- safely(transform)
  
  for ( i in 1:length(transFuncts)){
    ncfs_trans <- safe_transform(ncfs, transFuncts[[i]])
    if(is.null(ncfs_trans$error)){
      flowSet.trans <- as.flowSet(ncfs_trans$result)
      break
    }else if(i == length(transFuncts)){
      cat("\nERROR: No transform can be estimated, exiting now\n")
      q()
    }
  }
}else if(args$transform == "arcsinh"){
  flowSet.trans <- read.flowSet.transVS(targets, file)
}else if(args$transform == "none"){
  flowSet.trans <- read.flowSet(file, transformation = F, truncate_max_range = F)
}else{
  cat("\nNo transform specified, exiting now\n")
  q()
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

RESULTS_DIR <- args$OUT
nrsFile <- paste(RESULTS_DIR, "nrsPanelFile.txt", sep = "")

write.table(targets, file = nrsFile, sep = "\t", quote = F, row.names = F)

# workspaceFile <- paste(RESULTS_DIR, "NonRedundancyScoreComputationWorkspace.Rdata", sep = "")
# 
# save.image(file = workspaceFile)