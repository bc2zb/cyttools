#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
inspectr.R (-h | --help | --version)
inspectr.R DIR

Description:   This script perfroms QC analysis on spectral flow data, calulating similarity scores and spillover spreading
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

##########################################################################
############################ R code goes here ############################
##########################################################################


# capture all files in the directory
all_fcs_files <- list.files(args$DIR, full.names = T, pattern = "\\.fcs$")

# read in fcs files
ncfs <- read.ncdfFlowSet(all_fcs_files)

# transform flowSet
chnls <- colnames(ncfs)[grep("SSC|FSC|Time|\\-H", colnames(ncfs), invert = T)]
safe_estimate_logicle <- safely(estimateLogicle)
transFuncts <- fsApply(ncfs, safe_estimate_logicle, channels = paste0("^", chnls)) %>%
  modify_depth(1, 1) %>%
  discard(is_null)

safe_transform <- safely(transform)

for ( i in 1:length(transFuncts)){
  ncfs_trans <- safe_transform(ncfs, transFuncts[[i]])
  if(is.null(ncfs_trans$error)){
    ncfs_trans <- ncfs_trans$result
    break
  }else if(i == length(transFuncts)){
    cat("\nERROR: No transform can be estimated, exiting now\n")
    q()
  }
}

# preprocess flowFrames

if(args$clean == T){
  ncfs_trans <- fsApply(ncfs_trans, safe_preprocess_frame) %>%
    modify_depth(1, 1) %>%
    discard(is_null) %>%
    flowSet()
}

dir.create(paste0(RESULTS_DIR, "TRANSFORMED_FCS/"),
           showWarnings = F)

for( files in all_fcs_files){
  out.fcs.file <- paste0(RESULTS_DIR, "TRANSFORMED_FCS/logicle_transformed_", basename(files))
  out_fcs_frame <- ncfs_trans[basename(files)][[1]]
  write.FCS(out_fcs_frame, out.fcs.file)
}
##########################################################################
############################     End code     ############################
##########################################################################

