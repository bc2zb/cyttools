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

# order detectors
detector_order <- c("V1-A",  "V2-A",  "V3-A",  "V4-A",  "V5-A",  "V6-A",  "V7-A",  "V8-A",  "V9-A",  "V10-A", "V11-A", "V12-A", "V13-A", "V14-A", "V15-A", "V16-A",
                    "B1-A", "B2-A", "B3-A", "B4-A", "B5-A", "B6-A", "B7-A", "B8-A", "B9-A", "B10-A", "B11-A", "B12-A", "B13-A", "B14-A",
                    "R1-A", "R2-A", "R3-A", "R4-A", "R5-A", "R6-A", "R7-A", "R8-A")

# capture all files in the directory
all_fcs_files <- list.files(args$DIR, full.names = T)

# find files labelled as "Reference Group"
ref_group_labelled_fcs_files <- all_fcs_files[str_detect(all_fcs_files, "Reference Group|Super Bright")]

# read in reference group fcs files
ncfs <- read.ncdfFlowSet(ref_group_labelled_fcs_files)

# transform flowSet
chnls <- colnames(ncfs)[grep("SSC|FSC|Time|\\-H", colnames(ncfs), invert = T)]
transFuncts <- estimateLogicle(ncfs[[1]], channels = paste0("^", chnls))
ncfs_trans <- transform(ncfs, transFuncts)

# preprocess flowFrames
preprocess_frame <- function(flow_frame){
  # gate out debris
  nondebris_gate <- gate_mindensity(flow_frame, "FSC-A")
  # split flow_frame into positive and negative frames
  nondebris_set <- split(flow_frame,
                         flowCore::filter(flow_frame, nondebris_gate))
  # beads and cells have different responses to this gate, grab whichever has the higher number of events
  num_pos_events <- nrow(nondebris_set$`+`)
  num_neg_events <- nrow(nondebris_set$`-`)
  if(num_pos_events > num_neg_events){
    nondebris_frame <- nondebris_set$`+`
  }else{
    nondebris_frame <- nondebris_set$`-`
  }
  # gate singlets
  singlets_gate <- gate_singlet(nondebris_frame, maxit = 20, wider_gate = T)
  singlets_set <- split(nondebris_frame,
                        flowCore::filter(nondebris_frame, singlets_gate))
  # retrive singlets only
  singlets_frame <- singlets_set$`singlet+`
  # perform ellipsoid gating to tighten up population
  ellipsoid_gate <- gate_flowClust_2d(singlets_frame,
                                      xChannel = "FSC-A",
                                      yChannel = "SSC-A")
  ellipsoid_set <- split(singlets_frame,
                         flowCore::filter(singlets_frame, ellipsoid_gate))
  # retrieve events from within ellipsoid
  preprocessed_frame <- ellipsoid_set$`defaultEllipsoidGate+`
  return(preprocessed_frame)
}
safe_preprocess_frame <- safely(preprocess_frame)

preprocessed_set <- fsApply(ncfs_trans, safe_preprocess_frame)

preprocessed_set %>% modify_depth(1, 1) %>% 
##########################################################################
############################     End code     ############################
##########################################################################

