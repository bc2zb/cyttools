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
all_fcs_files <- list.files(args$DIR, full.names = T, pattern = "\\.fcs$")

# find files labelled as "Reference Group"
ref_group_labelled_fcs_files <- all_fcs_files[str_detect(all_fcs_files, "Reference Group|Super Bright")]

# read in reference group fcs files
ncfs <- read.ncdfFlowSet(ref_group_labelled_fcs_files)

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
preprocessed_set <- fsApply(ncfs_trans, safe_preprocess_frame) %>%
  modify_depth(1, 1) %>%
  discard(is_null) %>%
  flowSet()

# are there a substantial number of raw detectors in the flowSet
if(length(which(colnames(preprocessed_set) %in% detector_order == T)) == length(detector_order)){
  # if yes, calculate similarity scores
  # find channel with MAX intensity
  median_values <- fsApply(preprocessed_set, find_peak_emission_detector)
  # create tidy data frame to perform correlation and euclidean distance calculations
  cor_data <- median_values %>%
    as.data.frame() %>%
    select(one_of(detector_order)) %>%
    rownames_to_column("file_names") %>%
    mutate(sample_id = gsub("Reference Group_|\\_[0-9]*\\_[0-9]*\\.fcs$|\\.fcs$", "", file_names),
           experiment_id = rep(args$DIR, nrow(median_values)))
  
  write_csv(cor_data,
            path = paste0(RESULTS_DIR,
                          "median-value-data-experiment-id-",
                          make.names(cor_data$experiment_id[1]),
                          ".csv"))
  
  # make matrix to perform correlation coefficient computation
  cor_matrix <- cor_data %>%
    select(one_of(detector_order)) %>%
    as.matrix()
  
  # make data frame containing label information for each fluorophore
  cor_matrix_row_info <- cor_data %>%
    select(-one_of(detector_order))
  
  # calculate the correlation coefficient between all fluorophores
  cor_coefficient_matrix <- cor(t(cor_matrix)) %>%
    as.data.frame() %>%
    rownames_to_column("row_id") %>%
    gather(col_id,
           correlation_coefficient,
           -row_id)
  
  # calculate euclidean distance between all fluorophores
  euc_distance_matrix <- dist(cor_matrix) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("row_id") %>%
    gather(col_id,
           euclidean_distance,
           -row_id) %>%
    mutate(col_id = paste0("V", col_id))
  
  # create tidy dataframe for correlation coefficient and euclidean distance between all pairs of fluorophores
  similarity_data <- cor_coefficient_matrix %>%
    left_join(euc_distance_matrix) %>%
    left_join(cor_matrix_row_info %>%
                setNames(paste0("col_", colnames(.))) %>%
                rownames_to_column("row_id") %>%
                mutate(col_id = paste0("V", row_id),
                       row_id = NULL)) %>%
    left_join(cor_matrix_row_info %>%
                setNames(paste0("row_", colnames(.))) %>%
                rownames_to_column("row_id"))
    write_csv(similarity_data,
              path = paste0(RESULTS_DIR,
                     "similarity-data-experiment-id-",
                     make.names(similarity_data$col_experiment_id[1]),
                     ".csv"))
}else{
  # if no, calculate spillover spreading
  # separate out unstained and stained frames
  if(any(str_detect(sampleNames(preprocessed_set), regex("Unstained", ignore_case = T)))){
    unstained_set <- preprocessed_set[str_detect(sampleNames(preprocessed_set), regex("Unstained", ignore_case = T))]
  }
  
  stained_set <- preprocessed_set[!str_detect(sampleNames(preprocessed_set), regex("Unstained", ignore_case = T))]
  # iterate through stained flowFrames, calculating spillover spreading
  spreading_statistics <- fsApply(stained_set, function(flow_frame){
    fluorophore <- str_remove_all(chnls, "\\-A$") %>%
      paste0("[\\ |\\_]", ., " \\(") %>%
      str_detect(flow_frame@description$FILENAME, .) %>%
      chnls[.]
    cat(fluorophore, "\n\n")
    positive_gate <- gate_mindensity(flow_frame, fluorophore)
    gated_set <- split(flow_frame, flowCore::filter(flow_frame, positive_gate))
    stained_frame <- gated_set$`+`
    
    if(exists("unstained_set")){
      refernc_pops_set <- rbind2(unstained_set, gated_set$`-`)
    }else{
      refernc_pops_set <- gated_set$`-`
    }
    
    spreading_statistics <- fsApply(refernc_pops_set, function(refernc_pop){
      res <- exprs(refernc_pop)[,chnls[!chnls %in% fluorophore]] %>%
        data.frame() %>%
        gather(Fluorophore,
               Intensity) %>%
        group_by(Fluorophore) %>%
        summarise(ref_width = (quantile(Intensity, .84) - quantile(Intensity, .5))) %>%
        left_join(exprs(stained_frame)[,chnls[!chnls %in% fluorophore]] %>%
                    data.frame() %>%
                    gather(Fluorophore,
                           Intensity) %>%
                    group_by(Fluorophore) %>%
                    summarise(stn_width = (quantile(Intensity, .84) - quantile(Intensity, .5)))) %>%
        mutate(spillover_spreading_value = sqrt((stn_width^2) - (ref_width^2)),
               stat = spillover_spreading_value/(sqrt(((median(exprs(stained_frame)[,fluorophore])) - (median(exprs(refernc_pop)[,fluorophore]))))),
               Primary_Detector = rep(fluorophore, nrow(.)),
               Primary_Detector_FileNames = rep(stained_frame@description$FILENAME, nrow(.)),
               Spillover_Channel_FileNames = rep(refernc_pop@description$FILENAME, nrow(.)))
        return(res)
    
  })
    return(spreading_statistics %>% bind_rows())
  })
  spreading_statistics %>%
    bind_rows() %>%
    mutate(Experiment = rep(args$DIR, nrow(.))) %>%
    write_csv(path = paste0(RESULTS_DIR,
                            "spillover-data-experiment-id-",
                            make.names(args$DIR),
                            ".csv"))
}

##########################################################################
############################     End code     ############################
##########################################################################

