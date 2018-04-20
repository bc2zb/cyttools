#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
phenotype_gridpoints.R (-h | --help | --version)
phenotype_gridpoints.R DIR

Description:   This script phenotypes grid points from FlowSOM with FlowType and performs consensus clustering
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

# read in panel design
targets <- read.delim(args$PANEL)
colsToCheck <- c("Ignore", "TransformCofactor", "Lineage", "Functional", "NRS")
if(checkDesignCols(targets, colsToCheck)){
  missingCols <- colsToCheck[which(colsToCheck %in% colnames(targets) == F)]
  cat("\n\nERROR: PANEL file does not include required columns.
      \n\nMissing Columns:", missingCols,
      "\n\nPlease run cyttools.R --makePanelBlank and cyttools.R --computeNRS to generate compatible panel file.\n\nStopping cyttools.R\n\n")
  q()
}

# read in clusterd FCS files
cluster_dir <- args$CLUSTERDIR # grabs directory from initial cyttools call
file <- list.files(cluster_dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory

lineage_markers <- targets$name[targets$Lineage == 1]
functional_markers <- targets$name[targets$Functional == 1]

cluster.flowSet.trans <- read.flowSet(file, transformation = F, truncate_max_range = F)

# read in phenotyped FCS files
pheno_dir <- args$PHENODIR # grabs directory from initial cyttools call

pheno_file <- list.files(pheno_dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory
pheno.flowSet.trans <- read.flowSet(pheno_file, transformation = F, truncate_max_range = F)

mappings <- cluster.flowSet.trans %>%
  fsApply(function(x){return(as.data.frame(exprs(x)))}, simplify = F) %>%
  bind_rows(.id = "FileNames") %>%
  select(Mapping) %>%
  bind_cols(pheno.flowSet.trans %>%
              fsApply(function(x){return(as.data.frame(exprs(x)))}, simplify = F) %>%
              bind_rows(.id = "FileNames") %>%
              select(contains("Phenotype_")))

map_counts <- mappings %>% 
  gather(PhenoCodes,
         Phenotype,
         -Mapping) %>% 
  group_by(Mapping, PhenoCodes, factor(Phenotype)) %>%
  tally()

diff_counts <- map_counts %>%
  spread(`factor(Phenotype)`,
         `n`,
         fill = NA)

if(str_detect(colnames(diff_counts)[length(colnames(diff_counts))], "3")){
  diff_counts <- diff_counts
}else{
  diff_counts <- diff_counts %>%
    ungroup() %>%
    mutate(`3` = rep(as.double(NA), nrow(.)))
}

reduced_phenocodes <- diff_counts %>%
  mutate(MarkerCode = case_when(is.na(`3`) & (`2`/(`1` + `2`)) > 0.75 ~ 2,
                                is.na(`3`) & (`1`/(`1` + `2`)) > 0.75 ~ 1,
                                is.na(`3`) == F & (`3`/(`3` + `2` + `1`)) > 0.75 ~ 3,
                                is.na(`3`) == F & (`2`/(`3` + `2` + `1`)) > 0.75 ~ 2,
                                is.na(`3`) == F & (`1`/(`3` + `2` + `1`)) > 0.75 ~ 1,
                                TRUE ~ 0)) %>%
  select(Mapping, PhenoCodes, MarkerCode) %>%
  #mutate(MarkerCode = factor(MarkerCode)) %>%
  spread(PhenoCodes, MarkerCode) %>%
  select(starts_with("Phenotype_"),
         Mapping)

phenocode_matrix <- lapply(seq_along(colnames(reduced_phenocodes %>%
                                                ungroup() %>%
                                                select(-Mapping))),
                   function(x){
                     marker <- gsub("Phenotype_", "", colnames(reduced_phenocodes)[x])
                     phenocode_vector <- reduced_phenocodes[,x]
                     if(max(as.numeric(as.character(phenocode_vector[[1]]))) == 3){
                       phenocode_vector[phenocode_vector == 3] <- paste0(marker, "hi")
                       phenocode_vector[phenocode_vector == 2] <- paste0(marker, "lo")
                       phenocode_vector[phenocode_vector == 1] <- paste0(marker, "-")
                       phenocode_vector[phenocode_vector == 0] <- ""
                     }else{
                       phenocode_vector[phenocode_vector == 2] <- paste0(marker, "+")
                       phenocode_vector[phenocode_vector == 1] <- paste0(marker, "-")
                       phenocode_vector[phenocode_vector == 0] <- ""
                     }
                     return(phenocode_vector)
                             }) %>%
  bind_cols()

lineage_lists <- phenocode_matrix %>%
  apply(1, function(x){
    reduced_x <- x[grep("^$", x, invert = T)]
    return(lapply(seq_along(1:length(reduced_x)), function(y){
      combn(reduced_x, y, simplify = T) %>%
        apply(2, paste, collapse = " ") %>%
        return()
      }))
    })


annotated_mappings <- data_frame(Mapping = 1:length(lineage_lists),
                                 Immunophenotypes = lineage_lists)

annotated_mappings_file <- paste(RESULTS_DIR, "annotated_mappings.txt", sep = "")
distinct_mappings_file <- paste(RESULTS_DIR, "distinct_mappings.txt", sep = "")


annotated_mappings %>%
  unnest() %>% 
  unnest() %>% 
  group_by(Immunophenotypes) %>%
  summarise(Mappings = paste0(Mapping, collapse = ",")) %>%
  mutate(Num_Markers = str_count(Immunophenotypes, "\\+|\\-|lo|hi"),
         Num_GridPoints = str_count(Mappings, ",") + 1) %>%
  arrange(desc(Num_Markers)) %>%
  distinct(Mappings,
           .keep_all = T) %>% 
  write_tsv(distinct_mappings_file)

annotated_mappings %>%
  unnest() %>% 
  unnest() %>% 
  write_tsv(annotated_mappings_file)
  
  
# immunophenotypes_order <- order(str_count(row.names(filtered_count_table)), decreasing = T)
# 
# unique_filtered_count_table <- unique(filtered_count_table[immunophenotypes_order,])
# uniq_immunophenotypes <- row.names(filtered_count_table)[immunophenotypes_order][!duplicated(filtered_count_table[immunophenotypes_order,])]

  
##########################################################################
############################     End code     ############################
##########################################################################

# workspaceFile <- paste(RESULTS_DIR, "phenotype_gridpoints.Workspace.Rdata", sep = "")

# save.image(file = workspaceFile)
