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

if(args$transform == T){
  cluster.flowSet.trans <- read.flowSet.transVS(targets, file)
}else{
  cluster.flowSet.trans <- read.flowSet(file)
}

# read in phenotyped FCS files
pheno_dir <- args$PHENODIR # grabs directory from initial cyttools call

pheno_file <- list.files(pheno_dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory

if(args$transform == T){
  pheno.flowSet.trans <- read.flowSet.transVS(targets, pheno_file)
}else{
  pheno.flowSet.trans <- read.flowSet(pheno_file)
}

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
         fill = NA) %>%
  mutate(`1` = if_else(is.na(`1`),
                       as.double(0),
                       as.double(`1`)),
         `2` = if_else(is.na(`2`),
                       as.double(0),
                       as.double(`2`)),
          diff = if_else(is.na(`3`),
                        (`2` - `1`)/(`2` + `1`),
                        (`3` - (`2` + `1`))/(`3` + `2` + `1`)))

reduced_phenocodes <- diff_counts %>%
  mutate(MarkerCode = case_when(diff < 0.75 & diff > -0.75 ~ 0,
                                is.na(`3`) & diff >= 0.75 ~ 2,
                                is.na(`3`) & diff <= -0.75 ~ 1,
                                is.na(`3`) == F & diff >= 0.75 ~ 3,
                                is.na(`3`) == F & diff <= -0.75 & ((`2` - (`3` + `1`))/(`3` + `2` + `1`)) >= -0.75 ~ 2,
                                is.na(`3`) == F & diff <= -0.75 & ((`1` - (`3` + `2`))/(`3` + `2` + `1`)) >= -0.75 ~ 1)) %>%
  select(Mapping, PhenoCodes, MarkerCode) %>%
  mutate(MarkerCode = factor(MarkerCode)) %>%
  spread(PhenoCodes, MarkerCode) %>%
  select(starts_with("Phenotype_"),
         Mapping)

medClusterCount <- cluster.flowSet.trans %>%
  fsApply(function(x){return(as.data.frame(exprs(x)))}, simplify = F) %>%
  bind_rows(.id = "FileNames") %>%
  group_by(Mapping, FileNames) %>%
  summarise(n()) %>%
  mutate(rare_pop_score = `n()`/100) %>%
  ungroup() %>%
  group_by(Mapping) %>%
  summarise(max_rare_pop = max(rare_pop_score),
            med_rare_pop = median(rare_pop_score),
            min_rare_pop = min(rare_pop_score))

### MERGING FUNCTION GOES HERE?

phenocode_matrix <- lapply(seq_along(colnames(reduced_phenocodes %>%
                                                ungroup() %>%
                                                select(-Mapping))),
                   function(x){
                     marker <- gsub("Phenotype_", "", colnames(reduced_phenocodes)[x])
                     phenocode_vector <- reduced_phenocodes[,x]
                     if(max(as.numeric(phenocode_vector[[1]])) == 3){
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

annotated_mappings <- apply(phenocode_matrix, 1, paste, collapse = " ") %>%
  trimws() %>%
  gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", ., perl = T)
  
##########################################################################
############################     End code     ############################
##########################################################################

workspaceFile <- paste(RESULTS_DIR, "phenotype_gridpoints.Workspace.Rdata", sep = "")

save.image(file = workspaceFile)
