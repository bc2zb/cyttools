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

phenocodes <- read_tsv(args$PHENOCODES)
flowTypeResults <- read_tsv(paste0(args$OUT, "BatchFlowTypeDataMergeResultsTable.txt"))
flowTypeCounts <- read.delim(paste0(args$OUT, "nodeCountFeatureTable.txt"), row.names = 1)
pop_median_counts <- apply(flowTypeCounts, 1, median)

flowTypeCounts <- flowTypeCounts[pop_median_counts > 1000,]

nodeExprTable <- cluster.flowSet.trans %>%
  fsApply(function(x){return(as.data.frame(exprs(x)))}, simplify = F) %>%
  bind_rows(.id = "FileNames") %>%
  group_by_at(vars(FileNames, Mapping)) %>%
  summarise_at(colnames(.)[colnames(.) %in% targets$name[targets$Functional == 1 | targets$Lineage == 1]],
               median) %>%
  gather(Metal,
         Intensity,
         -Mapping, 
         -FileNames) %>%
  spread(FileNames,
         Intensity) %>%
  left_join(targets, by = c("Metal" = "name"))

mappings <- cluster.flowSet.trans %>%
  fsApply(function(x){return(as.data.frame(exprs(x)))}, simplify = F) %>%
  bind_rows(.id = "FileNames") %>% 
  select(Mapping) %>%
  bind_cols(flowTypeResults %>%
              select(Mapping, NodeNames) %>%
              setNames(c("PhenoCodes", "Names")))

map_counts <- mappings %>%
  group_by(Mapping, PhenoCodes, Names) %>%
  summarise(n())

map_counts %>%
  filter(Mapping == 1) %>%
  ungroup() %>%
  mutate(PhenoCodes = factor(PhenoCodes))

map_counts <- map_counts %>%
  mutate(sepCodes = gsub("(.)", "\\1 ", PhenoCodes),
         sepCodes = trimws(sepCodes)) %>%
  separate(sepCodes,
           str_split(trimws(gsub("\\+|\\-", " ", map_counts$Names[1])), " ", simplify = T),
           " ")

map_counts <- map_counts %>%
  group_by(Mapping) %>%
  select(-PhenoCodes,
         -Names,
         -`n()`) %>%
  gather(Marker,
         PhenoCode,
         -Mapping) %>%
  group_by(Mapping, Marker, PhenoCode) %>%
  summarise(n()) %>%
  spread(PhenoCode,
         `n()`,
         fill = 0) %>%
  mutate(diff = (`2` - `1`)/(`2` + `1`))

reduced_phenocodes <- map_counts %>%
  mutate(MarkerCode = case_when(diff >= 0.75 ~ 2,
                                diff <= -0.75 ~ 1,
                                diff < 0.75 & diff > -0.75 ~ 0)) %>%
  select(Mapping, Marker, MarkerCode) %>%
  mutate(MarkerCode = factor(MarkerCode)) %>%
  spread(Marker, MarkerCode) %>%
  select(CD41,CD117,CD13,CD123,CD45RA,CD45,CD36,CD235a,CD61,CD34,CD9,CD38,CD71,Mapping)

reduced_phenocodes$PhenoCodes <- apply(reduced_phenocodes[1:13], 1, paste, collapse = "")

reduced_phenocodes <- reduced_phenocodes %>% 
  left_join(phenocodes)



##########################################################################
############################     End code     ############################
##########################################################################

workspaceFile <- paste(RESULTS_DIR, "phenotype_gridpoints.Workspace.Rdata", sep = "")

save.image(file = workspaceFile)
