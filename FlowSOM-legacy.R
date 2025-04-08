#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
FLowSOM.R (-h | --help | --version)
FLowSOM.R DIR

Description:   This script is a template for making docopts compatible Rscripts
Options:
--version       Show the current version.

Arguments:

DIR    Provide directory for cytools.args.Rdata to be found

" -> doc


args <- docopt(doc)

ARGS_DIR <- args$DIR

cat("\nLoading arguments from", ARGS_DIR, "\n")

load(paste(ARGS_DIR, "cyttools.args.Rdata", sep = ""))

RESULTS_DIR <- args$OUT

source("cyttoolsFunctions.R")

targets <- read.delim(args$PANEL)
colsToCheck <- c("Ignore", "TransformCofactor", "Lineage", "Functional")
if(checkDesignCols(targets, colsToCheck)){
  missingCols <- colsToCheck[which(colsToCheck %in% colnames(targets) == F)]
  cat("\n\nERROR: PANEL file does not include required columns.
      \n\nMissing Columns:", missingCols,
      "\n\nPlease run cyttools.R --makePanelBlank and cyttools.R --computeNRS to generate compatible panel file.\n\nStopping cyttools.R\n\n")
  q()
}

dir <- args$DIR # grabs directory from initial cyttools call
file <- list.files(dir ,pattern='.fcs$', full=TRUE) # captures all FCS files in the directory

targets <- targets %>% mutate(Lineage = if_else(desc %in% c("Time",
                                                 "Event_length",
                                                 "Viability",
                                                 "DNA") | desc == name,
                                        0,
                                        1))

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

fsom <- ReadInput(flowSet.trans, transform = FALSE, scale = FALSE)

colsToUse <- which(targets$name %in% lineage_markers[lineage_markers %in% targets$name[targets$Ignore == 0]] == T)

som <- BuildSOM(fsom,
                colsToUse = targets$name[colsToUse],
                xdim = 23,
                ydim = 23)

flowSOM.res <- BuildMST(som)

ResultsTable <- as.data.frame(flowSOM.res$data)
fileLabels <- flowSOM.res$metaData %>% unlist() %>% matrix(ncol = 2, byrow = T)
fileNames <- vector(length = nrow(ResultsTable))
for(i in 1:nrow(fileLabels)){
  fileNames[fileLabels[i,1]:fileLabels[i,2]] <- file[i]
}

ResultsTable$FileNames <- fileNames
ResultsTable$Mapping <- flowSOM.res$map$mapping[,1]
ResultsTable$DistToNode <- flowSOM.res$map$mapping[,2]

ResultsTable <- ResultsTable %>% left_join(flowSOM.res$MST$l %>%
                                             as.data.frame() %>%
                                             setNames(c("cyttools_dim_x", "cyttools_dim_y")) %>%
                                             rownames_to_column("Mapping") %>%
                                             mutate(Mapping = as.numeric(Mapping)))

#### perform kmeans gating ####

kmeans_filtered_0_results <- lapply(colsToUse, function(column_index){
  median_column_values <- data.frame(msi = ResultsTable[,column_index],
                                     cell_id = factor(c(1:nrow(ResultsTable))))
  filtered_column_values <- median_column_values %>%
    filter(msi > 0)
  kmeans_filter_results <- kmeans(filtered_column_values$msi, 2)
  filtered_column_values <- filtered_column_values %>%
    mutate(cluster = kmeans_filter_results$cluster)
  median_column_values <- median_column_values %>%
    left_join(filtered_column_values) %>%
    tibble() %>%
    mutate(cluster = if_else(is.na(cluster), 0, cluster))
  return(median_column_values)
}) %>%
  lapply(`[[`, "cluster") %>%
  as.data.frame() %>%
  tibble() %>%
  setNames(paste0("kmeans_", targets$name[colsToUse]))

#### find thresholds for different populations ####

threshold_df <- ResultsTable %>%
  select(-c(Time, Event_length, FileNames, cyttools_dim_x, cyttools_dim_y)) %>%
  setNames(paste0("msi_", colnames(.))) %>%
  bind_cols(kmeans_filtered_0_results) %>%
  mutate(map = c(1:nrow(.))) %>%
  pivot_longer(cols = -map,
               names_to = c(".value", "tag"),
               names_pattern = "(.*)_(.*)") %>%
  filter(!is.na(kmeans)) %>%
  group_by(tag, kmeans) %>%
  summarise(mean_msi = mean(msi),
            max_msi = max(msi),
            min_msi = min(msi),
            median_msi = median(msi)) %>%
  ungroup() %>%
  group_by(tag) %>%
  arrange(desc(median_msi)) %>%
  mutate(pos_gate_code = row_number()) %>%
  summarize(threshold = mean(min_msi[pos_gate_code == 1],
                             max_msi[pos_gate_code == 2]))
#### create phenotype matrix and phenotype table ####

phenotype_matrix <- ResultsTable %>%
  select(-c(Time, Event_length, FileNames, cyttools_dim_x, cyttools_dim_y)) %>%
  setNames(paste0("msi_", colnames(.))) %>%
  bind_cols(kmeans_filtered_0_results) %>%
  mutate(map = c(1:nrow(.))) %>%
  pivot_longer(cols = -map,
               names_to = c(".value", "tag"),
               names_pattern = "(.*)_(.*)") %>%
  filter(!is.na(kmeans)) %>%
  left_join(targets %>%
              select(name,desc),
            by = c("tag" = "name")) %>%
  left_join(threshold_df) %>%
  tibble() %>%
  transmute(map = map,
            tag = desc,
            phenotype = if_else(msi > threshold,
                                1L,
                                0L)) %>%
  pivot_wider(names_from = tag,
              values_from = phenotype) %>%
  column_to_rownames("map") %>%
  as.matrix()

phenotype_table <- ResultsTable %>%
  select(-c(Time, Event_length, FileNames, cyttools_dim_x, cyttools_dim_y)) %>%
  setNames(paste0("msi_", colnames(.))) %>%
  bind_cols(kmeans_filtered_0_results) %>%
  mutate(map = c(1:nrow(.))) %>%
  pivot_longer(cols = -map,
               names_to = c(".value", "tag"),
               names_pattern = "(.*)_(.*)") %>%
  filter(!is.na(kmeans)) %>%
  left_join(targets %>%
              select(name,desc),
            by = c("tag" = "name")) %>%
  left_join(threshold_df) %>%
  tibble() %>%
  transmute(map = map,
            tag = desc,
            phenotype = if_else(msi > threshold,
                                paste0(desc, "+"),
                                paste0(desc, "-"))) %>%
  pivot_wider(names_from = tag,
              values_from = phenotype)

#### perform immunophenotyping ####

immunophenotype_list <- read_csv("immunophenotypes-profiling-database-2024-may-08.csv") %>%
  select(-notes) %>%
  distinct(cell_name, lineage_markers) %>%
  filter(cell_name != "Root") %>%
  separate_longer_delim(lineage_markers, " ") %>%
  transmute(cell_type = cell_name,
    cell_type_index = 1,
    Marker = str_remove(lineage_markers, "\\+|\\-"),
    Marker_status = lineage_markers) %>%
  mutate(Marker = str_replace_all(Marker, "HLA-DR", "HLADR"),
         Marker_status = str_replace_all(Marker_status, "HLA-DR", "HLADR")) %>%
  filter(Marker %in% colnames(phenotype_table)) %>%
  mutate(population_desc = paste(cell_type, cell_type_index, sep = "_")) %>%
  split(.$population_desc) %>%
  lapply(pivot_wider, names_from = Marker,
                     values_from = Marker_status) %>%
  lapply(function(df){
    col_labels <- colnames(df)
    colnames(df)[1] <- df$cell_type[1]
    return(df %>%
             select(-c(cell_type_index, population_desc)))
  })

phenotyped_list <- immunophenotype_list %>%
  lapply(left_join, phenotype_table)

phenotyped_table <- phenotyped_list %>%
  lapply(transmute, map = map, val = 1L) %>%
  bind_rows(.id = "cell_desc") %>%
  pivot_wider(names_from = cell_desc,
              values_from = val,
              values_fill = 0L)

profile_list <- read_csv("immunophenotypes-profiling-database-2024-may-08.csv") %>%
  select(-notes) %>%
  distinct(CellSubset, cell_name, functional_markers) %>%
  filter(cell_name != "Root") %>%
  separate_longer_delim(functional_markers, " ") %>%
  transmute(cell_type = cell_name,
    CellSubset = CellSubset,
    cell_type_index = 1,
    Marker = str_remove(functional_markers, "\\hi|\\lo") %>%
        str_remove("\\_"),
    Marker_status = str_replace(functional_markers, "hi", "+") %>%
        str_replace("lo", "-")) %>%
  mutate(Marker = str_replace_all(Marker, "HLA-DR", "HLADR"),
         Marker_status = str_replace_all(Marker_status, "HLA-DR", "HLADR")) %>%
  mutate(population_desc = paste(cell_type, cell_type_index, sep = "_")) %>%
  split(.$CellSubset) %>%
  lapply(pivot_wider, names_from = Marker,
                     values_from = Marker_status) %>%
  lapply(function(df){
    col_labels <- colnames(df)
    colnames(df)[1] <- df$cell_type[1]
    return(df %>%
             select(-c(cell_type_index, population_desc)))
  })

profiled_table <- lapply(phenotyped_list, function(phenotyped_df){
    cell_type_name <- phenotyped_df[1,1]
    profile_cell_type_names <- lapply(profile_list, function(profile_df){
        profile_df[,1][[1]][1]
    }) %>%
    unlist()
    return(profile_list[names(profile_cell_type_names[profile_cell_type_names %in% cell_type_name[,1][[1]][1]])] %>%
        lapply(left_join, phenotyped_df) %>%
        lapply(transmute, map = map, val = 1L) %>%
        bind_rows(.id = "cell_desc")) 
}) %>%
    bind_rows() %>%
    filter(!is.na(map)) %>%
    pivot_wider(names_from = cell_desc,
                        values_from = val,
                        values_fill = 0L)

compartment_table <- phenotyped_table %>%
    select(map, `B Cell_1`, `T Cell_1`, `Myeloid_1`, `NK Cell_1`, `Granulocyte Basophil_1`)

#### write out results ####
dir.create(paste0(RESULTS_DIR, "CLUSTERED_FCS/"),
           showWarnings = F)
# for (files in file)
lapply(file, function(files){
  cat("processing ", files, "\n")
  rawFCS <- read.FCS(files, transformation = F)
  clusterData <- ResultsTable %>%
    mutate(map = c(1:nrow(.))) %>%
    dplyr::filter(FileNames == files) %>%
    select(map, Mapping, DistToNode, cyttools_dim_x, cyttools_dim_y) %>%
    mutate(root_unassigned = if_else(map %in% phenotyped_table$map, 0L, 1L)) %>%
    left_join(phenotyped_table, by = "map") %>%
    left_join(phenotype_matrix |>
      as.data.frame() |>
      setNames(paste(colnames(phenotype_matrix), "gate", sep = ".")) |>
      mutate(map = c(1:nrow(.))),
      by = "map") |>
    select(-map) %>%
    mutate(across(all_of(colnames(phenotyped_table)[-1]),
      ~if_else(is.na(.x), 0, .x)))
  
  # assignment.csv (Rows are channels, columns are population descriptions, values are median signal intensity.)

  ResultsTable %>%
    mutate(map = c(1:nrow(.))) %>%
    tibble() %>%
    dplyr::filter(FileNames == files) %>%
    select(map, all_of(colsToUse)) %>%
    mutate(root_unassigned = if_else(map %in% phenotyped_table$map, 0L, 1L)) %>%
    left_join(phenotyped_table, by = "map") %>%
    select(-map) %>%
    mutate(across(all_of(colnames(phenotyped_table)[-1]),
      ~if_else(is.na(.x), 0, .x))) %>%
    pivot_longer(ends_with("Di"),
        names_to = "name",
        values_to = "MSI"
        ) %>%
    pivot_longer(-c(name, MSI),
        names_to = "CellSubset",
        values_to = "subset_status") %>%
    filter(subset_status == 1) %>%
    group_by(name, CellSubset) %>%
    summarise(MSI = median(MSI)) %>%
    ungroup() %>%
    left_join(targets %>%
        select(name, descOriginal) %>%
        setNames(c("name", "Channel")),
        by = "name") %>%
    select(-name) %>%
    mutate(CellSubset = str_remove(CellSubset, "\\_1$")) %>%
    pivot_wider(names_from = CellSubset,
        values_from = "MSI") %>%
    write_csv(paste0(RESULTS_DIR,
                     "CLUSTERED_FCS/clustered_",
                     str_replace(basename(files),
                                 "\\.fcs$",
                                 "-assignment.csv")))
  cat("assignment.csv written successfully \n")
# cell_counts_assignment.csv (Rows are cell subsets, column is cell counts)

  clusterData %>%
    select(-c(Mapping:cyttools_dim_y)) %>%
    select(-ends_with(".gate")) |>
    colSums() %>%
    as.data.frame() %>%
    setNames("count") %>%
    rownames_to_column("CellSubset") %>%
    transmute(CellSubset = str_remove(CellSubset, "\\_1$"),
      N = count) %>%
    write_csv(paste0(RESULTS_DIR,
                     "CLUSTERED_FCS/clustered_",
                     str_replace(basename(files),
                                 "\\.fcs$",
                                 "-cell_counts_assignments.csv")))
  cat("cell_counts_assignment.csv written successfully \n")
# cell_counts_compartment.csv (Rows are cell compartments, column is cell counts)

  ResultsTable %>%
    mutate(map = c(1:nrow(.))) %>%
    tibble() %>%
    dplyr::filter(FileNames == files) %>%
    select(map, all_of(colsToUse)) %>%
    #mutate(root_unassigned = if_else(map %in% compartment_table$map, 0L, 1L)) %>%
    left_join(compartment_table, by = "map") %>%
    select(-map) %>%
    mutate(across(all_of(colnames(compartment_table)[-1]),
      ~if_else(is.na(.x), 0, .x)),
      root_unassigned = rowSums(across(all_of(colnames(compartment_table)[-1])))) %>%
    select(c(root_unassigned, all_of(colnames(compartment_table)[-1]))) %>%
    mutate(root_unassigned = if_else(root_unassigned == 0, 1, 0)) |>
    colSums() %>%
    as.data.frame() %>%
    setNames("count") %>%
    rownames_to_column("CellSubset") %>%
    transmute(CellSubset = str_remove(CellSubset, "\\_1$"),
      N = count) %>%
    write_csv(paste0(RESULTS_DIR,
                     "CLUSTERED_FCS/clustered_",
                     str_replace(basename(files),
                                 "\\.fcs$",
                                 "-cell_counts_compartment.csv")))
  cat("cell_counts_compartment.csv written successfully \n")
# compartment.csv (Rows are channels, columns are compartment descriptions, values are median signal intensity.)

  ResultsTable %>%
    mutate(map = c(1:nrow(.))) %>%
    tibble() %>%
    dplyr::filter(FileNames == files) %>%
    select(map, all_of(colsToUse)) %>%
    mutate(root_unassigned = if_else(map %in% compartment_table$map, 0L, 1L)) %>%
    left_join(compartment_table, by = "map") %>%
    select(-map) %>%
    mutate(across(all_of(colnames(compartment_table)[-1]),
      ~if_else(is.na(.x), 0, .x))) %>%
    pivot_longer(ends_with("Di"),
        names_to = "name",
        values_to = "MSI"
        ) %>%
    pivot_longer(-c(name, MSI),
        names_to = "CellSubset",
        values_to = "subset_status") %>%
    filter(subset_status == 1) %>%
    group_by(name, CellSubset) %>%
    summarise(MSI = median(MSI)) %>%
    ungroup() %>%
    left_join(targets %>%
        select(name, descOriginal) %>%
        setNames(c("name", "Channel")),
        by = "name") %>%
    select(-name) %>%
    mutate(CellSubset = str_remove(CellSubset, "\\_1$")) %>%
    pivot_wider(names_from = CellSubset,
        values_from = "MSI") %>%
    write_csv(paste0(RESULTS_DIR,
                     "CLUSTERED_FCS/clustered_",
                     str_replace(basename(files),
                                 "\\.fcs$",
                                 "-compartment.csv")))
  cat("compartment.csv written successfully \n")
# profiling.csv (Rows are channels, columns are profiled cell subsets, values are median signal intensity)

  ResultsTable %>%
    mutate(map = c(1:nrow(.))) %>%
    tibble() %>%
    dplyr::filter(FileNames == files) %>%
    select(map, all_of(colsToUse)) %>%
    mutate(root_unassigned = if_else(map %in% profiled_table$map, 0L, 1L)) %>%
    left_join(profiled_table, by = "map") %>%
    select(-map) %>%
    mutate(across(all_of(colnames(profiled_table)[-1]),
      ~if_else(is.na(.x), 0, .x))) %>%
    pivot_longer(ends_with("Di"),
        names_to = "name",
        values_to = "MSI"
        ) %>%
    pivot_longer(-c(name, MSI),
        names_to = "CellSubset",
        values_to = "subset_status") %>%
    filter(subset_status == 1) %>%
    group_by(name, CellSubset) %>%
    summarise(MSI = median(MSI)) %>%
    ungroup() %>%
    left_join(targets %>%
        select(name, descOriginal) %>%
        setNames(c("name", "Channel")),
        by = "name") %>%
    select(-name) %>%
    mutate(CellSubset = str_remove(CellSubset, "\\_1$")) %>%
    pivot_wider(names_from = CellSubset,
        values_from = "MSI") %>%
    write_csv(paste0(RESULTS_DIR,
                     "CLUSTERED_FCS/clustered_",
                     str_replace(basename(files),
                                 "\\.fcs$",
                                 "-profiling.csv")))
  cat("profiling.csv written successfully \n")
# cell_counts_profiling.csv (Rows are cell population profiles, column is cell counts)
  ResultsTable %>%
    mutate(map = c(1:nrow(.))) %>%
    tibble() %>%
    dplyr::filter(FileNames == files) %>%
    select(map, all_of(colsToUse)) %>%
    mutate(root_unassigned = if_else(map %in% profiled_table$map, 0L, 1L)) %>%
    left_join(profiled_table, by = "map") %>%
    select(-map) %>%
    mutate(across(all_of(colnames(profiled_table)[-1]),
      ~if_else(is.na(.x), 0, .x))) %>%
    select(all_of(colnames(profiled_table)[-1])) %>%
    colSums() %>%
    as.data.frame() %>%
    setNames("count") %>%
    rownames_to_column("CellSubset") %>%
    transmute(CellSubset = str_remove(CellSubset, "\\_1$"),
      N = count) %>%
    write_csv(paste0(RESULTS_DIR,
                     "CLUSTERED_FCS/clustered_",
                     str_replace(basename(files),
                                 "\\.fcs$",
                                 "-cell_counts_profiling.csv")))
  cat("cell_counts_profiling.csv written successfully \n")
# Granulocytes should be the combination of eosinophils, neutrophils, and basophils (which we have definitions for)
# Other is leftovers

  clusterFCS <- fr_append_cols(rawFCS, as.matrix(clusterData))
  row.names(pData(parameters(clusterFCS))) <- paste0("$P", c(1:nrow(pData(parameters(clusterFCS)))))
  out.fcs.file <- paste0(RESULTS_DIR, "CLUSTERED_FCS/clustered_", basename(files))
  write.FCS(clusterFCS, out.fcs.file)
})

# # write out results
# ResultsTableFile <- paste(RESULTS_DIR, "FlowSOMResultsTable.txt", sep = "")
# nodeExprTableFile <- paste(RESULTS_DIR, "nodeExpressionFeatureTable.txt", sep = "")
# nodeAbndncFeatureTableFile <- paste(RESULTS_DIR, "nodeAbundanceFeatureTable.txt", sep = "")
# nodeCountFeatureTableFile <- paste(RESULTS_DIR, "nodeCountFeatureTable.txt", sep = "")
# nodeMedianFeatureTableFile <- paste(RESULTS_DIR, "nodeMedianFeatureTable.txt", sep = "")
# 
# write.table(ResultsTable, ResultsTableFile, sep = "\t", quote = F, row.names = F)
# write.table(nodeExprTable, nodeExprTableFile, sep = "\t", quote = F, row.names = T)
# write.table(props_table, nodeAbndncFeatureTableFile, sep = "\t", quote = F, row.names = T)
# write.table(countTable, nodeCountFeatureTableFile, sep = "\t", quote = F, row.names = T)
# write.table(flowSOM.res$map$medianValues, nodeMedianFeatureTableFile, sep = "\t", quote = F, row.names = T)

# workspaceFile <- paste(RESULTS_DIR, "FlowSOMWorkspace.Rdata", sep = "")
# 
# save.image(file = workspaceFile)
