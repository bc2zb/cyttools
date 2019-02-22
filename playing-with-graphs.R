graph <- flowSOM.res$MST$graph

load("/Volumes/Home04/CapaldoBj/cyttools_run_FEB_2018/cyttoolsClusteringResults/phenotype_gridpoints.Workspace.Rdata")

map_counts %>%
  arrange(desc(n))

mappings %>%
  group_by_all() %>%
  count() %>%
  arrange(desc(n)) %>%
  head()

unique_filtered_count_table %>%
  rownames_to_column("immunophenotype") %>%
  rowsum()

parent_child <- data.frame(immunophenotype = row.names(unique_filtered_count_table)) %>%
  mutate(phenotype_depth = str_count(immunophenotype, " ")) %>%
  arrange(phenotype_depth)

network_frame <- data.frame()

for(depth in unique(parent_child$phenotype_depth)){
  current_depth <- depth
  deeper_depth <- depth + 1
  if(deeper_depth > max(unique(parent_child$phenotype_depth))){
    break
  }
  
  sep_col_names <- paste0("marker_", c(1:deeper_depth))
  
  parent_pops <- parent_child %>%
    filter(phenotype_depth == current_depth) %>%
    transmute(parent_immunophenotype = immunophenotype) %>%
    separate(parent_immunophenotype, into = sep_col_names, sep = " ", remove = F) %>%
    gather(marker_number,
           immunophenotype,
           starts_with("marker_")) %>%
    select(-marker_number)
  
  sep_col_names <- paste0("marker_", c(0:deeper_depth))
  
  child_pops <- parent_child %>%
    filter(phenotype_depth == deeper_depth) %>%
    transmute(child_immunophenotype = immunophenotype) %>%
    separate(child_immunophenotype, into = sep_col_names, sep = " ", remove = F) %>%
    gather(marker_number,
           immunophenotype,
           starts_with("marker_")) %>%
    select(-marker_number)
  
  parent_child_linked <- parent_pops %>%
    left_join(child_pops) %>%
    count(parent_immunophenotype, child_immunophenotype) %>%
    filter(n == deeper_depth)
 
  network_frame <- bind_rows(network_frame, parent_child_linked) 
}

library(igraph)

test_graph <- graph_from_data_frame(network_frame)  

ggraph(test_graph, 'circlepack', weight = "size") +
  geom_node_circle(aes(fill = n)) +
  coord_fixed()

ggraph(test_graph, layout = 'kk') +
  geom_node_point() + 
  geom_edge_link(color = "grey", alpha = 0.25, arrow = arrow(length = unit(4, 'mm'))) +
  geom_label_repel()

gr <- test_graph

ggraph(test_graph, layout = 'partition') + 
  geom_node_tile()

ggraph(gr, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal() + 
  geom_node_point(aes(filter = leaf)) + 
  coord_fixed()
