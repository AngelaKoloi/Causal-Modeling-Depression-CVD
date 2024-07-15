library(igraph)

convert_to_igraph <- function(graph) {
  adj_matrix <- as(graph@graph, "matrix")
  igraph_obj <- graph_from_adjacency_matrix(adj_matrix, mode = "directed")
  return(igraph_obj)
}

compare_graph_edges <- function(original, bootstraps, description) {
  edges_original <- get.edgelist(original)
  comparison_results <- lapply(bootstraps, function(bootstrap) {
    edges_bootstrap <- get.edgelist(bootstrap)
    common_edges <- intersect(paste(edges_original[, 1], edges_original[, 2], sep = "-"), 
                              paste(edges_bootstrap[, 1], edges_bootstrap[, 2], sep = "-"))
    data.frame(
      Total_Edges_Original = nrow(edges_original),
      Total_Edges_Bootstrap = nrow(edges_bootstrap),
      Common_Edges = length(common_edges),
      Common_Percentage = length(common_edges) / nrow(edges_original) * 100
    )
  })
  comparison_summary <- do.call(rbind, comparison_results)
  comparison_summary$Comparison <- description
  return(comparison_summary)
}

perform_path_analysis <- function(graph, start_node, end_node) {
  paths <- all_simple_paths(graph, from = start_node, to = end_node)
  if (length(paths) == 0) {  # Check if no paths were found
    return(list(
      Shortest_Path = NA,
      Most_Frequent_Path = NA,
      Most_Frequent_Nodes = NA
    ))
  }
  
  path_strings <- lapply(paths, function(path) paste(V(graph)[path]$name, collapse = "->"))
  
  shortest_path <- which.min(sapply(paths, length))
  shortest_path_string <- path_strings[[shortest_path]]
  
  path_table <- table(unlist(path_strings))
  most_frequent_path <- names(which.max(path_table))
  
  node_frequency <- table(unlist(lapply(paths, function(path) V(graph)[path]$name)))
  most_frequent_nodes <- names(sort(node_frequency, decreasing = TRUE)[1:3])
  
  return(list(
    Shortest_Path = shortest_path_string,
    Most_Frequent_Path = most_frequent_path,
    Most_Frequent_Nodes = most_frequent_nodes
  ))
}

load("causal_graph_0.05_with_whitelist.RData")
original_igraph <- convert_to_igraph(original_graph)

load("bootstrap_graph.RData")
bootstrap_igraph_0 <- lapply(results_0, convert_to_igraph)

load("bootstrap_graph_10.RData")
bootstrap_igraph_10 <- lapply(results_10, convert_to_igraph)

load("bootstrap_graph_20.RData")
bootstrap_igraph_20 <- lapply(results_20, convert_to_igraph)

# Compare edges
edge_comparison_results_0 <- compare_graph_edges(original_igraph, bootstrap_igraph_0, "Original vs. Bootstrap 0%")
edge_comparison_results_10 <- compare_graph_edges(original_igraph, bootstrap_igraph_10, "Original vs. Bootstrap 10%")
edge_comparison_results_20 <- compare_graph_edges(original_igraph, bootstrap_igraph_20, "Original vs. Bootstrap 20%")
all_comparisons <- rbind(edge_comparison_results_0, edge_comparison_results_10, edge_comparison_results_20)
write.csv(all_comparisons, "edge_comparison_results.csv", row.names = FALSE)


path_analysis1 <- perform_path_analysis(original_igraph, 'acidep09', 'ecidep09')
path_analysis2 <- perform_path_analysis(original_igraph, 'acidep09', 'amet_syn2')
path_analysis3 <- perform_path_analysis(original_igraph, 'acidep09', 'emet_syn2')
path_analysis4 <- perform_path_analysis(original_igraph, 'amet_syn2', 'emet_syn2')

# Combine results into a single data frame
combined_results <- rbind(data.frame(path='depression(B)_to_depression(FU)', path_analysis1),
                          data.frame(path='depression(B)_to_metS(B)', path_analysis2),
                          data.frame(path='depression(B)_to_metS(FU)', path_analysis3),
                          data.frame(path='metS(B)_to_metS(FU)', path_analysis4))

write.csv(combined_results, "path_analysis_results.csv", row.names = FALSE)

