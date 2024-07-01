library(visNetwork)
library(glue)
library(igraph)
library(htmlwidgets)

load("causal_graph_test.RData")
adj_matrix <- as(graph@graph, "matrix")
g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed")
nodes <- data.frame(id = V(g)$name, label = V(g)$name)
print(graph)
tier1 <- c('Age', 'Sexe')
tier2 <- c('aedu', 'asmokstat', 'AIPMETO2', 'aauditsc', 'aIRSsum9', 'abaiscal', 'aids', 'acidep09', 'amet_syn2', 'ams_waist', 'ams_hpt', 'ams_trig2', 'ams_hdl2', 'ams_gluc2', 'atri_med', 'ahdl_med', 'asbp_med', 'adbp_med', 'agluc_med', 'ahsCRP', 'aIL6', 'aApoB', 'aHDL_C', 'aTotFA', 'aSerum_TG', 'aGp', 'aIle')
tier3 <- c('eage', 'sex')
tier4 <- c('eipmeto2', 'eauditsc', 'eIRSsum9', 'ebaiscal', 'eids', 'ecidep09', 'emet_syn2', 'ems_waist', 'ems_hpt', 'ems_trig2', 'ems_hdl2', 'ems_gluc2', 'etri_med', 'ehdl_med', 'esbp_med', 'edbp_med', 'egluc_med', 'eHSCRP', 'eIL6', 'eApoB', 'eHDLC', 'eTotFA', 'eSerumTG', 'eGp', 'eIle')

renamed_labels <- list(
  'Age' = 'Age', 'eage' = 'Age',
  'Sexe' = 'Gender', 'sex' = 'Gender',
  'aedu' = 'Education',
  'asmokstat' = 'Smoking status',
  'aauditsc' = 'Alcohol consumption', 'eauditsc' = 'Alcohol consumption',
  'AIPMETO2' = 'Physical activity', 'eipmeto2' = 'Physical activity',
  'aIRSsum9' = 'Sleep pattern', 'eIRSsum9' = 'Sleep pattern',
  'acidep09' = 'Major depression', 'ecidep09' = 'Major depression',
  'amet_syn2' = 'MetS', 'emet_syn2' = 'MetS',
  'ams_waist' = 'Obesity', 'ems_waist' = 'Obesity',
  'ams_hpt' = 'Hypertension', 'ems_hpt' = 'Hypertension',
  'ams_trig2' = 'Hypertriglyceridemia', 'ems_trig2' = 'Hypertriglyceridemia',
  'ams_hdl2' = 'Low HDL cholesterol', 'ems_hdl2' = 'Low HDL cholesterol',
  'ams_gluc2' = 'Hyperglycemia', 'ems_gluc2' = 'Hyperglycemia',
  'atri_med' = 'TG', 'etri_med' = 'TG',
  'ahdl_med' = 'HDL', 'ehdl_med' = 'HDL',
  'asbp_med' = 'Systolic BP', 'esbp_med' = 'Systolic BP',
  'adbp_med' = 'Diastolic BP', 'edbp_med' = 'Diastolic BP',
  'agluc_med' = 'Glucose', 'egluc_med' = 'Glucose',
  'ahsCRP' = 'hs-CRP', 'eHSCRP' = 'hs-CRP',
  'aIL6' = 'IL-6', 'eIL6' = 'IL-6',
  'aApoB' = 'ApoB', 'eApoB' = 'ApoB',
  'aHDL_C' = 'HDLC', 'eHDLC' = 'HDLC',
  'aTotFA' = 'Total FA', 'eTotFA' = 'Total FA',
  'aSerum_TG' = 'Serum TG', 'eSerumTG' = 'Serum TG',
  'aGp' = 'AGP', 'eGp' = 'AGP',
  'aIle' = 'Ile', 'eIle' = 'Ile',
  'abaiscal' = 'Anxiety', 'ebaiscal' = 'Anxiety',
  'aids' = 'Depression severity', 'eids' = 'Depression severity'
)

# Rename node labels
nodes$label <- sapply(nodes$label, function(x) renamed_labels[[x]])

# Assign node shapes
baseline_variables <- c('Age', 'Sexe', 'aedu', 'asmokstat', 'AIPMETO2', 'aauditsc', 'aIRSsum9', 'aids', 
                        'abaiscal', 'acidep09', 'amet_syn2', 'ams_waist', 'ams_hpt', 'ams_trig2', 
                        'ams_hdl2', 'ams_gluc2', 'atri_med', 'ahdl_med', 'asbp_med', 'adbp_med', 
                        'agluc_med', 'ahsCRP', 'aIL6', 'aApoB', 'aHDL_C', 'aTotFA', 'aSerum_TG', 
                        'aGp', 'aIle')
follow_up_variables <- c('eage', 'sex', 'eipmeto2', 'eauditsc', 'eIRSsum9', 'ebaiscal', 'eids', 
                         'ecidep09', 'emet_syn2', 'ems_waist', 'ems_hpt', 'ems_trig2', 'ems_hdl2', 
                         'ems_gluc2', 'etri_med', 'ehdl_med', 'esbp_med', 'edbp_med', 'egluc_med', 
                         'eHSCRP', 'eIL6', 'eApoB', 'eHDLC', 'eTotFA', 'eSerumTG', 'eGp', 'eIle')

# Assign shapes based on baseline or follow-up
nodes$shape <- ifelse(nodes$id %in% baseline_variables, 'square', 'dot')

# Define colors for different variable groups
color_map <- list(
  'Background variables' = '#ADD7F6',
  'Conditional variables' = '#2667FF',
  'Metabolites' = '#90EE90',
  'Intermediate variables' = '#ffb703',
  'Diagnostic variables' = '#FF6347'
)

variable_groups <- list(
  'Background variables' = c('Age', 'Gender'),
  'Metabolites' = c('ApoB', 'Total FA', 'Serum TG', 'AGP', 'Ile', 'HDLC'),
  'Diagnostic variables' = c('Major depression', 'MetS'),
  'Conditional variables' = c('Smoking status', 'Physical activity', 'Education', 'Alcohol consumption', 'Sleep pattern', 'Anxiety', 'Depression severity')
)

all_variables <- unique(unlist(renamed_labels))
assigned_variables <- unlist(variable_groups)
variable_groups$`Intermediate variables` <- setdiff(all_variables, assigned_variables)

nodes$color <- sapply(nodes$label, function(x) {
  group <- names(variable_groups)[sapply(variable_groups, function(v) x %in% v)]
  if (length(group) > 0) {
    color_map[[group]]
  } else {
    '#FFD700' 
  }
})

assign_positions <- function(nodes, max_nodes_per_row = 8, row_height = 150) {
  positions <- data.frame(id = nodes$id, x = NA, y = NA, stringsAsFactors = FALSE)
  current_y <- 0
  total_nodes <- nrow(nodes)
  
  for (i in 1:ceiling(total_nodes / max_nodes_per_row)) {
    row_nodes <- nodes$id[((i - 1) * max_nodes_per_row + 1):min(i * max_nodes_per_row, total_nodes)]
    num_nodes <- length(row_nodes)
    
    x_start <- -(num_nodes - 1) * 75
    
    for (j in 1:num_nodes) {
      node_id <- row_nodes[j]
      positions[positions$id == node_id, "x"] <- x_start + (j - 1) * 150
      positions[positions$id == node_id, "y"] <- current_y
    }
    
    current_y <- current_y + row_height
  }
  
  return(positions)
}
positions <- assign_positions(nodes)

nodes <- merge(nodes, positions, by = "id")
edges <- data.frame(from = as.character(as_edgelist(g)[, 1]),
                    to = as.character(as_edgelist(g)[, 2]))
print(graph)

# Define legend for shapes and colors
legend_nodes <- data.frame(
  label = c('B', 'FU', names(color_map)),
  shape = c('square', 'dot', rep('dot', length(color_map))),
  color = c(NA, NA, unlist(color_map)),
  stringsAsFactors = FALSE
)

net <- visNetwork(nodes, edges) %>%
  visEdges(arrows = 'to') %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = FALSE) %>%
  visLayout(randomSeed = 123) %>%
  visNodes(fixed = TRUE) %>%
  visLegend(addNodes = legend_nodes, useGroups = FALSE)

saveWidget(net, file = "demo.html", selfcontained = TRUE)
  