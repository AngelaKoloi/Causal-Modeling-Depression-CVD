library(R6)

NetworkAnalysis <- R6Class("NetworkAnalysis",
    public = list(
        dataset = NULL,
        n_boots = NULL,
        tensor = NULL,
        network = NULL,
        save_location = NULL,
        image_width = NULL,
        image_height = NULL,
        initialize = function(dataset, save_location, n_boots = 1000, seed = 149,
                              image_width = 1200, image_height = 800) {
            self$dataset <- as.data.frame(dataset)
            self$n_boots <- n_boots
            self$get_imports()
            self$save_location <- save_location

            self$image_width <- image_width
            self$image_height <- image_height

            set.seed(seed)

            # self$build_coef()
        },
        set_coef = function(coef) {
            self$tensor <- coef
        },
        get_imports = function() {
            packages <- c(
                "glmnet", "boot", "progress", "qgraph", "bootnet",
                "ggplot2", "forcats", "dplyr", "gridExtra", "stringr",
                "igraph"
            )
            for (pkg in packages) {
                if (!require(pkg, character.only = TRUE)) {
                    install.packages(pkg)
                    library(pkg, character.only = TRUE)
                }
            }
        },
        get_basic_descriptive = function(round_decimal = 2) {
            df <- data.frame(do.call(cbind, lapply(self$dataset, summary)))
            df <- t(df)
            df <- round(df, 2)

            N_vector <- data.frame(N = rep(dim(self$dataset)[1], dim(df)[1]))
            df <- cbind(df, N_vector)

            temp_data <- self$dataset[complete.cases(self$dataset), ]
            modes <- sapply(temp_data, private$helper_find_mode)
            mode_vector <- data.frame(Mode = modes)

            sd <- sapply(temp_data, sd)
            sd_vector <- data.frame(SD = sd)

            df <- cbind(df, mode_vector, sd_vector)
            df <- df[
                order(rownames(df)),
                c("N", "1st Qu.", "3rd Qu.", "Min.", "Max.", "Mean", "Median", "Mode", "SD")
            ]

            return(df)
        },
        build_coef = function() {
            private$start_time()

            timestep_list <- private$get_timesteps_list()
            ncol <- ncol(self$dataset) / length(timestep_list)
            nrow <- ncol
            n_timesteps <- length(timestep_list) - 1

            matrix_list <- lapply(1:n_timesteps, function(timestep) {
                private$calc_coef(timestep)
            })

            self$tensor <- array(unlist(matrix_list), dim = c(nrow, ncol, length(matrix_list)))
            private$stop_time()

            return(self$tensor)
        },
        plot_functions = function() {
            self$build_coef()
            self$plot_networks()
            self$plot_centralities()
            self$plot_boot_weights()
            self$plot_boot_centralities()
        },
        plot_descriptive = function() {
            df <- get_basic_descriptive()
            table <- tableGrob(df)
            ggsave(paste0(self$save_location, "descriptive.png"),
                plot = table, width = self$image_width, height = self$image_height, dpi = 300
            )
        },
        plot_networks = function() {
            for (i in 1:dim(self$tensor)[3]) {
                graph <- self$tensor[, , i]
                network <- private$calc_network(graph, i)
                for (self_edge in c(TRUE, FALSE)) {
                    if (self_edge) {
                        colnames(network$graph) <- sapply(
                            network$colNames,
                            function(x) sub("_.*", "", x)
                        )
                        network$edge_type <- " Network "
                        image_name <- paste0("network_", network$years[1], network$years[2], ".png")
                    } else {
                        network <- private$remove_self_edge(network)
                        network$edge_type <- " Network Simplified "
                        image_name <- paste0("network_", network$years[1], network$years[2], "_simplified.png")
                    }
                    private$helper_plot_network(network, image_name, edge_width = 8)
                    return(network)
                }
            }
        },
        plot_centralities = function(decreasing = TRUE,
                                     legend_name = "Centrality metric for each timestep") {
            network_list <- list()
            for (self_edge in c(TRUE, FALSE)) {
                for (i in 1:dim(self$tensor)[3]) {
                    graph <- self$tensor[, , i]
                    network <- private$calc_network(graph, i)

                    if (self_edge) {
                        colnames(network$graph) <- sapply(
                            network$colNames,
                            function(x) sub("_.*", "", x)
                        )

                        indicator <- 2
                        network$edge_type <- ""
                    } else {
                        indicator <- 0
                        network <- private$remove_self_edge(network)
                        network$edge_type <- "simplified"
                    }
                    network$name <- paste(network$year[1], "→", network$year[2], network$edge_type)
                    network_list[[i + indicator]] <- network
                }
            }
            private$helper_plot_centrality(network_list, legend_name)
        },
        plot_boot_weights = function(n_cores = 8) {
            for (self_edge in c(TRUE, FALSE)) {
                for (i in 1:dim(self$tensor)[3]) {
                    graph <- self$tensor[, , i]
                    network <- private$calc_network(graph, i)

                    if (self_edge) {
                        colnames(network$graph) <- sapply(
                            network$colNames,
                            function(x) sub("_.*", "", x)
                        )
                        network$edge_type <- " Network "
                        image_name <- paste0("bootWeights_", network$years[1], network$years[2], "_selfEdges.png")
                    } else {
                        network <- private$remove_self_edge(network)
                        network$edge_type <- " Network Simplified "
                        image_name <- paste0("bootWeights_", network$years[1], network$years[2], "_NoselfEdges.png")
                    }
                    private$helper_plot_boot_weights(network, image_name, n_cores)
                }
            }
        },
        plot_boot_centralities = function() {
            for (self_edge in c(TRUE, FALSE)) {
                for (i in 1:dim(self$tensor)[3]) {
                    graph <- self$tensor[, , i]
                    network <- private$calc_network(graph, i)

                    if (self_edge) {
                        colnames(network$graph) <- sapply(
                            network$colNames,
                            function(x) sub("_.*", "", x)
                        )
                        network$edge_type <- " Network "
                        image_name <- paste0("bootCentralities_", network$years[1], network$years[2], "_selfEdges.png")
                    } else {
                        network <- private$remove_self_edge(network)
                        network$edge_type <- " Network Simplified "
                        image_name <- paste0("bootCentralities_", network$years[1], network$years[2], "_NoselfEdges.png")
                    }
                    private$helper_plot_boot_centralities(network, image_name, n_cores)
                    return(network)
                }
            }
        }
    ),
    private = list(
        full_name_labels = c(
            "Sadness", "Pessimism", "Past Failure", "Loss Of Pleasure", "Guilty Feelings",
            "Punishment Feelings", "Self Dislike", "Self Criticalness", "Suicidal Thought Or Wishes",
            "Crying", "Agitation", "Loss Of Interest", "Indecisiveness", "Worthlessness",
            "Loss Of Energy", "Changes In Sleep Pattern", "Irritability", "Changes In Appetite",
            "Concentration Difficulty", "Tiredness Or Fatigue", "Loss Of Interest In Sex",
            "Acetate", "Apoprotein", "CRP", "diastolic KV", "GLUK", "Cholesterol HDL",
            "INSU", "Cholesterol LDL", "Systolic Blood Pressure", "Cholesterol Total", "Triglycerides"
        ),
        var_start_time = NULL,
        var_stop_time = NULL,
        start_time = function() {
            private$var_start_time <- Sys.time()
        },
        stop_time = function() {
            private$var_stop_time <- Sys.time()
            cat("time spent: ", private$var_stop_time - private$var_start_time, "\n")
        },
        get_timesteps_list = function() {
            column_names <- colnames(self$dataset)
            timesteps_list <- unique(as.numeric(
                regmatches(column_names, regexpr("\\d{4}", column_names))
            ))

            timesteps_list <- sort(as.numeric(timesteps_list))
            return(timesteps_list)
        },
        get_subset_by_timestep = function(data, timestep) {
            subset_data <- data[, grep(paste0(timestep), names(data), value = TRUE)]
            return(subset_data)
        },
        get_min_and_max_edge = function(tensor) {
            #' Function to calculate the min and max values in the tensor
            #'
            #' @param tensor A tensor (a multi-dimensional array) of numeric values.
            #' @return A vector of length two, where the first element is the minimum value in the tensor and the second element is the maximum value.

            return(c(min(tensor), max(tensor)))
        },
        get_suffix = function(name) {
            str_extract(name, "_[^_]+$")
        },
        remove_self_edge = function(network, remove_self_edge = FALSE) {
            colnames(network$graph) <- network$colNames
            row_suffixes <- sapply(rownames(network$graph), private$get_suffix)
            col_suffixes <- sapply(colnames(network$graph), private$get_suffix)

            for (i in 1:nrow(network$graph)) {
                for (j in 1:ncol(network$graph)) {
                    if (i == j && remove_self_edge) {
                        break
                    }
                    if (row_suffixes[i] == col_suffixes[j]) {
                        network$graph[i, j] <- 0
                    }
                }
            }
            colnames(network$graph) <- sapply(
                network$colNames,
                function(x) sub("_.*", "", x)
            )
            return(network)
        },
        calc_coef = function(timestep) {
            timesteps_list <- private$get_timesteps_list()
            current_timestep <- timesteps_list[timestep]
            next_timestep <- timesteps_list[timestep + 1]

            pb <- progress_bar$new(
                format = paste0(
                    "Timestep: ", current_timestep, "-", next_timestep,
                    " Bootstrapping [:bar] :percent in :elapsed"
                ),
                total = self$n_boots,
                clear = FALSE,
                width = 60
            )

            boot_result <- boot(self$dataset, statistic = function(data, indices) {
                if (!pb$finished) {
                    pb$tick()
                }
                private$boot_fn(data, indices, current_timestep, next_timestep)
            }, R = self$n_boots)

            booted_coefficients_matrix <- apply(boot_result$t, 2, mean)

            return(matrix(booted_coefficients_matrix,
                nrow = ncol(private$get_subset_by_timestep(self$dataset, current_timestep))
            ))
        },
        calc_network = function(graph, i) {
            node_adverbiations <- c(
                "Sad", "Pes", "PaF", "LOP", "GuF", "PuF", "SDl", "SCr", "Sui",
                "Cry", "Agi", "LOI", "Ind", "Wor", "LOE", "CSP", "Irr", "CIA",
                "Cod", "TOF", "LIS", "Ace", "Apo", "CRP", "DKV", "GLU", "CHDL",
                "INS", "CLDL", "SBP", "CHT", "Tri"
            )

            timestep_list <- private$get_timesteps_list()
            year_t1 <- timestep_list[i]
            year_t2 <- timestep_list[i + 1]

            data_t1 <- private$get_subset_by_timestep(self$dataset, year_t1)
            data_t2 <- private$get_subset_by_timestep(self$dataset, year_t2)

            network <- estimateNetwork(graph, default = "IsingFit")
            network$labels <- node_adverbiations
            network$graph <- graph
            network$rowNames <- colnames(data_t1)
            network$colNames <- colnames(data_t2)
            network$years <- c(year_t1, year_t2)
            
            colnames(data_t1) <- sapply(
              colnames(data_t1),
              function(x) sub("_.*", "", x)
            )
            colnames(data_t2) <- sapply(
              colnames(data_t2),
              function(x) sub("_.*", "", x)
            )
            network$data <- rbind(data_t1, data_t2)
            rownames(network$graph) <- network$rowNames
            colnames(network$graph) <- network$colNames
            return(network)
        },
        boot_fn = function(data, indices, current_timestep, next_timestep) {
            boot_data <- data[indices, ]
            timestep_t0 <- private$get_subset_by_timestep(boot_data, current_timestep)
            timestep_t1 <- private$get_subset_by_timestep(boot_data, next_timestep)

            row_names <- colnames(timestep_t0)
            col_names <- colnames(timestep_t1)

            timestep_t0 <- as.matrix(timestep_t0)
            timestep_t1 <- as.matrix(timestep_t1)

            coefficients_matrix <- apply(timestep_t1, 2, function(column_vector) {
                fit <- cv.glmnet(timestep_t0, column_vector,
                    family = "gaussian",
                    type.measure = "mse", nfolds = 20
                )
                lambda <- fit$lambda.min
                as.numeric(coef(fit, s = lambda)[-1])
            })

            return(as.vector(coefficients_matrix))
        },
        helper_find_mode = function(x) {
            ux <- unique(x)
            ux[which.max(tabulate(match(x, ux)))]
        },
        helper_add_to_indices = function(vector, min = 0, max = 1, by = 0.5) {
            increment <- matrix(seq(min, max, by = by), ncol = 1)
            increment <- c(increment, increment[2:3])
            for (i in 1:length(vector)) {
                vector[i] <- vector[i] + increment[(i - 1) %% length(increment) + 1]
            }
            return(vector)
        },
        helper_plot_network = function(network, image_name,
                                       edge_threshold = 0.3, node_size = 12, bold = TRUE,
                                       edge_width = 0.2, arrow_width = 2) {
            #' Helper function to plot a network graph using the igraph package.
            #' Network is based on a bipartite graph layout with some custom adjustments to make the graph look better.
            #' The edge width is based on the highest edge weight in the tensor.
            #' Network is exported to the specified save location.
            #'
            #' @param network An object containing the adjacency matrix of the network graph.
            #' @param image_name A character string specifying the name of the output image file.
            #' @param edge_threshold A numeric value specifying the threshold for edge weights. Edges with weights below this threshold will be removed from the plot. Default is 0.3.
            #' @param node_size A numeric value specifying the size of the nodes in the plot. Default is 12.
            #' @param bold A logical value indicating whether node labels should be in bold font. Default is TRUE.
            #' @param edge_width A numeric value specifying the width of edges in the plot. Default is 0.2.
            #' @param arrow_width A numeric value specifying the width of arrowheads on directed edges. Default is 2.

            traits <- c(rep("Depression", 21), rep("CVD", 11))
            max_edge <- private$get_min_and_max_edge(self$tensor)[2]

            png(
                filename = paste0(self$save_location, image_name),
                width = self$image_width, height = self$image_height
            )

            # Create igraph object from adjacency matrix
            igraph <- graph_from_adjacency_matrix(network$graph, weighted = TRUE, mode = "directed")
            V(igraph)$type <- c(rep(TRUE, 21), rep(FALSE, 11))
            V(igraph)$name <- network$labels
            igraph <- delete_edges(igraph, E(igraph)[abs(weight) < edge_threshold])

            # Create bipartite layout
            layout <- layout_as_bipartite(igraph, vgap = 0.4)
            layout <- layout[, c(2, 1)]
            max_size <- 50
            layout <- cbind(layout, c(rep(1, 21), rep(0, 11)))
            layout[layout[, 3] == 1, 2] <- seq(0, max_size, length.out = sum(layout[, 3] == 1))
            layout[layout[, 3] == 0, 2] <- seq(0, max_size, length.out = sum(layout[, 3] == 0))
            layout[, 1] <- private$helper_add_to_indices(layout[, 1], max = 0.3, by = 0.1)

            # Custom adjustments to make the graph look better
            layout[22, 1] <- layout[22, 1] - 0.3
            layout[2, 2] <- layout[2, 2] - 7

            if (bold) {
                font <- 2
            } else {
                font <- 1
            }

            plot(igraph,
                layout = layout,
                vertex.color = c(rep("royalblue", 21), rep("orange", 11)),
                vertex.size = node_size,
                vertex.label.family = "sans",
                vertex.label = V(igraph)$name,
                vertex.label.font = font,
                edge.color = ifelse(E(igraph)$weight > 0, "blue", "red"),
                edge.width = ((abs(E(igraph)$weight) / max_edge) * edge_width),
                edge.arrow.width = arrow_width,
                asp = 0.6
            )

            legend_labels <- c("Depression", "CVD")
            legend_colors <- c("royalblue", "orange")
            legend("bottomright",
                legend = legend_labels,
                col = legend_colors, pch = 22, pt.bg = legend_colors,
                pt.cex = 2.5, cex = 1.5
            )

            title(main = paste0(network$years[1], " → ", network$years[2], network$edge_type), font.main = 2)
            dev.off()
        },
        helper_plot_centrality = function(network_list, legend_name, decreasing = TRUE) {
            df <- data.frame()
            for (i in seq_along(network_list)) {
                network <- network_list[[i]]
                centrality_df <- centralityTable(network)
                centrality_df$graph <- paste0("graph ", i)
                centrality_df$name <- network$name
                df <- rbind(df, centrality_df)
            }
            df$type <- "type 1"


            plot <- df %>%
                mutate(
                    graph = case_when(graph %in% unique(df$graph) ~ df$name),
                    graph = as.factor(graph),
                    node = as.factor(node)
                ) %>%
                mutate(node = fct_reorder(node, value, .desc = decreasing)) %>%
                ggplot(aes(x = node, y = value, group = graph, color = graph)) +
                geom_line(aes(linetype = I("solid")), size = 1) +
                labs(x = "", y = "") +
                scale_color_manual(name = legend_name, values = c("royalblue", "orange", "forestgreen", "firebrick")) +
                scale_linetype_manual(name = legend_name, values = rep("solid", length(network_list))) +
                guides(linetype = "none") +
                coord_flip() +
                facet_grid(~measure) +
                theme_bw()

            image_name <- paste0(self$save_location, "centrality_plot", network$edge_type, ".png")
            ggsave(image_name, plot = plot, width = 20, height = 10)
        },
        helper_plot_boot_weights = function(network, image_name, n_cores) {
            png(
                filename = paste0(self$save_location, image_name),
                width = self$image_width, height = self$image_height
            )
            bootnet <- bootnet(network,
                default = "EBICglasso",
                nBoots = self$n_boots, nCores = n_cores
            )
            plot(bootnet, labels = FALSE, order = "sample")
            title(main = paste0(network$years[1], " → ", network$years[2], network$edge_type), font.main = 2)
            dev.off()
        },
        helper_plot_boot_centralities = function(network, image_name, n_cores) {
            png(
                filename = paste0(self$save_location, image_name),
                width = self$image_width, height = self$image_height
            )
            bootnet_case_dropping <- bootnet(network,
                nBoots = self$n_boots,
                type = "case",
                nCores = n_cores,
                statistics = c(
                    "outStrength",
                    "inStrength",
                    "closeness",
                    "betweenness",
                    "outExpectedInfluence",
                    "inExpectedInfluence"
                )
            )

            plot(bootnet_case_dropping, "all")
            title(main = paste0(network$years[1], " → ", network$years[2], network$edge_type), font.main = 2)
            dev.off()
        }
    )
)

library(readr)
data <- read_csv("D:\\Programming\\uva\\2023-2024\\scriptie\\data\\dataset.csv")

coef <- obj$build_coef()

obj <- NetworkAnalysis$new(
    dataset = data[, -1], n_boots = 10,
    save_location = "D:\\Programming\\uva\\2023-2024\\scriptie\\images\\"
)
obj$set_coef(coef)
network <- obj$plot_networks()
obj$plot_boot_weights()
obj$plot_functions()

# plot_functions = function() {
#   self$build_coef()
#   self$plot_networks()
#   self$plot_centralities()
#   self$plot_boot_weights()
#   self$plot_boot_centralities()
# },
