#' NetworkAnalysis Class
#'
#' The class provides a framework for performing network analysis on datasets.
#' It includes functionalities for initializing the analysis with a dataset,
#' setting coefficients tensor, and managing required packages for the analysis.
#' Additional functionalities such as bootstrapping and network visualization
#' can be added to enhance its capabilities.
#'
#' Made by Tycho Stam
#' https://github.com/kingilsildor
#' Assumptions about the data are noted at the end of the script.


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
        #' Function to set the coefficients tensor for the NetworkAnalysis class.
        #' Used to set the coefficients tensor after it has been calculated for debuging uses.
        #' @param coef A tensor containing the coefficients for each timestep.
        #' @examples
        #' set_coef(coef)
        set_coef = function(coef) {
            self$tensor <- coef
        },
        #' Function to get the required packages for the NetworkAnalysis class.
        #' @examples
        #' get_imports()
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
        #' Function to get the basic descriptive statistics of the dataset.
        #' The function calculates the number of observations, quartiles, minimum, maximum, mean, median, mode, and standard deviation.
        #' The function also removes missing values from the dataset before calculating the statistics.
        #' @param round_decimal An integer specifying the number of decimal places to round the results to. Default is 2.
        #' @return A data frame containing the basic descriptive statistics of the dataset.
        #' @examples
        #' get_basic_descriptive()
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
        #' Function to build the tensor of coefficients for each timestep in the dataset.
        #' The function calculates the coefficients for each timestep using the `calc_coef` function.
        #' Time is depended on the number of boots selected when creating the object.
        #' @return A tensor containing the coefficients for each timestep.
        #' @examples
        #' build_coef()
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
        #' Interface function to plot all the different functions in the NetworkAnalysis class.
        #' Also calls the calculate coefficient function.
        #' @examples
        #' plot_functions()
        plot_functions = function() {
            self$build_coef()
            self$plot_networks()
            self$plot_centralities()
            self$plot_boot_weights()
        },
        #' Interface function to plot the descriptive statistics of the dataset.
        #' Saves the table as an image in the specified save location.
        #' @examples
        #' plot_descriptive()
        plot_descriptive = function() {
            df <- get_basic_descriptive()
            table <- tableGrob(df)
            ggsave(paste0(self$save_location, "descriptive.png"),
                plot = table, width = self$image_width, height = self$image_height, dpi = 300
            )
        },
        #' Interface function to plot the network graph for each timestep in the tensor.
        #' Makes a distinction between the full network and the reduced network.
        #' @param legend A logical value indicating whether a legend should be included in the plot. Default is TRUE.
        #' @return A network object, mostly to use in manual testing.
        #' @examples
        #' plot_networks()
        plot_networks = function(legend = TRUE) {
            for (i in 1:dim(self$tensor)[3]) {
                graph <- self$tensor[, , i]
                network <- private$calc_network(graph, i)
                for (group_edge in c(TRUE, FALSE)) {
                    if (group_edge) {
                        # Remove the suffixes from the column names
                        colnames(network$graph) <- sapply(
                            network$colNames,
                            function(x) sub("_.*", "", x)
                        )
                        network$edge_type <- " Network "
                        image_name <- paste0("network_", network$years[1], network$years[2], ".png")
                    } else {
                        network <- private$remove_group_edge(network)
                        network$edge_type <- " Network reduced "
                        image_name <- paste0("network_", network$years[1], network$years[2], "_reduced.png")
                    }
                    private$helper_plot_network(network, image_name, legend, edge_width = 8)
                }
            }
            return(network)
        },
        #' Interface function to plot the different centrality metrics for each network in the tensor.
        #' Current centrality metrics used are: betweenness, closeness, in & out strenght,and in & out expected influence.
        #' @param decreasing A logical value indicating whether the centrality metrics should be sorted in decreasing order. Default is TRUE.
        #' @param legend_name A character string specifying the name of the legend in the plot. Default is "Centrality metric for each timestep".
        #' @examples
        #' plot_centralities()
        plot_centralities = function(decreasing = TRUE,
                                     legend_name = "Centrality metric for each timestep") {
            network_list <- list()
            for (group_edge in c(TRUE, FALSE)) {
                for (i in 1:dim(self$tensor)[3]) {
                    graph <- self$tensor[, , i]
                    network <- private$calc_network(graph, i)

                    if (group_edge) {
                        # Remove the suffixes from the column names
                        colnames(network$graph) <- sapply(
                            network$colNames,
                            function(x) sub("_.*", "", x)
                        )
                        indicator <- 2
                        network$edge_type <- ""
                    } else {
                        indicator <- 0
                        network <- private$remove_group_edge(network)
                        network$edge_type <- "reduced"
                    }
                    network$name <- paste(network$year[1], "→", network$year[2], network$edge_type)
                    network_list[[i + indicator]] <- network
                }
            }
            private$helper_plot_centrality(network_list, legend_name)
        },
        #' Interface function to plot the bootstrapped weights for each network in the tensor.
        #' Helper function to plot the bootstrapped weight is called for each network in the tensor.
        #' @param n_cores The number of cores to use for parallel processing. Default is 8.
        #' @examples
        #' plot_boot_weights()
        plot_boot_weights = function(n_cores = 8) {
            for (group_edge in c(TRUE, FALSE)) {
                for (i in 1:dim(self$tensor)[3]) {
                    graph <- self$tensor[, , i]
                    network <- private$calc_network(graph, i)

                    # Plot the bootstrapped weights for the network
                    if (group_edge) {
                        # Remove the suffixes from the column names
                        colnames(network$graph) <- sapply(
                            network$colNames,
                            function(x) sub("_.*", "", x)
                        )
                        network$edge_type <- " Network "
                        image_name <- paste0("bootWeights_", network$years[1], network$years[2], ".png")
                    } else {
                        network <- private$remove_group_edge(network)
                        network$edge_type <- " Network reduced "
                        image_name <- paste0("bootWeights_", network$years[1], network$years[2], "_reduced.png")
                    }

                    private$helper_plot_boot_weights(network, image_name, n_cores)
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
            "Acetate", "Apoprotein", "C-reactive Protein", "Diastolic Blood Pressure", "Glucose", "Cholesterol HDL",
            "Insulin", "Cholesterol LDL", "Systolic Blood Pressure", "Cholesterol Total", "Triglycerides"
        ),
        var_start_time = NULL,
        var_stop_time = NULL,
        #' Function to start the timer.
        #' @examples
        #' start_time()
        start_time = function() {
            private$var_start_time <- Sys.time()
        },
        #' Function to stop the timer and print the time spent.
        #' @examples
        #' stop_time()
        stop_time = function() {
            private$var_stop_time <- Sys.time()
            cat("time spent: ", private$var_stop_time - private$var_start_time, "\n")
        },
        #' Function to get all the timesteps from a given dataset.
        #' @return A numeric vector of sorted unique timesteps.
        #' @examples
        #' get_timesteps_list()
        get_timesteps_list = function() {
            column_names <- colnames(self$dataset)
            timesteps_list <- unique(as.numeric(
                regmatches(column_names, regexpr("\\d{4}", column_names))
            ))

            timesteps_list <- sort(as.numeric(timesteps_list))
            return(timesteps_list)
        },
        #' Function to create a subsets of the input data based on the specified timestep.
        #' @param data A data frame from which to subset.
        #' @param timestep A string indicating the timestep to use for subsetting.
        #' @return A data frame containing only the columns from the input data that contain the specified timestep in their names.
        #' @examples
        #' get_subset_by_timestep(data, timestep)
        get_subset_by_timestep = function(data, timestep) {
            subset_data <- data[, grep(paste0(timestep), names(data), value = TRUE)]
            return(subset_data)
        },
        #' Function to calculate the min and max values in the tensor
        #' @param tensor A tensor (a multi-dimensional array) of numeric values.
        #' @return A vector of length two, where the first element is the minimum value in the tensor and the second element is the maximum value.
        get_min_and_max_edge = function(tensor) {
            return(c(min(tensor), max(tensor)))
        },
        get_suffix = function(name) {
            str_extract(name, "_[^_]+$")
        },
        #' Function to remove the edges between nodes in the same group from the network graph.
        #' Extra option to also remove self-edges.
        #' @param network A list containing the network information, including the graph.
        #' @return The input network with group-edges removed from the graph.
        #' @examples
        #' remove_group_edge(network)
        remove_group_edge = function(network) {
            colnames(network$graph) <- network$colNames
            row_suffixes <- sapply(rownames(network$graph), private$get_suffix)
            col_suffixes <- sapply(colnames(network$graph), private$get_suffix)

            for (i in 1:nrow(network$graph)) {
                for (j in 1:ncol(network$graph)) {
                    if (row_suffixes[i] == col_suffixes[j] && i != j) {
                        network$graph[i, j] <- 0
                    }
                }
            }
            # Remove the suffixes from the column names
            colnames(network$graph) <- sapply(
                network$colNames,
                function(x) sub("_.*", "", x)
            )
            return(network)
        },
        #' Function to calculates the coefficients for a given timestep using bootstrap sampling.
        #' The bootstrap sampling is performed on the dataset using a custom statistic function that calls the `boot_fn` function.
        #' The mean of the bootstrap results is then calculated for each column to obtain the bootstrapped coefficients matrix.
        #' @param timestep An index indicating the current timestep.
        #' @return A matrix of bootstrapped coefficients for the current timestep.
        #' @examples
        #' calc_coef(timestep)
        calc_coef = function(timestep) {
            timesteps_list <- private$get_timesteps_list()
            current_timestep <- timesteps_list[timestep]
            next_timestep <- timesteps_list[timestep + 1]

            # Create a progress bar
            pb <- progress_bar$new(
                format = paste0(
                    "Timestep: ", current_timestep, "-", next_timestep,
                    " Bootstrapping [:bar] :percent in :elapsed"
                ),
                total = self$n_boots,
                clear = FALSE,
                width = 60
            )

            # Calculate the bootstrapped coefficients for each column of the next timestep
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
        #' Function that creates a network based on the input graph and a specified timestep.
        #' The function also modifies the column names of the data subsets for the two timesteps.
        #' @param graph A graph object on which the network will be estimated.
        #' @param i An index indicating the current timestep. The function will use this and the next timestep to subset the data.
        #' @return A network object containing the adjacency matrix, to use in future functions.
        #' @examples
        #' calc_network(graph, i)
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

            # Create empty network object and fill it with the calculated values
            network <- estimateNetwork(graph, default = "IsingFit")
            network$labels <- node_adverbiations
            network$graph <- graph
            network$rowNames <- colnames(data_t1)
            network$colNames <- colnames(data_t2)
            network$years <- c(year_t1, year_t2)

            # Remove the suffixes from the column names
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
        #' Function to performs bootstrap sampling on the input data
        #' and calculates the coefficients of a generalized linear model for each column of the next timestep using the glmnet package.
        #' The model is trained using cross-validation to find the optimal lambda value that minimizes the mean squared error.
        #' @param data A data frame containing the data to be bootstrapped.
        #' @param indices A vector of indices specifying the bootstrap sample from the data.
        #' @param current_timestep A value indicating the current timestep.
        #' @param next_timestep A value indicating the next timestep.
        #' @return A matrix with all the glmnet coefficients.
        #' @examples
        #' boot_fn(data, indices, current_timestep, next_timestep)
        boot_fn = function(data, indices, current_timestep, next_timestep) {
            boot_data <- data[indices, ]
            timestep_t0 <- private$get_subset_by_timestep(boot_data, current_timestep)
            timestep_t1 <- private$get_subset_by_timestep(boot_data, next_timestep)

            row_names <- colnames(timestep_t0)
            col_names <- colnames(timestep_t1)

            timestep_t0 <- as.matrix(timestep_t0)
            timestep_t1 <- as.matrix(timestep_t1)

            # Get the coefficients for each column of the next timestep
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
        #' Helper function to find the mode in a given vector. If there are multiple modes, it returns the first one.
        #' @param x A numeric or character vector from which the mode is to be found.
        #' @return The mode of the vector.
        #' @examples
        #' helper_find_mode(c(1, 2, 2, 3, 3, 3, 4))
        helper_find_mode = function(x) {
            ux <- unique(x)
            ux[which.max(tabulate(match(x, ux)))]
        },
        #' Helper function to add a sequence of increments to each element of a vector.
        #' The increments help with the visability of the network graph.
        #' @param vector A numeric vector to which the increments will be added.
        #' @param min The starting value of the sequence of increments. Default is 0.
        #' @param max The ending value of the sequence of increments. Default is 1.
        #' @param by The increment of the sequence. Default is 0.5.'
        #' @return A numeric vector with the increments added to the original elements.
        #' @examples
        #' helper_add_to_indices(c(1, 2, 3))
        helper_add_to_indices = function(vector, min = 0, max = 1, by = 0.5) {
            increment <- matrix(seq(min, max, by = by), ncol = 1)
            increment <- c(increment, increment[2:3])
            for (i in 1:length(vector)) {
                vector[i] <- vector[i] + increment[(i - 1) %% length(increment) + 1]
            }
            return(vector)
        },
        #' Helper function to plot a network graph using the igraph package.
        #' Network is based on a bipartite graph layout with some custom adjustments to make the graph look better.
        #' The edge width is based on the highest edge weight in the tensor.
        #' Network is exported to the specified save location.
        #' @param network An object containing the adjacency matrix of the network graph.
        #' @param image_name A character string specifying the name of the output image file.
        #' @param edge_threshold A numeric value specifying the threshold for edge weights.
        #' Edges with weights below this threshold will be removed from the plot. Default is 0.3.
        #' @param node_size A numeric value specifying the size of the nodes in the plot. Default is 12.
        #' @param bold A logical value indicating whether node labels should be in bold font. Default is TRUE.
        #' @param edge_width A numeric value specifying the width of edges in the plot. Default is 0.2.
        #' @param arrow_width A numeric value specifying the width of arrowheads on directed edges. Default is 2.
        #' @examples
        #' helper_plot_network(network, "network.png")
        helper_plot_network = function(network, image_name, legend,
                                       edge_threshold = 0.3, node_size = 12, bold = TRUE,
                                       edge_width = 0.2, arrow_width = 2) {
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
                vertex.color = c(rep("cornflowerblue", 21), rep("orange", 11)),
                vertex.size = node_size,
                vertex.label.family = "sans",
                vertex.label = V(igraph)$name,
                vertex.label.font = font,
                vertex.label.cex = 1.8,
                edge.color = ifelse(E(igraph)$weight > 0, "blue", "red"),
                edge.width = ((abs(E(igraph)$weight) / max_edge) * edge_width),
                edge.arrow.width = arrow_width,
                asp = 0.6
            )

            if (legend) {
                legend_labels <- c("Depression", "CVD")
                legend_colors <- c("cornflowerblue", "orange")
                legend("bottomright",
                    legend = legend_labels,
                    col = legend_colors, pch = 22, pt.bg = legend_colors,
                    pt.cex = 2.5, cex = 1.5
                )
                title(main = paste0(network$years[1], " → ", network$years[2], network$edge_type), font.main = 2)
            }

            dev.off()
        },
        #' Helper function to plot the centrality metrics for each network in the tensor.
        #' Plot all the centrality metrics for each network in the tensor in a single plot.
        #' @param network_list A list containing the network objects for each timestep.
        #' @param legend_name A character string specifying the name of the legend in the plot.
        #' @param decreasing A logical value indicating whether the centrality metrics should be sorted in decreasing order. Default is TRUE.
        #' @examples
        #' helper_plot_centrality(network_list, "Centrality metric for each timestep")
        helper_plot_centrality = function(network_list, legend_name, decreasing = TRUE) {
            df <- data.frame()
            for (i in seq_along(network_list)) {
                network <- network_list[[i]]
                colnames(network$graph) <- private$full_name_labels

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
                scale_color_manual(name = legend_name, values = c("cornflowerblue", "orange", "forestgreen", "indianred")) +
                scale_linetype_manual(name = legend_name, values = rep("solid", length(network_list))) +
                guides(linetype = "none") +
                coord_flip() +
                facet_grid(~measure) +
                theme_bw()

            image_name <- paste0(self$save_location, "centrality_plot", network$edge_type, ".png")
            ggsave(image_name, plot = plot, width = 20, height = 10)
        },
        #' Helper function to plot the bootstrapped stability of edge weight for each network in the tensor.
        #' Stability is calculated using the bootnet package.
        #' @param network A network object containing the adjacency matrix of the network graph.
        #' @param image_name A character string specifying the name of the output image file.
        #' @param n_cores The number of cores to use for parallel processing. Default is 8.
        #' @examples
        #' helper_plot_boot_weights(network, "bootWeights.png")
        helper_plot_boot_weights = function(network, image_name, n_cores) {
            bootnet <- bootnet(network,
                default = "EBICglasso",
                nBoots = self$n_boots, nCores = n_cores
            )

            ggsave(
                filename = paste0(self$save_location, image_name),
                width = 10, height = 7, units = "in",
                plot = plot(bootnet, labels = FALSE, order = "sample"),
            )
        }
    )
)

#' The script is build on the assumption that the data is already preprocessed and saved as a CSV file.
#' The output of the preprocessing is a dataframe with floating types expect for the patient ID column.

# Output of the preprocessing.
data <- read_csv("..\\dataset.csv")
# Remove patient ID column
data <- data[, -1]

# Create the network object
obj <- NetworkAnalysis$new(
    dataset = data, n_boots = 10,
    save_location = "..\\images\\"
)

obj$plot_functions()