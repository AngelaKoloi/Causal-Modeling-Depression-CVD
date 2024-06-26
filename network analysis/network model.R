
####################
library(qgraph)
library(mgm)
library(Rmisc)
library('ggplot2')
library(tidygraph)
library(dplyr)
library(ggraph)
library(RColorBrewer)
library(magrittr)
library(tidyverse)



#Import data 
data_resid <- read.csv("residualized_data2007SA.csv")

# Randomly drop 10%, 20%, 30% sample
data_resid <- data_resid %>% 
  sample_frac(0.7, replace = FALSE)


apply(data_resid, 2, function(x) length(unique(x)))

type_list <- c(rep("g", 74))

level_list <- c(rep(1, 74))

# Run the MGM model and do prediction
########## 
set.seed(123)
fit_mgm <- mgm(data = data.matrix(data_resid), 
               type = type_list,
               level = level_list,
               k = 2,
               #k = 3,
              # lambdaSel = "EBIC", 
               lambdaFolds = 10, #comment out when doing EBIC
               
               #lambdaGam = 0., #comment out when doing , = 0 meaning effectively allowing for all possible interactions in the model.
               ruleReg = "OR")

# do model prediction
pred_mgm <- predict(object = fit_mgm,
                    data = data_resid,
                    errorCon = c("RMSE", "R2"), 
                    errorCat = c("CC", "nCC"))

fit_mgm$pairwise$wadj

# Convert adjacency matrix to data frame
adjacency_df <- as.data.frame(fit_mgm$pairwise$wadj)

View(adjacency_df)

rownames(adjacency_df) <- colnames(data_resid)
colnames(adjacency_df) <- colnames(data_resid)




# Create a graph from the adjacency matrix with edge weights
graph <- graph_from_adjacency_matrix(as.matrix(adjacency_df), weighted = TRUE)

# Set the edge weights from the adjacency_df dataframe
#E(graph)$weight <- adjacency_df[which(adjacency_df != 0)]

# Export the graph to a GraphML file
#write_graph(graph, "graph2007.graphml", format = "graphml")

# Export adjacency matrix with node names
write.table(adjacency_df, "adjacency_matrix_2007SA.csv", sep = ",", col.names = FALSE, row.names = FALSE)




adjacency_df[1:49, 1:49] <- 0.0 #metabolites
adjacency_df[50:70, 50:70] <- 0.0 #depressive symptoms 
adjacency_df[71:74, 71:74] <- 0.0 #cvd risk factors

adjacency_df

library(igraph)

# Create a graph from the adjacency matrix with edge weights
graph <- graph_from_adjacency_matrix(as.matrix(adjacency_df), weighted = TRUE)

# Set the edge weights from the adjacency_df dataframe
#E(graph)$weight <- adjacency_df[which(adjacency_df != 0)]

# Export the graph to a GraphML file
write_graph(graph, "graph2007intra.graphml", format = "graphml")

# Export adjacency matrix with node names
write.table(adjacency_df, "adjacency_matrix_2007intra.csv", sep = ",", col.names = FALSE, row.names = FALSE)



#BOOTSTRAPPING
set.seed(123)
res_obj <- resample(object = fit_mgm,
                    data = data.matrix(data_resid),
                    nB = 100,
                   quantiles = c(.05,.95)
                    )



# Convert bootParameters to a data frame
boot_params_df <- as.data.frame(res_obj$bootParameters)

# Export adjacency matrix with node names
write.table(boot_params_df, "boot_params_df.csv", sep = ",", col.names = TRUE, row.names = TRUE)

#stability 

library(ggplot2)

# Set the default font to "sans" for Windows
if (Sys.info()["sysname"] == "Windows") {
  windowsFonts(A = windowsFont("sans"))
}

#stable connections

rownrs <-c( 4,6,7,38,39,40,40) #adjust
colnrs <- c(61,54,65,67,69,52,65) #adjust

rownrs <-c(4,6,7,7,38,39,39,40)
colnrs <- c(74,74,73,74,73,72,74,74)

#make lists to select important links
link_list_left <- c()
link_list_right <- c()
link_list <- c()
ci_upper <- c()
ci_mean <- c()
ci_lower <- c()
sd <- c()
zeros <- c()


# Replace these lines to map indices to node labels
node_labels <- c('Alanine',
                 'Glutamine',
                 'Histidine',
                 'Isoleucine',
                 'Leucine',
                 'Phenylalanine',
                 'Tyrosine',
                 'Valine',
                 'Apolipoprotein A-I',
                 'Apolipoprotein B',
                 'Esterified cholesterol',
                 'Free cholesterol',
                 'Total cholesterol HDL',
                 'Total cholesterol HDL2',
                 'Total cholesterol HDL3',
                 'Total cholesterol LDL',
                 'Remnant cholesterol',
                 'Serum total cholesterol',
                 'Total cholesterol VLDL',
                 'Linoleic acid',
                 'Docosahexaenoic acid',
                 'Linoleic acid R18:2',
                 'Monounsaturated fatty acids',
                 'Omega-3 fatty acids',
                 'Omega-6 fatty acids',
                 'Polyunsaturated fatty acids',
                 'Saturated fatty acids',
                 'Albumin signal area',
                 'Creatinine',
                 'Phosphatidylcholine',
                 'Serum total triglycerides',
                 'Sphingomyelins',
                 'Total cholines',
                 'Total phosphoglycerides',
                 'Triglycerides in HDL',
                 'Triglycerides in LDL',
                 'Triglycerides in VLDL',
                 'Citrate',
                 'Glucose',
                 'Lactate',
                 'Glycoprotein acetyls',
                 '3-Hydroxybutyrate',
                 'Acetate',
                 'Acetoacetate',
                 'HDL particles diameter',
                 'LDL particles diameter',
                 'VLDL particles diameter',
                 'Degree of unsaturation',
                 'Total fatty acids',
                 
                 'Sadness',
                 'Pessimism',
                 'Past failure',
                 'Loss of pleasure',
                 'Guilty feelings',
                 'Punishment feelings',
                 'Self-dislike',
                 'Self-criticalness',
                 'Suicidal thought or wishes',
                 'Crying',
                 'Agitation',
                 'Loss of interest',
                 'Indecisiveness',
                 'Worthlessness',
                 'Loss of energy',
                 'Changes in sleep pattern',
                 'Irritability',
                 'Changes in appetite',
                 'Concentration difficulty',
                 'Tiredness or fatigue',
                 'Loss of interest in sex',
                 'Hypertension', 'Systolic blood pressure',
                 'Diastolic blood pressure',
                 'Intima-media thickness')  


for (i in 1:length(rownrs)){
  link_list[i] <- paste( node_labels[rownrs[i]], '\n', node_labels[colnrs[i]])
  ci_upper[i] <- CI(res_obj$bootParameters[rownrs[i],colnrs[i],1:100], ci=0.95)['upper']
  ci_mean[i] <- CI(res_obj$bootParameters[rownrs[i],colnrs[i],1:100], ci=0.95)['mean']
  ci_lower[i] <- CI(res_obj$bootParameters[rownrs[i],colnrs[i],1:100], ci=0.95)['lower']
  sd[i] <- sd(res_obj$bootParameters[rownrs[i],colnrs[i],1:100])
  zeros[i] <- round(1- (sum(res_obj$bootParameters[rownrs[i],colnrs[i],1:100]==0)/100),digits = 2)
}

stability <- data.frame(
  nodes = link_list,
  lower = ci_lower,
  mean = ci_mean,
  upper = ci_upper,
  sd = sd,
  zeros = zeros
)

#stability['sd'] <- stability['mean'] - stability['lower']



# Create the plot
p <- ggplot(stability, aes(x=reorder(nodes, +mean), y=mean)) +
  geom_point(shape=4, alpha=1) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.3) +
  geom_text(label=zeros, nudge_x = 0.35, nudge_y = 0.0) +
  coord_flip() +
  labs(title="Average link strength and 95% CI", y='Link strength',
       x="Node pairs") +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size=9),
        axis.text = element_text(size=7))


# Print and display the plot
print(p)

