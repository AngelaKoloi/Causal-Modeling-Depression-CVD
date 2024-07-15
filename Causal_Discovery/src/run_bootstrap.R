library(tpc)
library(micd)
library(parallel)

data <- read.csv("data/discrete_data.csv")

tier1 <- c('Age', 'Sexe')
tier2 <- c('aedu', 'asmokstat', 'AIPMETO2', 'aauditsc', 'aIRSsum9', 'abaiscal', 'aids', 'acidep09', 'amet_syn2', 'ams_waist', 'ams_hpt', 'ams_trig2', 'ams_hdl2', 'ams_gluc2', 'atri_med', 'ahdl_med', 'asbp_med', 'adbp_med', 'agluc_med', 'ahsCRP', 'aIL6', 'aApoB', 'aHDL_C', 'aTotFA', 'aSerum_TG', 'aGp', 'aIle')
tier3 <- c('eage', 'sex')
tier4 <- c('eipmeto2', 'eauditsc', 'eIRSsum9', 'ebaiscal', 'eids', 'ecidep09', 'emet_syn2', 'ems_waist', 'ems_hpt', 'ems_trig2', 'ems_hdl2', 'ems_gluc2', 'etri_med', 'ehdl_med', 'esbp_med', 'edbp_med', 'egluc_med', 'eHSCRP', 'eIL6', 'eApoB', 'eHDLC', 'eTotFA', 'eSerumTG', 'eGp', 'eIle')

tiered_variables <- c(tier1, tier2, tier3, tier4)
data <- data[, tiered_variables]

# Create the tiers vector based on the reordered data
tiers <- c(rep(1, length(tier1)), rep(2, length(tier2)), rep(3, length(tier3)), rep(4, length(tier4)))

# Blacklist
bl <- c("sex", "Age", "eage", "Sexe")

# Prepare forbidden edges
forbEdges <- matrix(FALSE, ncol = length(colnames(data)), nrow = length(colnames(data)), 
                    dimnames = list(colnames(data), colnames(data)))

# Populate the forbidden edges matrix
for (var in bl) {
  if (var %in% colnames(data)) {
    forbEdges[, var] <- TRUE 
  }
}

# Whitelist
wl_tier1 <- rbind(
  data.frame(from = "Age", to = c("aSerum_TG", "aHDL_C", "aApoB")),
  data.frame(from = "Sexe", to = c("aSerum_TG", "aHDL_C", "aApoB")),
  data.frame(from = "Sexe", to = c("asmokstat", "acidep09"))
)

wl_tier2 <- rbind(
  data.frame(from = "asmokstat", to = c("aSerum_TG", "aHDL_C", "aApoB")),
  expand.grid(from = c("aSerum_TG", "aHDL_C", "aApoB"), to = c("atri_med")),
  expand.grid(from = c("aSerum_TG", "aHDL_C", "aApoB"), to = c("ahdl_med")),
  data.frame(from = "ams_gluc2", to = c("agluc_med")),
  expand.grid(from = "ams_hpt", to = c("asbp_med", "adbp_med")),
  data.frame(from = "acidep09", to = c("aedu")),
  data.frame(from = "aauditsc", to = c("asmokstat")),
  expand.grid(from = "asmokstat", to = c("eSerumTG", "eHDLC", "eApoB")),
  expand.grid(from = c("ahsCRP"), to = c("eIL6")),
  data.frame(from = "acidep09", to = c("emet_syn2")),
  data.frame(from = "aedu", to = c("eauditsc"))
)

wl <- rbind(wl_tier1, wl_tier2)
context_tier1 <- unique(wl_tier1$from)
context_tier2 <- unique(wl_tier2$from)

# Create context.tier vector
context_tier <- rep(NA, ncol(data))
names(context_tier) <- colnames(data)
context_tier[names(data) %in% context_tier1] <- names(data)[names(data) %in% context_tier1]
context_tier[names(data) %in% context_tier2] <- names(data)[names(data) %in% context_tier2]
context_tier <- context_tier[!is.na(context_tier)]
print(context_tier)

bootstrap_analysis <- function(data, drop_rate) {
  suff_stat_list <- list()
  for (i in 1:20) {
    print(paste("Starting bootstrap iteration:", i))
    boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
    final_data <- boot_data[sample(nrow(boot_data), nrow(boot_data) * (1 - drop_rate)), ]
    
    suff_stat <- getSuff(final_data, test = "flexCItest")
    graph <- tpc(
      suffStat = suff_stat,
      indepTest = flexCItest,
      skel.method = "stable.parallel",
      labels = as.character(colnames(final_data)),
      alpha = 0.05,
      tiers = tiers,
      forbEdges = forbEdges,
      numCores = detectCores() - 1,
      verbose = FALSE,
      context.tier = context_tier
    )
    suff_stat_list[[i]] <- graph
  }
  return(suff_stat_list)
}

set.seed(123) 

results_0 <- bootstrap_analysis(data, 0.0)
save(results_0, file = "bootstrap_graph.RData")

results_10 <- bootstrap_analysis(data, 0.1)
save(results_10, file = "bootstrap_graph_10.RData")

results_20 <- bootstrap_analysis(data, 0.2)
save(results_20, file = "bootstrap_graph_20.RData")