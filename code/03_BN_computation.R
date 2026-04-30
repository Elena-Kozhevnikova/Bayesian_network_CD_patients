# ==============================================================================
# Script 03: Bayesian Network (BN) Construction - TWO NETWORKS
# Project: Bayes_26
# ==============================================================================

suppressWarnings({
  suppressPackageStartupMessages({
    library(bnlearn)
    library(tidyverse)
    library(here)
    library(parallel)
  })
})

# ==============================================================================
# SETTINGS (CHOOSE VARIANT)
# ==============================================================================
# current_hubs_file <- "hubs_CLASSIC_original.csv"
current_hubs_file <- "hubs_5_modules_20_genes.csv"
# current_hubs_file <- "hubs_10_modules_10_genes.csv"
# current_hubs_file <- "hubs_20_modules_5_genes.csv"

# R_BOOTSTRAP <- 100
R_BOOTSTRAP <- 10
prefix <- gsub(".csv", "", current_hubs_file)

# ==============================================================================
cat("\n--- Running Script 03: Bayesian Networks ---\n")
cat("Using hub genes file:", current_hubs_file, "\n")

# 1. Data Loading
datExpr <- readRDS(here("results", "02_datExpr.rds"))
metadata <- readRDS(here("results", "cleaned_metadata.rds")) %>% rownames_to_column("sample_id")
hub_genes <- read_csv(here("results", current_hubs_file), show_col_types = FALSE)

# DATASET 1: Gene expression only
bn_data_genes <- as.data.frame(datExpr[, colnames(datExpr) %in% hub_genes$gene])

# Prepare Calprotectin data (handling missing values with group medians)
group_medians <- metadata %>%
  mutate(Calprotectin_ug_g = as.numeric(ifelse(Calprotectin_ug_g %in% c("NA", "na"), NA, Calprotectin_ug_g))) %>%
  group_by(Patient_group) %>%
  summarise(Cal_med = median(Calprotectin_ug_g, na.rm = TRUE), .groups = 'drop')

meta_clean <- metadata %>%
  mutate(Calprotectin_ug_g = as.numeric(ifelse(Calprotectin_ug_g %in% c("NA", "na"), NA, Calprotectin_ug_g))) %>%
  left_join(group_medians, by = "Patient_group") %>%
  mutate(Calprotectin_ug_g = ifelse(is.na(Calprotectin_ug_g), Cal_med, Calprotectin_ug_g))

# Synchronize samples
shared_samples <- intersect(rownames(bn_data_genes), meta_clean$sample_id)
bn_data_genes <- bn_data_genes[shared_samples, ]
meta_subset <- meta_clean[match(rownames(bn_data_genes), meta_clean$sample_id), ]

# DATASET 2: Genes + Calprotectin phenotype
bn_data_cal <- cbind(bn_data_genes, Calprotectin_ug_g = meta_subset$Calprotectin_ug_g) %>%
  dplyr::select(-any_of("X"))

write_csv(bn_data_cal %>% rownames_to_column("sample_id"), here("results", paste0("BN_training_data_", prefix, ".csv")))

# Configure parallel processing cluster
n_cores <- max(1, detectCores() - 1)
cl <- makeCluster(n_cores)

# ==============================================================================
# NETWORK 1: GENE-ONLY NETWORK (NO PHENOTYPE)
# ==============================================================================
cat("\n[1/2] Computing Gene-Only Network...\n")
clusterExport(cl, c("bn_data_genes"))

boot_genes <- boot.strength(data = bn_data_genes, R = R_BOOTSTRAP, algorithm = "hc", 
                            algorithm.args = list(score = "bic-g"), cluster = cl)

sig_edges_genes <- boot_genes[boot_genes$strength > 0.5 & boot_genes$direction > 0.5, ]
cons_bn_genes <- averaged.network(sig_edges_genes, threshold = 0.5)

if (!directed(cons_bn_genes)) cons_bn_genes <- cextend(cons_bn_genes)

fit_bn_genes <- bn.fit(cons_bn_genes, data = bn_data_genes)

saveRDS(cons_bn_genes, here("results", paste0("BN_GENES_ONLY_structure_", prefix, ".rds")))
saveRDS(fit_bn_genes, here("results", paste0("BN_GENES_ONLY_fitted_", prefix, ".rds")))
saveRDS(boot_genes, here("results", paste0("BN_GENES_ONLY_boot_", prefix, ".rds")))

# ==============================================================================
# NETWORK 2: CALPROTECTIN-AUGMENTED NETWORK
# ==============================================================================
cat("\n[2/2] Computing Calprotectin-Augmented Network...\n")
# Define blacklist: Calprotectin cannot be a parent to gene nodes
blacklist <- data.frame(from = "Calprotectin_ug_g", to = setdiff(colnames(bn_data_cal), "Calprotectin_ug_g"))
clusterExport(cl, c("bn_data_cal", "blacklist"))

boot_cal <- boot.strength(data = bn_data_cal, R = R_BOOTSTRAP, algorithm = "hc", 
                          algorithm.args = list(blacklist = blacklist, score = "bic-g"), cluster = cl)

sig_edges_cal <- boot_cal[boot_cal$strength > 0.5 & boot_cal$direction > 0.5, ]
cons_bn_cal <- averaged.network(sig_edges_cal, threshold = 0.5)

if (!directed(cons_bn_cal)) cons_bn_cal <- cextend(cons_bn_cal)

fit_bn_cal <- bn.fit(cons_bn_cal, data = bn_data_cal)

saveRDS(cons_bn_cal, here("results", paste0("BN_CALPROTECTIN_structure_", prefix, ".rds")))
saveRDS(fit_bn_cal, here("results", paste0("BN_CALPROTECTIN_fitted_", prefix, ".rds")))
saveRDS(boot_cal, here("results", paste0("BN_CALPROTECTIN_boot_", prefix, ".rds")))

stopCluster(cl)

# ==============================================================================
# RESULTS SUMMARY
# ==============================================================================
cat("\n--- BN CONSTRUCTION SUMMARY ---\n")
cat("1. Gene-Only Network: nodes -", length(nodes(cons_bn_genes)), "| arcs -", length(arcs(cons_bn_genes))/2, "\n")
cat("2. Calprotectin-Augmented Network: nodes -", length(nodes(cons_bn_cal)), "| arcs -", length(arcs(cons_bn_cal))/2, "\n")

parents <- cons_bn_cal$nodes$Calprotectin_ug_g$parents
cat("\nParent nodes (predictors) for Calprotectin:\n")
if (length(parents) > 0) {
  print(parents)
} else {
  cat("No significant predictors found for Calprotectin.\n")
}

cat("\nScript 03 completed successfully!\n")