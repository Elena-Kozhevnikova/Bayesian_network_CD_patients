# ==============================================================================
# Script 05: GO Enrichment Analysis (Biological Interpretation)
# Project: Bayes_26
# ==============================================================================

suppressWarnings({
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(tidyverse)
    library(here)
  })
})

# Settings
current_prefix <- "hubs_CLASSIC_original"
network_type <- "GENES_ONLY" 

cat("--- Running Script 05: GO Analysis ---\n")

# 1. Data Loading
hub_genes_info <- read_csv(here("results", paste0(current_prefix, ".csv")), show_col_types = FALSE)
boot_res <- readRDS(here("results", paste0("BN_", network_type, "_boot_", current_prefix, ".rds")))

# ==============================================================================
# PART A: Genes from the Strict Bayesian Network (cutoff 0.7)
# ==============================================================================
cat("Extracting genes from the network with strength threshold 0.7...\n")
strong_edges <- boot_res[boot_res$strength > 0.7 & boot_res$direction > 0.5, ]
bn_strict_genes <- unique(c(strong_edges$from, strong_edges$to))
bn_strict_genes <- setdiff(bn_strict_genes, "Calprotectin_ug_g") # Remove phenotype marker from gene list

cat("Number of genes in the strict network for GO analysis:", length(bn_strict_genes), "\n")

ego_bn <- enrichGO(gene          = bn_strict_genes,
                   keyType       = "SYMBOL",
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC", # Cellular component
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

if (!is.null(ego_bn) && nrow(ego_bn) > 0) {
  pdf(here("plots", paste0("05_GO_BN_strict_0.7_", current_prefix, ".pdf")), width = 10, height = 8)
  print(dotplot(ego_bn, showCategory = 15, title = "GO: Bayesian Network (strength > 0.7)"))
  dev.off()
  write_csv(as.data.frame(ego_bn), here("results", paste0("GO_BN_strict_", current_prefix, ".csv")))
} else {
  cat("No significant GO terms found for the network genes.\n")
}

# ==============================================================================
# PART B: Genes by WGCNA Modules
# ==============================================================================
cat("\nRunning GO analysis by WGCNA modules...\n")
modules <- unique(hub_genes_info$module)

pdf(here("plots", paste0("05_GO_WGCNA_modules_", current_prefix, ".pdf")), width = 12, height = 10)

for (mod in modules) {
  mod_genes <- hub_genes_info$gene[hub_genes_info$module == mod]
  
  ego_mod <- enrichGO(gene = mod_genes, keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
  
  if (!is.null(ego_mod) && nrow(ego_mod) > 0) {
    print(dotplot(ego_mod, showCategory = 10, title = paste("GO Module:", mod)))
    # Save the enrichment table
    write_csv(as.data.frame(ego_mod), here("results", paste0("GO_Module_", mod, "_", current_prefix, ".csv")))
  }
}
dev.off()

cat("Script 05 completed successfully! Plots and tables have been saved.\n")