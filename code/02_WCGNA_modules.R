# ==============================================================================
# Script 02: WGCNA (Weighted Gene Co-expression Network Analysis) and Hub Gene Identification
# Project: Bayes_26
# ==============================================================================

suppressWarnings({
  suppressPackageStartupMessages({
    library(WGCNA)
    library(tidyverse)
    library(here)
  })
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# 1. Data Loading
cat("Loading data...\n")
deg_counts <- readRDS(here("results", "deg_counts_for_wgcna.rds"))
metadata <- readRDS(here("results", "cleaned_metadata.rds"))

datExpr <- t(deg_counts)
# Keep genes expressed in at least 50% of samples
keep_genes <- rowSums(datExpr > 1) >= 0.5 * ncol(datExpr)
datExpr <- datExpr[, keep_genes]

# Prepare Traits for correlation analysis
traitData <- metadata[rownames(datExpr), ] %>%
  mutate(
    CD_status = as.numeric(Patient_group == "Israel_CD"),
    Male = as.numeric(Gender == "male"),
    CRP_mg_L = as.numeric(ifelse(CRP_mg_L == "NA" | CRP_mg_L == "na", NA, CRP_mg_L)),
    Calprotectin_ug_g = as.numeric(ifelse(Calprotectin_ug_g == "NA" | Calprotectin_ug_g == "na", NA, Calprotectin_ug_g))
  )

# Create dummy variables for Age
age_dummies <- model.matrix(~ Age - 1, data = traitData)
colnames(age_dummies) <- gsub("Age", "Age_", colnames(age_dummies))

# Build Trait Matrix and impute missing values with mean
traitMatrix <- cbind(traitData %>% dplyr::select(CD_status, Male, CRP_mg_L, Calprotectin_ug_g) %>%
                       mutate(across(where(is.numeric), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) %>%
                       as.matrix(), age_dummies)
traitMatrix <- apply(traitMatrix, 2, as.numeric)
rownames(traitMatrix) <- rownames(datExpr)


# ==============================================================================
# FUNCTION TO GENERATE DIFFERENT NETWORK VARIANTS
# ==============================================================================
generate_hubs <- function(min_size, cut_height, num_mods, genes_per_mod, prefix, rank_by_calprotectin = TRUE) {
  cat(sprintf("\n--- Building network: %s (minSize=%d) ---\n", prefix, min_size))
  
  net <- blockwiseModules(datExpr, power = 5, networkType = "signed", TOMType = "signed",
                          minModuleSize = min_size, mergeCutHeight = cut_height, deepSplit = 2,
                          pamRespectsDendro = FALSE, verbose = 0)
  
  moduleColors <- net$colors
  MEs <- net$MEs
  adjacency <- adjacency(datExpr, power = 5, type = "signed", corFnc = "bicor")
  kIN <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = TRUE)
  
  if(rank_by_calprotectin) {
    # Rank modules by correlation with Calprotectin
    moduleTraitCor <- cor(MEs, traitMatrix, use = "complete.obs")
    mod_sig <- abs(moduleTraitCor[, "Calprotectin_ug_g"])
    mod_sig <- mod_sig[names(mod_sig) != "MEgrey"]
    
    actual_mods <- min(num_mods, length(mod_sig))
    top_mods_colors <- gsub("ME", "", names(sort(mod_sig, decreasing = TRUE)[1:actual_mods]))
  } else {
    # Select ALL active modules (excluding grey)
    top_mods_colors <- setdiff(unique(moduleColors), "grey")
  }
  
  hub_gene_df <- do.call(rbind, lapply(top_mods_colors, function(mod) {
    mod_genes <- names(moduleColors)[moduleColors == mod]
    mod_connectivity <- kIN[mod_genes, "kWithin"]
    mod_kME <- cor(datExpr[, mod_genes], MEs[, paste0("ME", mod)])
    gs <- as.data.frame(cor(datExpr[, mod_genes], traitMatrix, use = "p"))
    
    # Calculate integrated hub score
    hub_score <- scale(mod_connectivity) + scale(abs(mod_kME)) + scale(rowMeans(abs(gs)))
    top_genes <- mod_genes[order(hub_score, decreasing = TRUE)[1:min(genes_per_mod, length(mod_genes))]]
    
    data.frame(
      gene = top_genes, 
      module = mod, 
      kWithin = kIN[top_genes, "kWithin"],
      kME = cor(datExpr[, top_genes], MEs[, paste0("ME", mod)]),
      GS = rowMeans(abs(cor(datExpr[, top_genes], traitMatrix)))
    )
  }))
  
  # Filter out immunoglobulin and specific transporter genes
  hub_gene_df <- hub_gene_df %>% filter(!grepl("^(IG|IGH|IGK|IGL|SLC)", gene))
  write_csv(hub_gene_df, here("results", paste0(prefix, ".csv")))
  cat("Saved file:", paste0(prefix, ".csv"), "\n")
  
  return(net) # Return network object to save the primary one
}

# ==============================================================================
# EXECUTION OF ALL VARIANTS
# ==============================================================================

# 1. ORIGINAL NETWORK AND 5x20 VARIANT (Standard WGCNA parameters)
net_base <- generate_hubs(min_size = 20, cut_height = 0.20, num_mods = NULL, genes_per_mod = 20, 
                          prefix = "hubs_CLASSIC_original", rank_by_calprotectin = FALSE)

generate_hubs(min_size = 20, cut_height = 0.20, num_mods = 5, genes_per_mod = 20, 
              prefix = "hubs_5_modules_20_genes", rank_by_calprotectin = TRUE)

# 2. 10x10 VARIANT (Higher granularity)
generate_hubs(min_size = 10, cut_height = 0.15, num_mods = 10, genes_per_mod = 10, 
              prefix = "hubs_10_modules_10_genes", rank_by_calprotectin = TRUE)

# 3. 20x5 VARIANT (Small clusters / Micro-modules)
generate_hubs(min_size = 5, cut_height = 0.10, num_mods = 20, genes_per_mod = 5, 
              prefix = "hubs_20_modules_5_genes", rank_by_calprotectin = TRUE)


# ==============================================================================
# PLOTTING (For the primary base network)
# ==============================================================================
moduleColors <- net_base$colors
MEs <- net_base$MEs
moduleTraitCor <- cor(MEs, traitMatrix, use = "complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

pdf(here("plots", "02_WGCNA_dendrogram.pdf"), width = 10, height = 6)
par(mar = c(6, 8.5, 3, 3))
plotDendroAndColors(net_base$dendrograms[[1]], moduleColors[net_base$blockGenes[[1]]], 
                    "Module Colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

pdf(here("plots", "02_WGCNA_module_trait_heatmap.pdf"), width = 12, height = 8)
age_cols <- grep("^Age_", colnames(traitMatrix))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(traitMatrix), yLabels = names(MEs), ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = ifelse(moduleTraitPvalue < 0.05, paste0(round(moduleTraitCor, 2), "*"), round(moduleTraitCor, 2)),
               setStdMargins = FALSE, cex.text = 0.7, zlim = c(-1, 1), main = "Module-Trait Relationships",
               xLabelsAngle = 45, xColorLabels = ifelse(colnames(traitMatrix) %in% colnames(traitMatrix)[age_cols], "blue", "black"))
dev.off()

saveRDS(datExpr, here("results", "02_datExpr.rds"))
saveRDS(moduleColors, here("results", "02_moduleColors.rds"))

cat("\nScript 02 completed successfully! All 4 variants are saved in the results folder.\n")