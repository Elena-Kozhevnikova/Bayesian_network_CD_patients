# ==============================================================================
# Script 04: Bayesian Network Visualization (Rgraphviz Implementation)
# Project: Bayes_26
# ==============================================================================

suppressWarnings({
  suppressPackageStartupMessages({
    library(bnlearn)
    library(Rgraphviz)
    library(tidyverse)
    library(here)
    library(wesanderson)
  })
})

# ==============================================================================
# SETTINGS
# ==============================================================================
current_prefix <- "hubs_5_modules_20_genes"
network_type <- "CALPROTECTIN" # or "GENES_ONLY"

cat("Rendering network:", network_type, "for dataset:", current_prefix, "\n")

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================
cons_bn <- readRDS(here("results", paste0("BN_", network_type, "_structure_", current_prefix, ".rds")))
boot_res <- readRDS(here("results", paste0("BN_", network_type, "_boot_", current_prefix, ".rds")))
hub_genes_info <- read_csv(here("results", paste0(current_prefix, ".csv")), show_col_types = FALSE)

# ==============================================================================
# 2. COLOR CONFIGURATION
# ==============================================================================
all_nodes <- nodes(cons_bn)
AsteroidCity_colors <- wes_palette("AsteroidCity1", 5, type = "continuous")

# Map colors to WGCNA modules
module_palette <- c(
  "yellow"    = AsteroidCity_colors[2], 
  "blue"      = AsteroidCity_colors[1],
  "turquoise" = AsteroidCity_colors[4], 
  "brown"     = AsteroidCity_colors[3],
  "green"     = AsteroidCity_colors[5]
)

node_modules <- hub_genes_info$module[match(all_nodes, hub_genes_info$gene)]
node_fill_colors <- module_palette[node_modules]

# Highlight Calprotectin with a distinct color
node_fill_colors[all_nodes == "Calprotectin_ug_g"] <- "#FF4500" 
node_fill_colors[is.na(node_fill_colors)] <- "grey90"
names(node_fill_colors) <- all_nodes

# ==============================================================================
# 3. GRAPH RENDERING BLOCK
# ==============================================================================
cat("Preparing graph layout...\n")

# Initialize graph with strength/bootstrap values
g <- strength.plot(cons_bn, boot_res, threshold = 0.7, render = FALSE)

# Global layout settings with compression focus
g <- layoutGraph(g, layoutType = "dot", 
                 attrs = list(
                   graph = list(
                     rankdir = "TB", 
                     nodesep = "0.3",    # Minimal distance between nodes in a row
                     ranksep = "0.5",    # Minimal distance between levels
                     ratio = "compress", # Force layout compression
                     size = "15,15" 
                   ),
                   node = list(
                     shape = "ellipse", 
                     fixedsize = FALSE, 
                     fontsize = 20,
                     margin = "1.0,0.5"  # Optimal node padding
                   )
                 )
)

# Apply colors and styles (preserving geometric layout)
nodeRenderInfo(g)$fill <- node_fill_colors
nodeRenderInfo(g)$col <- "black"
nodeRenderInfo(g)$lwd <- 2

# ==============================================================================
# 4. EXPORT
# ==============================================================================
out_file <- here("plots", paste0("BN_final_", network_type, "_", current_prefix, ".pdf"))
cat("Saving plot to:", out_file, "\n")

pdf(out_file, width = 20, height = 20)
renderGraph(g)
dev.off()

cat("Script 04 completed successfully! Check the plots folder.\n")