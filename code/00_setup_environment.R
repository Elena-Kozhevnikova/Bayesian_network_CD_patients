# ==============================================================================
# Script 00: Environment Setup (Run ONLY once on a new machine)
# Project: Bayes_26
# ==============================================================================

cat("Starting installation of base CRAN packages...\n")

cran_packages <- c(
  "tidyverse", "here", "bnlearn", "igraph", 
  "scales", "tidyr", "tibble", "ggpubr", "lmtest", "wesanderson"
)

# Function to install missing CRAN packages
install_missing_cran <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
}
install_missing_cran(cran_packages)


cat("Starting installation of Bioconductor packages...\n")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c(
  "edgeR", 
  "WGCNA", 
  "org.Hs.eg.db", 
  "clusterProfiler",
  "Rgraphviz" # Rgraphviz officially resides in Bioconductor
)

# Function to install missing Bioconductor packages
install_missing_bioc <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new_pkgs)) BiocManager::install(new_pkgs, update = FALSE, ask = FALSE)
}
install_missing_bioc(bioc_packages)

cat("All packages installed successfully! You can now proceed to Script 01.\n")