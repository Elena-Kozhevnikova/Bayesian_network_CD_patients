suppressWarnings({
  suppressPackageStartupMessages({
    library(bnlearn)
    library(tidyverse)
    library(here)
    library(parallel)
  })
})

cat("--- Script 03b: BN bootstrap on all 1116 DEGs ---\n")

# 1. Загрузка данных
full_counts <- readRDS(here("results", "full_count_matrix.rds"))
sig_genes   <- read.csv(here("results", "significant_genes.csv"), row.names = 1)
metadata    <- readRDS(here("results", "cleaned_metadata.rds"))

# 2. Подготовка матрицы
deg_genes      <- rownames(sig_genes)
available      <- intersect(deg_genes, rownames(full_counts))
cat("DEGs available:", length(available), "\n")

expr_mat       <- t(full_counts[available, ])
common_samples <- intersect(rownames(expr_mat), rownames(metadata))
expr_mat       <- expr_mat[common_samples, ]

gene_vars <- apply(expr_mat, 2, var)
expr_mat  <- expr_mat[, gene_vars > 0]
cat("Genes after variance filter:", ncol(expr_mat), "\n")

bn_data <- as.data.frame(expr_mat)

# 3. Bootstrap
n_cores <- max(1, detectCores() - 1)
cat("Using", n_cores, "cores, R = 50 bootstraps\n")

cl <- makeCluster(n_cores)
clusterExport(cl, "bn_data")

time_start <- proc.time()

boot_res <- boot.strength(
  data          = bn_data,
  R             = 5,
  algorithm     = "hc",
  algorithm.args = list(score = "bic-g"),
  cluster       = cl
)

time_end <- proc.time()
stopCluster(cl)

elapsed <- time_end - time_start
cat(sprintf("\nВремя: %.1f сек (%.1f мин)\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))

# 4. Усреднённая сеть
sig_edges  <- boot_res[boot_res$strength > 0.5 & boot_res$direction > 0.5, ]
cat("Значимых рёбер (strength > 0.5):", nrow(sig_edges), "\n")

cons_bn <- averaged.network(sig_edges, threshold = 0.5)
cat("Узлов в итоговой сети:", length(nodes(cons_bn)), "\n")
cat("Рёбер в итоговой сети:", nrow(arcs(cons_bn)), "\n")

# 5. Сохраняем
saveRDS(cons_bn,  here("results", "BN_ALL_DEGs_boot50_structure.rds"))
saveRDS(boot_res, here("results", "BN_ALL_DEGs_boot50_boot.rds"))

write.csv(as.data.frame(arcs(cons_bn)),
          here("results", "BN_ALL_DEGs_boot50_arcs.csv"),
          row.names = FALSE)

cat("\nScript completed!\n")