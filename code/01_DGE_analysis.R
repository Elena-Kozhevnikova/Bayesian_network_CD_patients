# ==============================================================================
# Script 01: Data Preparation and Differential Gene Expression (DGE)
# Project: Bayes_26
# ==============================================================================

# 1. Load libraries
suppressWarnings({
  suppressPackageStartupMessages({
    library(edgeR)
    library(tidyverse)
    library(here) 
  })
})

# 2. Load and prepare counts
cat("Loading data...\n")
count_data <- read_csv(here("data", "israel_transcript_tpm.csv"), show_col_types = FALSE) %>%
  group_by(gene_name, sample_id) %>%
  summarise(tpm = sum(tpm), .groups = 'drop') %>%
  pivot_wider(names_from = sample_id, values_from = tpm, values_fill = 0) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

count_samples <- colnames(count_data)

# 3. Load and prepare metadata
metadata <- read.csv(here("data", "Israel_metadata.csv"), row.names = "pn_ID")

metadata <- metadata %>%
  mutate(
    Patient_group = as.factor(Patient_group),
    Gender = as.factor(Gender),
    Age_merged = case_when(
      Age %in% c("15_19", "20_29") ~ "15_29",
      Age %in% c("60_69", "70_79") ~ "60_79",
      TRUE ~ as.character(Age)
    ),
    Age = factor(Age_merged, levels = c("15_29", "30_39", "40_49", "50_59", "60_79"))
  )

# Clean NAs and synchronize samples
metadata[is.na(metadata)] <- "NA"
common_samples <- intersect(colnames(count_data), rownames(metadata))

count_data <- count_data[, common_samples]
metadata <- metadata[common_samples, ]

# 4. Differential Gene Expression (edgeR)
cat("Running edgeR...\n")
y <- DGEList(counts = count_data, group = metadata$Patient_group)

# Calculate means for output
mean_expr <- data.frame(
  gene_id = rownames(count_data),
  mean_Israel_CD = rowMeans(count_data[, metadata$Patient_group == "Israel_CD"]),
  mean_Israel_control = rowMeans(count_data[, metadata$Patient_group == "Israel_control"])
)

# Filtering and normalization
keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# Design matrix and calculations
design <- model.matrix(~0 + Patient_group + Age + Gender, data = metadata)
colnames(design) <- make.names(colnames(design))

y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design)

contrast <- makeContrasts("Patient_groupIsrael_CD - Patient_groupIsrael_control", levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

# 5. Assemble results
results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
results$gene_id <- rownames(results)
results <- merge(results, mean_expr, by = "gene_id")
rownames(results) <- results$gene_id
results <- results[order(results$PValue), ]

sig_genes <- results %>%
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  arrange(FDR)

# 6. Save results for downstream scripts
cat("Saving results...\n")

write.csv(results, here("results", "full_dge_results.csv"), row.names = TRUE)
write.csv(sig_genes, here("results", "significant_genes.csv"), row.names = TRUE)

deg_counts <- count_data[rownames(sig_genes), ]
saveRDS(deg_counts, here("results", "deg_counts_for_wgcna.rds"))
saveRDS(metadata, here("results", "cleaned_metadata.rds"))
saveRDS(count_data, here("results", "full_count_matrix.rds"))

cat("Script 01 completed successfully!\n")