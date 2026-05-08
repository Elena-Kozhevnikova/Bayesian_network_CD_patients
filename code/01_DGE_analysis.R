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
current_dataset <- "china_transcript_tpm.csv" # тут надо менять
count_data <- read_csv(here("data", current_dataset), show_col_types = FALSE) %>%
  group_by(gene_name, sample_id) %>%
  summarise(tpm = sum(tpm), .groups = 'drop') %>%
  pivot_wider(names_from = sample_id, values_from = tpm, values_fill = 0) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

count_samples <- colnames(count_data)

# 3. Load and prepare metadata
current_metadata <- "china_metadata.csv" # тут надо менять
needed_cols <- c("Patient_group", "Age", "Gender", "CRP_mg_L", "Calprotectin_ug_g")
metadata <- read.csv(here("data", current_metadata), row.names = "pn_ID")
metadata <- metadata[, intersect(names(metadata), needed_cols)]

metadata <- metadata %>%
  mutate(
    Patient_group = case_when(
      Patient_group == "China CD" ~ "China_CD",
      grepl("China", Patient_group) ~ "China_control", # Все остальные, где есть слово China
      TRUE ~ as.character(Patient_group)             
    ),
    Patient_group = as.factor(Patient_group),
    Gender = as.factor(Gender),
    Age_merged = case_when(
      Age %in% c("15-19", "20-29") ~ "15-29",  
      Age %in% c("60-69", "70-79") ~ "60-79",
      TRUE ~ as.character(Age)
    ),
    Age = factor(Age_merged, levels = c("15-29", "30-39", "40-49", "50-59", "60-79"))
  )

# Clean NAs and synchronize samples
metadata[is.na(metadata)] <- "NA"
common_samples <- intersect(colnames(count_data), rownames(metadata))

count_data <- count_data[, common_samples]
metadata <- metadata[common_samples, ]

metadata <- droplevels(metadata) # удаляем уровни, для которых не осталось данных после объединения

# 4. Differential Gene Expression (edgeR)
cat("Running edgeR...\n")
y <- DGEList(counts = count_data, group = metadata$Patient_group)

# Calculate means for output
mean_expr <- data.frame(
  gene_id = rownames(count_data),
  mean_China_CD = rowMeans(count_data[, metadata$Patient_group == "China_CD"]),
  mean_China_control = rowMeans(count_data[, metadata$Patient_group == "China_control"])
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

contrast <- makeContrasts("Patient_groupChina_CD - Patient_groupChina_control", levels = design)
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

# Определяем путь к папке
results_dir <- here("results", "China_results")

# Если папки нет, создаем её
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  cat("Created directory:", results_dir, "\n")
}

write.csv(results, file.path(results_dir, "01_full_dge_results.csv"), row.names = TRUE)
write.csv(sig_genes, file.path(results_dir, "01_significant_genes.csv"), row.names = TRUE)

deg_counts <- count_data[rownames(sig_genes), ]
saveRDS(deg_counts, file.path(results_dir, "01_deg_counts_for_wgcna.rds"))
saveRDS(metadata, file.path(results_dir, "01_cleaned_metadata.rds"))
saveRDS(count_data, file.path(results_dir, "01_full_count_matrix.rds"))

cat("Script 01 completed successfully!\n")