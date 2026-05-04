# ==============================================================================
# Script 06: Predictor Influence, China Cohort Validation, and Confusion Matrix
# Project: Bayes_26
# ==============================================================================

suppressWarnings({
  suppressPackageStartupMessages({
    library(bnlearn)
    library(tidyverse)
    library(edgeR)
    library(here)
    library(wesanderson)
    library(ggsignif) # For significance brackets and stars
    library(pROC)     # For classification threshold calculation
  })
})

current_prefix <- "hubs_CLASSIC_original"
cat("--- Running Script 06: Validation and Prediction ---\n")

# 1. Load trained model and training data (Israel)
fitted_bn_calprotectin <- readRDS(here("results", paste0("BN_CALPROTECTIN_fitted_", current_prefix, ".rds")))
train_data <- read_csv(here("results", paste0("BN_training_data_", current_prefix, ".csv")), show_col_types = FALSE) %>%
  column_to_rownames("sample_id")

AsteroidCity1_colors <- wes_palette("AsteroidCity1")

# ==============================================================================
# PART 1: PREDICTOR EFFECT ASSESSMENT AND TRAINING ACCURACY
# ==============================================================================
cat("\n[1/4] Calculating predictor gene effects...\n")
parents <- fitted_bn_calprotectin$Calprotectin_ug_g$parents

if (length(parents) > 0) {
  # 1.1 Effect Size Plot
  coefs <- fitted_bn_calprotectin$Calprotectin_ug_g$coefficients
  # Calculate standardized effects (SD change in Calprotectin per SD change in gene)
  sd_effects <- coefs[-1] * sapply(train_data[, parents, drop=FALSE], sd)
  
  coef_plot_df <- data.frame(
    Gene = factor(parents, levels = parents[order(abs(sd_effects))]),
    Effect = sd_effects
  )
  
  effect_size_plot <- ggplot(coef_plot_df, aes(x = Gene, y = Effect, fill = Effect > 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(AsteroidCity1_colors[1], AsteroidCity1_colors[3])) +
    coord_flip() +
    labs(title = "Gene Effects on Calprotectin Prediction",
         y = "Effect Size (SD units change in calprotectin per SD increase in gene)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave(here("plots", paste0("06_calprotectin_parent_effect_size_", current_prefix, ".pdf")), 
         plot = effect_size_plot, width = 8, height = 6, dpi = 300, bg = "white")
  
  # 1.2 Prediction Accuracy assessment on Training Data
  cat("Evaluating linear prediction model (Train Data)...\n")
  train_preds <- predict(fitted_bn_calprotectin, node = "Calprotectin_ug_g", 
                         data = train_data, method = "bayes-lw")
  
  plot_data <- data.frame(
    Observed = train_data$Calprotectin_ug_g,
    Predicted = as.numeric(train_preds)
  )
  
  # Statistics for the plot
  lm_fit <- lm(Predicted ~ Observed, data = plot_data)
  r2_val <- summary(lm_fit)$r.squared
  p_val <- anova(lm_fit)[1, "Pr(>F)"]
  stat_label <- sprintf("R-squared = %.3f\np-value = %.2e", r2_val, p_val)
  
  # Prediction Accuracy Plot
  accuracy_plot <- ggplot(plot_data, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.6, color = AsteroidCity1_colors[3], size = 2) +
    geom_smooth(method = "lm", color = AsteroidCity1_colors[1], se = TRUE) +
    annotate("text", x = min(plot_data$Observed, na.rm = TRUE), 
             y = max(plot_data$Predicted, na.rm = TRUE), 
             label = stat_label, hjust = 0, vjust = 1, size = 5, fontface = "bold") +
    labs(title = "Prediction Accuracy (Training Set)",
         x = "Observed Calprotectin (ug/g)",
         y = "Predicted Calprotectin (ug/g)") +
    theme_minimal()
  
  ggsave(here("plots", paste0("06_prediction_accuracy_", current_prefix, ".pdf")), 
         plot = accuracy_plot, width = 6, height = 5, dpi = 300)
}

# ==============================================================================
# PART 2: DGE CALCULATION FOR CHINA COHORT
# ==============================================================================
cat("\n[2/4] Preparing China cohort data and calculating DGE...\n")
# Using new filenames as defined
metadata_china <- read_csv(here("data", "china_metadata.csv"), show_col_types = FALSE) %>% as.data.frame()
expression_data <- read_csv(here("data", "china_transcript_tpm.csv"), show_col_types = FALSE) %>% as.data.frame()

# Aggregating transcript TPM to Gene level
gene_tpm <- expression_data %>%
  group_by(gene_name, id) %>%
  summarize(tpm = sum(tpm), .groups = 'drop') %>%
  pivot_wider(names_from = id, values_from = tpm, values_fill = 0)

count_matrix <- as.matrix(gene_tpm[, -1])
rownames(count_matrix) <- gene_tpm$gene_name

metadata_china_dge <- metadata_china %>% 
  filter(pn_ID %in% colnames(count_matrix)) %>% 
  column_to_rownames("pn_ID")

# Basic edgeR pipeline for China cohort
y <- DGEList(counts = count_matrix, group = metadata_china_dge$Patient_group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0 + Patient_group, data = metadata_china_dge)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)
contrast <- makeContrasts("Patient_groupChina_urban - Patient_groupChina_CD", levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

dge_results <- topTags(qlf, n = Inf)$table %>% tibble::rownames_to_column("gene_name")
write.table(dge_results, here("results", paste0("china_dge_results_", current_prefix, ".tsv")), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# PART 3: PREDICTION ON CHINA DATA AND SIGNIFICANCE BOXPLOT
# ==============================================================================
cat("\n[3/4] Predicting Calprotectin and hypothesis testing...\n")

common_ids <- intersect(metadata_china$pn_ID, expression_data$id)
merged_data <- metadata_china %>% 
  filter(pn_ID %in% common_ids) %>% 
  left_join(
    expression_data %>% filter(id %in% common_ids) %>% 
      group_by(id, gene_name) %>% summarise(tpm = sum(tpm), .groups = 'drop') %>% 
      pivot_wider(names_from = gene_name, values_from = tpm, values_fill = 0),
    by = c("pn_ID" = "id")
  )

# Select predictor genes present in China dataset
genes_to_select <- setdiff(nodes(fitted_bn_calprotectin), "Calprotectin_ug_g")
available_genes <- genes_to_select[genes_to_select %in% colnames(merged_data)]
prediction_data <- merged_data[, available_genes, drop = FALSE]

# Predict Calprotectin levels for China cohort using the Israel-trained BN
predictions <- predict(fitted_bn_calprotectin, node = "Calprotectin_ug_g", 
                       data = prediction_data, method = "bayes-lw")

results_china <- merged_data %>%
  dplyr::select(pn_ID, Patient_group, Age, Gender) %>%
  dplyr::mutate(Predicted_Calprotectin = as.numeric(predictions)) %>%
  mutate(Plot_group = factor(case_when(
    Patient_group %in% c("China_urban", "China_rural") ~ "China_Control",
    TRUE ~ as.character(Patient_group)
  ), levels = c("China_Control", "China_CD")))

# Boxplot with significance stars (ggsignif)
my_fill_colors <- c("China_Control" = "#E69F00", "China_CD" = "#56B4E9")

boxplot_signif <- ggplot(results_china, aes(x = Plot_group, y = Predicted_Calprotectin)) +
  geom_boxplot(aes(fill = Plot_group), outlier.shape = NA) +
  geom_jitter(aes(fill = Plot_group), shape = 21, color = "black", width = 0.2, alpha = 0.7, size = 3) +
  scale_fill_manual(values = my_fill_colors) +
  # Add significance bracket and stars via Wilcoxon test
  geom_signif(comparisons = list(c("China_Control", "China_CD")), 
              test = "wilcox.test", 
              map_signif_level = TRUE, 
              textsize = 6, vjust = 0.5, margin_top = 0.1) +
  labs(title = "Predicted Calprotectin by Patient Group (China Cohort)",
       x = "Patient Group",
       y = "Predicted Calprotectin (ug/g)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(here("plots", paste0("06_predicted_calprotectin_signif_", current_prefix, ".pdf")), 
       plot = boxplot_signif, width = 5, height = 6, dpi = 300)

# ==============================================================================
# PART 4: CONFUSION MATRIX FOR CHINA COHORT
# ==============================================================================
cat("\n[4/4] Calculating ROC curve and Confusion Matrix...\n")

# Build ROC to find optimal Calprotectin threshold separating CD and Control
roc_obj <- roc(results_china$Plot_group, results_china$Predicted_Calprotectin, 
               levels = c("China_Control", "China_CD"), direction = "<", quiet = TRUE)

optimal_coords <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE)
best_thresh <- optimal_coords$threshold[1]

cat(sprintf("Optimal Calprotectin threshold for classification: %.2f\n", best_thresh))

# Classify patients based on optimal threshold
results_china <- results_china %>%
  mutate(Predicted_Diagnosis = factor(ifelse(Predicted_Calprotectin >= best_thresh, "China_CD", "China_Control"),
                                      levels = c("China_Control", "China_CD")))

# Generate Confusion Matrix
cm_table <- table(Actual = results_china$Plot_group, Predicted = results_china$Predicted_Diagnosis)
cm_df <- as.data.frame(cm_table)

# Plot Confusion Matrix Heatmap
conf_matrix_plot <- ggplot(cm_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = Freq), color = "white", size = 10, fontface = "bold") +
  scale_fill_gradient(low = "#56B4E9", high = "#E69F00") +
  labs(title = "Confusion Matrix: Prediction vs Actual Diagnosis",
       subtitle = sprintf("Based on optimal Calprotectin cutoff (%.1f)", best_thresh),
       x = "Predicted Diagnosis (BN Model)",
       y = "Actual Diagnosis") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

ggsave(here("plots", paste0("06_confusion_matrix_", current_prefix, ".pdf")), 
       plot = conf_matrix_plot, width = 6, height = 5, dpi = 300)

cat("\nScript 06 completed successfully! Validation results and plots saved.\n")