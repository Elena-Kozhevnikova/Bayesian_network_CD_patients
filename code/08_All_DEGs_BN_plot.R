suppressWarnings({
  suppressPackageStartupMessages({
    library(bnlearn)
    library(igraph)
    library(tidyverse)
    library(here)
  })
})

# 1. Загружаем структуру
bn_struct <- readRDS(here("results", "BN_ALL_DEGs_structure.rds"))
arcs_df   <- as.data.frame(arcs(bn_struct))

cat("Всего узлов:", length(nodes(bn_struct)), "\n")
cat("Всего рёбер:", nrow(arcs_df), "\n")

# 2. Считаем degree каждого гена
g <- graph_from_data_frame(arcs_df, directed = TRUE)
deg <- degree(g, mode = "all")
in_deg  <- degree(g, mode = "in")
out_deg <- degree(g, mode = "out")

deg_df <- data.frame(
  gene    = names(deg),
  degree  = deg,
  in_deg  = in_deg[names(deg)],
  out_deg = out_deg[names(deg)]
) %>% arrange(desc(degree))

# Топ-20 хабов
cat("\nТоп-20 генов по degree:\n")
print(head(deg_df, 20))

write.csv(deg_df, here("results", "BN_ALL_DEGs_degree.csv"), row.names = FALSE)

# 3. Визуализация подсети топ-N генов
TOP_N <- 80  # можно менять: 50, 100, 150

top_genes  <- head(deg_df$gene, TOP_N)
sub_arcs   <- arcs_df %>% filter(from %in% top_genes & to %in% top_genes)
g_sub      <- graph_from_data_frame(sub_arcs, directed = TRUE, 
                                    vertices = top_genes)

# Цвет узлов по out-degree (много исходящих = регулятор)
pal        <- colorRampPalette(c("#cce5f0", "#01696f"))(100)
deg_norm   <- (out_deg[top_genes] - min(out_deg[top_genes])) /
  (max(out_deg[top_genes]) - min(out_deg[top_genes]) + 1)
node_cols  <- pal[ceiling(deg_norm * 99) + 1]

# Размер узла по общему degree
node_sizes <- 3 + 10 * (deg[top_genes] / max(deg[top_genes]))

pdf(here("plots", "BN_ALL_DEGs_top80.pdf"), width = 20, height = 20)
plot(g_sub,
     layout          = layout_with_fr(g_sub),
     vertex.size     = node_sizes,
     vertex.color    = node_cols,
     vertex.label    = V(g_sub)$name,
     vertex.label.cex= 0.6,
     vertex.label.color = "black",
     edge.arrow.size = 0.3,
     edge.color      = "grey60",
     main            = paste0("BN: Top ", TOP_N, " hub genes (by degree)"))
dev.off()

cat("График сохранён в plots/BN_ALL_DEGs_top80.pdf\n")