# =============================================================================
# 06_visualization.R
# Step 6: Network visualization and summary figures
# =============================================================================

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

# --- Load previous workspace ---
load("results/05_druggability_done.RData")

dir.create("figures", showWarnings = FALSE)

# --- Prepare tidygraph object ---
tg <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  left_join(target_df, by = c("name" = "gene"))

# --- FIGURE 1: Main network plot ---
# Color by module, size by degree, label hubs and core genes
set.seed(42)

# Genes to label: top 5 per module by degree + all core genes
label_genes <- target_df %>%
  group_by(community) %>%
  slice_max(degree, n = 3) %>%
  pull(gene)
label_genes <- unique(c(label_genes, core_genes$gene, bridges$gene[1:5]))

p_network <- ggraph(tg, layout = "fr") +
  geom_edge_link(alpha = 0.06, color = "grey60") +
  geom_node_point(aes(size = degree,
                      fill = factor(community),
                      shape = ifelse(evidence_count >= 3, "core", "other")),
                  color = "black", stroke = 0.3) +
  geom_node_text(
    aes(label = ifelse(name %in% label_genes, name, "")),
    size = 2.5, repel = TRUE, max.overlaps = 25,
    fontface = "bold"
  ) +
  scale_fill_brewer(palette = "Set2", name = "Module",
                    labels = c("1: Transcription",
                               "2: Cell Cycle",
                               "3: Metabolism",
                               "4: Immune",
                               "5: Cytoskeleton",
                               "6: Translation")) +
  scale_size_continuous(range = c(1, 10), name = "Degree") +
  scale_shape_manual(values = c("core" = 24, "other" = 21),
                     name = "Evidence",
                     labels = c("≥3 layers", "2 layers")) +
  guides(fill = guide_legend(override.aes = list(size = 5, shape = 21))) +
  labs(title = "CML TKI Resistance Network",
       subtitle = "184 nodes | 1106 edges | 6 functional modules | 5 evidence layers") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.position = "right"
  )

ggsave("figures/network_main.png", p_network, width = 16, height = 12, dpi = 300)
cat("Saved: figures/network_main.png\n")

# --- FIGURE 2: Hub genes bar chart with module colors ---
top30 <- target_df %>% head(30)

p_hubs <- ggplot(top30, aes(x = reorder(gene, degree), y = degree)) +
  geom_col(aes(fill = factor(community)), width = 0.7) +
  geom_point(aes(y = degree + 2, size = evidence_count), shape = 16, color = "black") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2", name = "Module") +
  scale_size_continuous(range = c(1, 4), name = "Evidence\nLayers") +
  labs(
    title = "Top 30 Hub Genes by Degree Centrality",
    x = NULL,
    y = "Degree"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 13))

ggsave("figures/hub_genes_barplot.png", p_hubs, width = 10, height = 8, dpi = 300)
cat("Saved: figures/hub_genes_barplot.png\n")

# --- FIGURE 3: Small-world comparison histograms ---
sw_df <- data.frame(
  clustering = random_clustering,
  path_length = random_path_length
)

p_sw1 <- ggplot(sw_df, aes(x = clustering)) +
  geom_histogram(bins = 40, fill = "grey70", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = real_clustering, color = "red", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = real_clustering, y = Inf,
           label = sprintf("Real = %.3f", real_clustering),
           vjust = 2, hjust = -0.1, color = "red", fontface = "bold", size = 3.5) +
  labs(title = "Clustering Coefficient", subtitle = "vs 1000 random networks",
       x = "Clustering Coefficient", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

p_sw2 <- ggplot(sw_df %>% filter(!is.na(path_length)), aes(x = path_length)) +
  geom_histogram(bins = 40, fill = "grey70", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = real_path, color = "red", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = real_path, y = Inf,
           label = sprintf("Real = %.3f", real_path),
           vjust = 2, hjust = -0.1, color = "red", fontface = "bold", size = 3.5) +
  labs(title = "Average Path Length", subtitle = "vs 1000 random networks",
       x = "Average Path Length", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

library(patchwork)
p_sw_combined <- p_sw1 + p_sw2 +
  plot_annotation(
    title = sprintf("Small-World Network Test (σ = %.2f)", sigma),
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave("figures/small_world_test.png", p_sw_combined, width = 12, height = 5, dpi = 300)
cat("Saved: figures/small_world_test.png\n")

# --- FIGURE 4: Evidence × Module heatmap ---
ev_mod <- target_df %>%
  group_by(community) %>%
  summarise(
    L1_Krishnan_scRNAseq = sum(L1_Krishnan_scRNAseq),
    L2_Awad_CRISPR_KO = sum(L2_Awad_CRISPR_KO),
    L3_Krishnan_ML = sum(L3_Krishnan_ML),
    L4_Sacco_RNAseq = sum(L4_Sacco_RNAseq),
    L5_Sacco_Phospho = sum(L5_Sacco_Phospho),
    .groups = "drop"
  ) %>%
  column_to_rownames("community")

colnames(ev_mod) <- c("Krishnan\nscRNA-seq", "Awad\nCRISPR KO",
                      "Krishnan\nML", "Sacco\nRNA-seq", "Sacco\nPhospho")
rownames(ev_mod) <- c("M1: Transcription", "M2: Cell Cycle",
                      "M3: Metabolism", "M4: Immune",
                      "M5: Cytoskeleton", "M6: Translation")

png("figures/evidence_module_heatmap.png", width = 8, height = 5, units = "in", res = 300)
pheatmap(
  as.matrix(ev_mod),
  color = colorRampPalette(c("white", "#2F5496"))(50),
  display_numbers = TRUE,
  number_format = "%d",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 10,
  fontsize_number = 10,
  main = "Evidence Layer Contribution per Module\n(number of genes from each layer)"
)
dev.off()
cat("Saved: figures/evidence_module_heatmap.png\n")

# --- Save final workspace ---
save.image(file = "results/06_final.RData")

cat("\n============================================================\n")
cat("ALL ANALYSIS AND VISUALIZATION COMPLETE\n")
cat("============================================================\n")
cat("\nFigures generated:\n")
cat("  figures/network_main.png           — Main network visualization\n")
cat("  figures/hub_genes_barplot.png      — Top 30 hub genes\n")
cat("  figures/small_world_test.png       — Random network comparison\n")
cat("  figures/evidence_module_heatmap.png — Evidence × Module breakdown\n")
cat("  figures/target_landscape.png       — Therapeutic target scatter\n")
cat("  figures/top_targets_barplot.png    — Top 20 target candidates\n")
cat("  figures/go_enrichment_*.png        — GO enrichment per module\n")
cat("  figures/upset_plot.png             — Evidence intersection\n")
cat("  figures/evidence_distribution.png  — Evidence count histogram\n")
cat("  figures/top_genes_heatmap.png      — Core genes heatmap\n")
cat("  figures/pairwise_overlap_heatmap.png — Layer overlap matrix\n")
cat("\nProject complete. Ready for README and GitHub.\n")