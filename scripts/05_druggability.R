# =============================================================================
# 05_druggability.R
# Step 5: Overlay druggability scores and identify therapeutic targets
# =============================================================================

library(tidyverse)

# --- Load previous workspace ---
load("results/04_enrichment_done.RData")

dir.create("figures", showWarnings = FALSE)

# --- Merge druggability with centrality ---
drug_by_gene <- druggability %>%
  group_by(gene_name) %>%
  summarise(
    druggability_score = max(final_score, na.rm = TRUE),
    drug_method = paste(unique(method), collapse = "; "),
    .groups = "drop"
  )

target_df <- centrality_df %>%
  left_join(drug_by_gene, by = c("gene" = "gene_name")) %>%
  mutate(
    druggability_score = replace_na(druggability_score, 0),
    is_hub = degree >= quantile(degree, 0.9),
    is_bridge = gene %in% bridges$gene,
    is_core = evidence_count >= 3,
    is_druggable = druggability_score > 0
  )

# --- Composite target score ---
# Combines network importance + evidence strength + druggability
target_df <- target_df %>%
  mutate(
    degree_norm = (degree - min(degree)) / (max(degree) - min(degree)),
    between_norm = (betweenness - min(betweenness)) / (max(betweenness) - min(betweenness)),
    evidence_norm = (evidence_count - 1) / 4,  # scale 1-5 to 0-1
    drug_norm = pmin(druggability_score / max(druggability_score[druggability_score > 0], na.rm = TRUE), 1),
    target_score = 0.3 * degree_norm + 0.2 * between_norm + 0.3 * evidence_norm + 0.2 * drug_norm
  ) %>%
  arrange(desc(target_score))

# --- Top therapeutic targets ---
cat("=== Top 25 Therapeutic Targets (Composite Score) ===\n")
cat(sprintf("%-12s %5s %6s %10s %5s %8s %8s\n",
            "Gene", "Deg", "Evid", "Betw", "Comm", "DrugScr", "TargScr"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (i in 1:min(25, nrow(target_df))) {
  r <- target_df[i, ]
  cat(sprintf("%-12s %5d %6d %10.4f %5d %8.2f %8.3f\n",
              r$gene, r$degree, r$evidence_count, r$betweenness,
              r$community, r$druggability_score, r$target_score))
}

# --- Categorize targets ---
cat("\n=== Target Categories ===\n")

# Category 1: High-confidence druggable hubs (the best targets)
cat1 <- target_df %>% filter(is_hub & is_druggable & evidence_count >= 2)
cat(sprintf("\n1. Druggable Hub Genes (top 10%% degree + druggable + ≥2 evidence): %d genes\n", nrow(cat1)))
if (nrow(cat1) > 0) {
  for (i in 1:nrow(cat1)) {
    cat(sprintf("   %s (degree=%d, evidence=%d/5, drug_score=%.2f, module=%d)\n",
                cat1$gene[i], cat1$degree[i], cat1$evidence_count[i],
                cat1$druggability_score[i], cat1$community[i]))
  }
}

# Category 2: Bridge genes (connect modules — potential for combination therapy)
cat2 <- target_df %>% filter(is_bridge)
cat(sprintf("\n2. Bridge Genes (connect different resistance modules): %d genes\n", nrow(cat2)))
if (nrow(cat2) > 0) {
  for (i in 1:nrow(cat2)) {
    cat(sprintf("   %s (degree=%d, betweenness=%.4f, druggable=%s, module=%d)\n",
                cat2$gene[i], cat2$degree[i], cat2$betweenness[i],
                ifelse(cat2$is_druggable[i], "YES", "no"), cat2$community[i]))
  }
}

# Category 3: Core multi-evidence genes (≥3 layers)
cat3 <- target_df %>% filter(is_core)
cat(sprintf("\n3. Core Multi-Evidence Genes (≥3/5 layers): %d genes\n", nrow(cat3)))
for (i in 1:nrow(cat3)) {
  cat(sprintf("   %s (evidence=%d/5, degree=%d, druggable=%s, module=%d)\n",
              cat3$gene[i], cat3$evidence_count[i], cat3$degree[i],
              ifelse(cat3$is_druggable[i], "YES", "no"), cat3$community[i]))
}

# --- Figure: Target landscape scatter plot ---
p_landscape <- ggplot(target_df, aes(x = degree, y = betweenness)) +
  geom_point(aes(size = evidence_count, 
                 color = druggability_score,
                 shape = is_core),
             alpha = 0.7) +
  geom_text(
    data = target_df %>% filter(target_score >= quantile(target_score, 0.92)),
    aes(label = gene),
    size = 3, nudge_y = 0.005, check_overlap = TRUE
  ) +
  scale_color_gradient2(low = "grey70", mid = "orange", high = "red",
                        midpoint = 2, name = "Druggability\nScore") +
  scale_size_continuous(range = c(1.5, 6), name = "Evidence\nLayers") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     name = "Core Gene\n(≥3 layers)") +
  labs(
    title = "Therapeutic Target Landscape",
    subtitle = "Node position = network importance; Color = druggability; Size = evidence support",
    x = "Degree Centrality",
    y = "Betweenness Centrality"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

ggsave("figures/target_landscape.png", p_landscape, width = 12, height = 8, dpi = 300)
cat("\nSaved: figures/target_landscape.png\n")

# --- Figure: Top 20 targets bar chart ---
top20 <- target_df %>% head(20)

p_bars <- ggplot(top20, aes(x = reorder(gene, target_score), y = target_score)) +
  geom_col(aes(fill = factor(community)), width = 0.7) +
  geom_text(aes(label = sprintf("%d/5", evidence_count)),
            hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2", name = "Module") +
  labs(
    title = "Top 20 Therapeutic Target Candidates",
    subtitle = "Composite score: network centrality + evidence support + druggability",
    x = NULL,
    y = "Target Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave("figures/top_targets_barplot.png", p_bars, width = 10, height = 7, dpi = 300)
cat("Saved: figures/top_targets_barplot.png\n")

# --- Figure: Module-resolved target heatmap ---
module_summary <- target_df %>%
  group_by(community) %>%
  summarise(
    n_genes = n(),
    n_druggable = sum(is_druggable),
    pct_druggable = 100 * mean(is_druggable),
    mean_degree = mean(degree),
    top_hub = gene[which.max(degree)],
    n_core = sum(is_core),
    .groups = "drop"
  )

cat("\n=== Module Druggability Summary ===\n")
cat(sprintf("%-10s %6s %10s %12s %10s %12s %6s\n",
            "Module", "Genes", "Druggable", "% Druggable", "Avg Degree", "Top Hub", "Core"))
for (i in 1:nrow(module_summary)) {
  r <- module_summary[i, ]
  cat(sprintf("%-10d %6d %10d %11.1f%% %10.1f %12s %6d\n",
              r$community, r$n_genes, r$n_druggable, r$pct_druggable,
              r$mean_degree, r$top_hub, r$n_core))
}

# --- Save results ---
write_csv(target_df, "results/target_analysis.csv")
write_csv(module_summary, "results/module_summary.csv")

save(g, target_df, module_summary, centrality_df, communities, bridges,
     go_results_list, kegg_results_list,
     real_clustering, mean_random_clustering, real_path, mean_random_path, sigma,
     evidence, consensus_2plus, consensus_3plus, core_genes,
     druggability, edges, nodes, network_genes,
     file = "results/05_druggability_done.RData")

cat("\n=== Files saved ===\n")
cat("  results/target_analysis.csv\n")
cat("  results/module_summary.csv\n")
cat("  results/05_druggability_done.RData\n")
cat("\nReady for 06_visualization.R\n")