# =============================================================================
# 04_enrichment.R
# Step 4: GO and KEGG pathway enrichment per module
# =============================================================================

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# --- Load previous workspace ---
load("results/03_network_analyzed.RData")

# --- Convert gene symbols to Entrez IDs ---
# clusterProfiler needs Entrez IDs for enrichment
gene_mapping <- bitr(
  centrality_df$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cat(sprintf("Mapped %d / %d gene symbols to Entrez IDs\n",
            nrow(gene_mapping), nrow(centrality_df)))

# Add Entrez IDs to centrality data
centrality_df <- centrality_df %>%
  left_join(gene_mapping, by = c("gene" = "SYMBOL"))

# --- Prepare gene lists per module ---
# Only analyze modules with ≥10 genes (skip Module 5 with 6 genes)
module_ids <- centrality_df %>%
  count(community) %>%
  filter(n >= 10) %>%
  pull(community)

cat(sprintf("\nAnalyzing %d modules (≥10 genes each)\n", length(module_ids)))

# --- GO Biological Process enrichment per module ---
cat("\n=== GO Biological Process Enrichment ===\n")

go_results_list <- list()

for (mod_id in module_ids) {
  mod_genes <- centrality_df %>%
    filter(community == mod_id & !is.na(ENTREZID)) %>%
    pull(ENTREZID)
  
  # Background: all genes in the network
  bg_genes <- centrality_df %>%
    filter(!is.na(ENTREZID)) %>%
    pull(ENTREZID)
  
  go_res <- enrichGO(
    gene = mod_genes,
    universe = bg_genes,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  
  if (!is.null(go_res) && nrow(go_res@result %>% filter(p.adjust < 0.05)) > 0) {
    go_results_list[[paste0("Module_", mod_id)]] <- go_res
    
    top_terms <- go_res@result %>%
      filter(p.adjust < 0.05) %>%
      head(5)
    
    cat(sprintf("\nModule %d — Top GO:BP terms:\n", mod_id))
    for (j in 1:nrow(top_terms)) {
      cat(sprintf("  %s (p.adj=%.2e, %d genes)\n",
                  top_terms$Description[j],
                  top_terms$p.adjust[j],
                  top_terms$Count[j]))
    }
  } else {
    cat(sprintf("\nModule %d — No significant GO:BP terms\n", mod_id))
  }
}

# --- KEGG pathway enrichment per module ---
cat("\n=== KEGG Pathway Enrichment ===\n")

kegg_results_list <- list()

for (mod_id in module_ids) {
  mod_genes <- centrality_df %>%
    filter(community == mod_id & !is.na(ENTREZID)) %>%
    pull(ENTREZID)
  
  bg_genes <- centrality_df %>%
    filter(!is.na(ENTREZID)) %>%
    pull(ENTREZID)
  
  kegg_res <- enrichKEGG(
    gene = mod_genes,
    universe = bg_genes,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  
  if (!is.null(kegg_res) && nrow(kegg_res@result %>% filter(p.adjust < 0.05)) > 0) {
    kegg_results_list[[paste0("Module_", mod_id)]] <- kegg_res
    
    top_kegg <- kegg_res@result %>%
      filter(p.adjust < 0.05) %>%
      head(5)
    
    cat(sprintf("\nModule %d — Top KEGG pathways:\n", mod_id))
    for (j in 1:nrow(top_kegg)) {
      cat(sprintf("  %s (p.adj=%.2e, %d genes)\n",
                  top_kegg$Description[j],
                  top_kegg$p.adjust[j],
                  top_kegg$Count[j]))
    }
  } else {
    cat(sprintf("\nModule %d — No significant KEGG pathways\n", mod_id))
  }
}

# --- Generate enrichment plots ---

# Dotplot for each module with GO results
for (mod_name in names(go_results_list)) {
  mod_id <- gsub("Module_", "", mod_name)
  
  p <- dotplot(go_results_list[[mod_name]], showCategory = 15) +
    ggtitle(paste0("Module ", mod_id, " — GO Biological Process")) +
    theme(plot.title = element_text(face = "bold", size = 12))
  
  ggsave(
    filename = sprintf("figures/go_enrichment_module_%s.png", mod_id),
    plot = p, width = 10, height = 8, dpi = 300
  )
  cat(sprintf("Saved: figures/go_enrichment_module_%s.png\n", mod_id))
}

# --- Comparative dotplot: all modules side by side ---
if (length(go_results_list) >= 2) {
  # Merge results for comparison
  combined_go <- merge_result(go_results_list)
  
  p_compare <- dotplot(combined_go, showCategory = 5) +
    ggtitle("GO Biological Process — All Modules") +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.text.y = element_text(size = 8)
    )
  
  ggsave("figures/go_enrichment_all_modules.png",
         plot = p_compare, width = 14, height = 10, dpi = 300)
  cat("Saved: figures/go_enrichment_all_modules.png\n")
}

# --- Save enrichment results as CSVs ---
for (mod_name in names(go_results_list)) {
  mod_id <- gsub("Module_", "", mod_name)
  write_csv(
    go_results_list[[mod_name]]@result %>% filter(p.adjust < 0.05),
    sprintf("results/go_enrichment_module_%s.csv", mod_id)
  )
}

for (mod_name in names(kegg_results_list)) {
  mod_id <- gsub("Module_", "", mod_name)
  write_csv(
    kegg_results_list[[mod_name]]@result %>% filter(p.adjust < 0.05),
    sprintf("results/kegg_enrichment_module_%s.csv", mod_id)
  )
}

# --- Save workspace ---
save(g, centrality_df, communities, bridges,
     real_clustering, mean_random_clustering,
     real_path, mean_random_path, sigma,
     go_results_list, kegg_results_list,
     evidence, consensus_2plus, consensus_3plus, core_genes,
     druggability, edges, nodes, network_genes,
     file = "results/04_enrichment_done.RData")

cat("\n=== Enrichment analysis complete ===\n")
cat("Ready for 05_druggability.R\n")