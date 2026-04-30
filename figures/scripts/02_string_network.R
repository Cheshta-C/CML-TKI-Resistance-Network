# =============================================================================
# 02_string_network.R
# Step 2: Query STRING-db and build the PPI network
# =============================================================================

library(tidyverse)
library(STRINGdb)

# --- Load previous workspace ---
load("results/01_data_loaded.RData")

# --- Initialize STRING-db ---
# species 9606 = Homo sapiens
# score_threshold 400 = medium confidence (range: 0-1000)
# version "12.0" = latest
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400,
  network_type = "full"    # "full" includes all interaction types
)

cat("STRING-db initialized\n")

# --- Prepare gene list for mapping ---
# STRING needs a data frame with a column called "gene" (or you specify the column)
consensus_df <- data.frame(gene = consensus_2plus, stringsAsFactors = FALSE)

cat(sprintf("Mapping %d consensus genes to STRING identifiers...\n", nrow(consensus_df)))

# --- Map gene symbols to STRING IDs ---
mapped <- string_db$map(consensus_df, "gene", removeUnmappedRows = TRUE)

cat(sprintf("Successfully mapped: %d / %d genes (%.1f%%)\n",
            nrow(mapped), length(consensus_2plus),
            100 * nrow(mapped) / length(consensus_2plus)))

# Check unmapped genes
unmapped <- setdiff(consensus_2plus, mapped$gene)
if (length(unmapped) > 0) {
  cat(sprintf("\nUnmapped genes (%d): %s\n", length(unmapped), paste(unmapped, collapse=", ")))
}

# --- Get interactions between mapped genes ---
interactions <- string_db$get_interactions(mapped$STRING_id)

cat(sprintf("\n=== Network Statistics ===\n"))
cat(sprintf("Interactions (edges): %d\n", nrow(interactions)))
cat(sprintf("Unique proteins in network: %d\n",
            length(unique(c(interactions$from, interactions$to)))))

# --- Convert STRING IDs back to gene symbols ---
# Create a lookup table: STRING_id -> gene symbol
id_to_gene <- setNames(mapped$gene, mapped$STRING_id)

edges <- interactions %>%
  mutate(
    gene1 = id_to_gene[from],
    gene2 = id_to_gene[to]
  ) %>%
  filter(!is.na(gene1) & !is.na(gene2)) %>%
  select(gene1, gene2, combined_score)

cat(sprintf("Edges after gene symbol mapping: %d\n", nrow(edges)))

# --- Build node attribute table ---
# Genes that appear in at least one interaction
network_genes <- unique(c(edges$gene1, edges$gene2))

nodes <- evidence %>%
  filter(gene %in% network_genes) %>%
  select(gene, everything())

cat(sprintf("Nodes in final network: %d\n", nrow(nodes)))

# --- Add druggability info to nodes ---
nodes <- nodes %>%
  left_join(
    druggability %>%
      group_by(gene_name) %>%
      summarise(druggability_score = max(final_score, na.rm = TRUE), .groups = "drop"),
    by = c("gene" = "gene_name")
  )

# --- Identify isolated consensus genes (mapped but no interactions) ---
mapped_but_isolated <- setdiff(mapped$gene, network_genes)
cat(sprintf("Mapped but isolated (no interactions): %d genes\n", length(mapped_but_isolated)))

# --- Save outputs ---
write_csv(edges, "results/network_edges.csv")
write_csv(nodes, "results/network_nodes.csv")
write_csv(mapped, "results/string_mapping.csv")

save(edges, nodes, mapped, network_genes, id_to_gene,
     evidence, consensus_2plus, consensus_3plus, core_genes,
     druggability,
     file = "results/02_network_built.RData")

cat("\n=== Files saved ===\n")
cat("  results/network_edges.csv\n")
cat("  results/network_nodes.csv\n")
cat("  results/string_mapping.csv\n")
cat("  results/02_network_built.RData\n")
cat("\nReady for 03_network_analysis.R\n")