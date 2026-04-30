# =============================================================================
# 03_network_analysis.R
# Step 3: Network topology analysis — centrality, communities, small-world test
# =============================================================================

library(tidyverse)
library(igraph)

# --- Load previous workspace ---
load("results/02_network_built.RData")

# --- Build igraph object ---
g <- graph_from_data_frame(
  edges %>% select(gene1, gene2),
  directed = FALSE,
  vertices = nodes
)

# Remove self-loops and multiple edges
g <- simplify(g)

cat("=== Raw Network ===\n")
cat(sprintf("Nodes: %d\n", vcount(g)))
cat(sprintf("Edges: %d\n", ecount(g)))

# --- Extract largest connected component ---
components <- components(g)
largest_cc <- which.max(components$csize)
g_main <- induced_subgraph(g, which(components$membership == largest_cc))

cat(sprintf("\n=== Largest Connected Component ===\n"))
cat(sprintf("Nodes: %d (%.1f%% of total)\n", vcount(g_main),
            100 * vcount(g_main) / vcount(g)))
cat(sprintf("Edges: %d\n", ecount(g_main)))
cat(sprintf("Isolated components removed: %d\n", components$no - 1))

# Work with the main component from here
g <- g_main

# --- Compute centrality metrics ---
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, normalized = TRUE)
V(g)$closeness <- closeness(g, normalized = TRUE)
V(g)$eigenvector <- eigen_centrality(g)$vector

# --- Global network metrics ---
cat(sprintf("\n=== Global Network Metrics ===\n"))
cat(sprintf("Density: %.4f\n", edge_density(g)))
cat(sprintf("Average path length: %.3f\n", mean_distance(g)))
cat(sprintf("Clustering coefficient: %.4f\n", transitivity(g, type = "global")))
cat(sprintf("Diameter: %d\n", diameter(g)))
cat(sprintf("Average degree: %.2f\n", mean(degree(g))))

# --- Community detection (Louvain) ---
set.seed(42)
communities <- cluster_louvain(g)
V(g)$community <- membership(communities)

cat(sprintf("\n=== Community Detection (Louvain) ===\n"))
cat(sprintf("Number of communities: %d\n", length(communities)))
cat(sprintf("Modularity: %.4f\n", modularity(communities)))

# Show community sizes
comm_sizes <- table(V(g)$community)
comm_sizes_sorted <- sort(comm_sizes, decreasing = TRUE)
cat("\nCommunity sizes:\n")
for (i in seq_along(comm_sizes_sorted)) {
  cat(sprintf("  Module %s: %d genes\n",
              names(comm_sizes_sorted)[i], comm_sizes_sorted[i]))
}

# --- Hub genes (top 20 by degree) ---
centrality_df <- data.frame(
  gene = V(g)$name,
  degree = V(g)$degree,
  betweenness = V(g)$betweenness,
  closeness = V(g)$closeness,
  eigenvector = V(g)$eigenvector,
  community = V(g)$community,
  evidence_count = V(g)$evidence_count,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(degree))

cat("\n=== Top 20 Hub Genes (by degree) ===\n")
cat(sprintf("%-15s %6s %12s %10s %10s %5s %5s\n",
            "Gene", "Degree", "Betweenness", "Closeness", "Eigenvec", "Comm", "Evid"))
for (i in 1:min(20, nrow(centrality_df))) {
  row <- centrality_df[i, ]
  cat(sprintf("%-15s %6d %12.4f %10.4f %10.4f %5d %5d\n",
              row$gene, row$degree, row$betweenness,
              row$closeness, row$eigenvector,
              row$community, row$evidence_count))
}

# --- Bridge genes (high betweenness but moderate degree) ---
# These connect different modules — often the most interesting targets
centrality_df <- centrality_df %>%
  mutate(
    degree_rank = rank(-degree),
    betweenness_rank = rank(-betweenness),
    bridge_score = degree_rank - betweenness_rank
    # Negative bridge_score = higher betweenness than degree would predict = bridge gene
  )

bridges <- centrality_df %>%
  filter(betweenness_rank <= 30 & degree_rank > 15) %>%
  arrange(betweenness_rank)

if (nrow(bridges) > 0) {
  cat("\n=== Bridge Genes (high betweenness, moderate degree) ===\n")
  for (i in 1:min(10, nrow(bridges))) {
    row <- bridges[i, ]
    cat(sprintf("  %s: degree=%d (rank %d), betweenness=%.4f (rank %d), community=%d\n",
                row$gene, row$degree, as.integer(row$degree_rank),
                row$betweenness, as.integer(row$betweenness_rank),
                row$community))
  }
}

# --- Small-world test: compare to random networks ---
cat("\n=== Small-World Network Test ===\n")
cat("Generating 1000 random networks for comparison...\n")

n_random <- 1000
random_clustering <- numeric(n_random)
random_path_length <- numeric(n_random)

set.seed(42)
for (i in 1:n_random) {
  g_random <- sample_gnm(vcount(g), ecount(g))
  random_clustering[i] <- transitivity(g_random, type = "global")
  if (is_connected(g_random)) {
    random_path_length[i] <- mean_distance(g_random)
  } else {
    random_path_length[i] <- NA
  }
}

real_clustering <- transitivity(g, type = "global")
real_path <- mean_distance(g)
mean_random_clustering <- mean(random_clustering, na.rm = TRUE)
mean_random_path <- mean(random_path_length, na.rm = TRUE)

# Small-world coefficient: sigma = (C/C_rand) / (L/L_rand)
# sigma > 1 indicates small-world properties
sigma <- (real_clustering / mean_random_clustering) / (real_path / mean_random_path)

cat(sprintf("\n  Real clustering coefficient: %.4f\n", real_clustering))
cat(sprintf("  Random clustering (mean):    %.4f\n", mean_random_clustering))
cat(sprintf("  Ratio C/C_rand:              %.2f\n", real_clustering / mean_random_clustering))
cat(sprintf("\n  Real avg path length:         %.3f\n", real_path))
cat(sprintf("  Random path length (mean):    %.3f\n", mean_random_path))
cat(sprintf("  Ratio L/L_rand:               %.2f\n", real_path / mean_random_path))
cat(sprintf("\n  Small-world coefficient (σ):  %.2f\n", sigma))

if (sigma > 1) {
  cat("  → Network exhibits SMALL-WORLD properties (σ > 1)\n")
  cat("  → High clustering + short path lengths = biologically meaningful structure\n")
} else {
  cat("  → Network does not show clear small-world properties\n")
}

# --- Genes per module and their evidence counts ---
cat("\n=== Module Composition Summary ===\n")
for (comm_id in names(comm_sizes_sorted)) {
  comm_genes <- centrality_df %>% filter(community == as.integer(comm_id))
  multi_evidence <- comm_genes %>% filter(evidence_count >= 2)
  top_hub <- comm_genes %>% slice(1)
  cat(sprintf("\nModule %s (%d genes):\n", comm_id, nrow(comm_genes)))
  cat(sprintf("  Genes with ≥2 evidence layers: %d\n", nrow(multi_evidence)))
  cat(sprintf("  Top hub: %s (degree=%d, evidence=%d/5)\n",
              top_hub$gene, top_hub$degree, top_hub$evidence_count))
  # List core genes (≥3 layers) in this module
  core_in_module <- comm_genes %>% filter(evidence_count >= 3)
  if (nrow(core_in_module) > 0) {
    cat(sprintf("  Core genes: %s\n", paste(core_in_module$gene, collapse=", ")))
  }
}

# --- Save results ---
write_csv(centrality_df, "results/network_centrality.csv")

save(g, centrality_df, communities, bridges,
     real_clustering, mean_random_clustering,
     real_path, mean_random_path, sigma,
     random_clustering, random_path_length,
     evidence, consensus_2plus, consensus_3plus, core_genes,
     druggability, edges, nodes, network_genes,
     file = "results/03_network_analyzed.RData")

cat("\n=== Files saved ===\n")
cat("  results/network_centrality.csv\n")
cat("  results/03_network_analyzed.RData\n")
cat("\nReady for 04_enrichment.R\n")