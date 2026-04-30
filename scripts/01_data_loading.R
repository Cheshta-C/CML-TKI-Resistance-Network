# =============================================================================
# 01_data_loading.R
# Multi-Evidence Network Analysis of CML TKI Resistance
# Step 1: Load and validate all evidence layers
# =============================================================================

library(tidyverse)

# --- Load the evidence matrix ---
evidence <- read_csv("data/evidence_matrix.csv")

# Rename the 'gene' column if needed
colnames(evidence)[1] <- "gene"

cat("=== Evidence Matrix ===\n")
cat(sprintf("Total genes: %d\n", nrow(evidence)))
cat(sprintf("Layers: %s\n", paste(colnames(evidence)[2:6], collapse=", ")))

# --- Evidence distribution ---
cat("\n=== Evidence Distribution ===\n")
for (i in 5:1) {
  n <- sum(evidence$evidence_count == i)
  cat(sprintf("  %d/5 layers: %d genes\n", i, n))
}

# --- Load consensus gene lists ---
consensus_2plus <- read_lines("data/consensus_genes_2plus.txt")
consensus_3plus <- read_lines("data/consensus_genes_3plus.txt")

cat(sprintf("\nConsensus genes (≥2 layers): %d\n", length(consensus_2plus)))
cat(sprintf("Core genes (≥3 layers): %d\n", length(consensus_3plus)))

# --- Load individual layer data for later annotation ---
layer1_full <- read_csv("data/krishnan_2023_layer1_full.csv")
layer2_full <- read_csv("data/awad_2024_layer2_full.csv")
layer3_full <- read_csv("data/krishnan_2023_layer3_classifier_genes.csv")
layer4_full <- read_csv("data/sacco_2025_layer4_patient_degs.csv")
layer5_full <- read_csv("data/sacco_2025_layer5_full.csv")
druggability <- read_csv("data/sacco_2025_druggability.csv")

cat("\n=== Layer data loaded ===\n")
cat(sprintf("  Layer 1 (Krishnan scRNA-seq): %d entries\n", nrow(layer1_full)))
cat(sprintf("  Layer 2 (Awad CRISPR KO): %d entries\n", nrow(layer2_full)))
cat(sprintf("  Layer 3 (Krishnan ML): %d entries\n", nrow(layer3_full)))
cat(sprintf("  Layer 4 (Sacco RNA-seq): %d entries\n", nrow(layer4_full)))
cat(sprintf("  Layer 5 (Sacco Phospho): %d entries\n", nrow(layer5_full)))
cat(sprintf("  Druggability scores: %d entries\n", nrow(druggability)))

# --- Show the 14 core genes ---
core_genes <- evidence %>%
  filter(evidence_count >= 3) %>%
  arrange(desc(evidence_count))

cat("\n=== Core Resistance Genes (≥3 layers) ===\n")
for (i in 1:nrow(core_genes)) {
  row <- core_genes[i, ]
  layers_hit <- colnames(evidence)[2:6][as.logical(row[2:6])]
  cat(sprintf("  %s [%d/5]: %s\n", 
              row$gene, row$evidence_count, 
              paste(layers_hit, collapse=", ")))
}

# --- Save workspace for next script ---
save(evidence, consensus_2plus, consensus_3plus, core_genes,
     layer1_full, layer2_full, layer3_full, layer4_full, layer5_full,
     druggability,
     file = "results/01_data_loaded.RData")

cat("\nWorkspace saved to results/01_data_loaded.RData\n")
cat("Ready for 02_string_network.R\n")