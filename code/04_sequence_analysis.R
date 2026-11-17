
# ============================================================
# Script: Main Bioinformatics Analysis
#
# Main methods:
# 1. Multiple Sequence Alignment
# 2. Conserved and Variable Site Analysis
# 3. Pairwise Sequence Similarity
# 4. Genetic Distance Calculation
# 5. Neighbor-Joining Phylogenetic Tree
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(Biostrings)
library(DECIPHER)
library(ape)
library(ggplot2)


# ------------------------------------------------------------
# 2. Define input and output directories
# ------------------------------------------------------------

fasta_file <- file.path(
  "data",
  "processed",
  "cleaned_16s_sequences.fasta"
)

metadata_file <- file.path(
  "data",
  "processed",
  "cleaned_sequence_metadata.csv"
)

figure_directory <- file.path(
  "results",
  "figures"
)

table_directory <- file.path(
  "results",
  "tables"
)

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Import the processed sequences and metadata
# ------------------------------------------------------------

if (!file.exists(fasta_file)) {
  stop("Processed FASTA file not found.")
}

if (!file.exists(metadata_file)) {
  stop("Processed metadata file not found.")
}

dna_sequences <- readDNAStringSet(
  fasta_file
)

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

if (length(dna_sequences) != nrow(metadata)) {
  stop("Number of sequences does not match metadata.")
}

cat("Processed DNA sequences imported successfully!\n")


# ------------------------------------------------------------
# 4. Match sequences with metadata
# ------------------------------------------------------------

fasta_headers <- names(dna_sequences)

fasta_accessions <- trimws(
  sub(
    ".*\\|[[:space:]]*",
    "",
    fasta_headers
  )
)

sequence_index <- match(
  metadata$accession,
  fasta_accessions
)

if (anyNA(sequence_index)) {
  stop("Unable to match all accession numbers.")
}

dna_sequences <- dna_sequences[
  sequence_index
]

# Use organism names as sequence identifiers

names(dna_sequences) <- metadata$organism

if (anyDuplicated(names(dna_sequences))) {
  stop("Duplicate organism names detected.")
}


# ------------------------------------------------------------
# 5. Perform Multiple Sequence Alignment
# ------------------------------------------------------------

cat("\n====================================\n")
cat("MULTIPLE SEQUENCE ALIGNMENT\n")
cat("====================================\n")

# Align the four DNA sequences

aligned_sequences <- DECIPHER::AlignSeqs(
  
  dna_sequences,
  
  iterations = 2,
  
  refinements = 1,
  
  processors = 1,
  
  verbose = TRUE
  
)

# Check whether the sequences have equal alignment lengths

aligned_lengths <- width(
  aligned_sequences
)

if (length(unique(aligned_lengths)) != 1) {
  stop("The aligned sequences do not have equal lengths.")
}

alignment_length <- aligned_lengths[1]

cat(
  "\nAlignment length:",
  alignment_length,
  "bp\n"
)


# ------------------------------------------------------------
# 6. Save the aligned DNA sequences
# ------------------------------------------------------------

aligned_fasta_file <- file.path(
  
  "data",
  
  "processed",
  
  "aligned_16s_sequences.fasta"
  
)

writeXStringSet(
  
  aligned_sequences,
  
  filepath = aligned_fasta_file,
  
  format = "fasta",
  
  width = 80
  
)

cat("Aligned FASTA file saved successfully!\n")


# ------------------------------------------------------------
# 7. Convert alignment into a character matrix
# ------------------------------------------------------------

alignment_matrix <- do.call(
  
  rbind,
  
  strsplit(
    as.character(aligned_sequences),
    split = "",
    fixed = TRUE
  )
  
)

rownames(alignment_matrix) <- metadata$organism

number_of_sequences <- nrow(
  alignment_matrix
)

alignment_length <- ncol(
  alignment_matrix
)


# ------------------------------------------------------------
# 8. Identify conserved and variable sites
# ------------------------------------------------------------

# Only positions containing A, T, G and C
# in every sequence are used for this analysis.

canonical_bases <- c(
  "A", "T", "G", "C"
)

fully_resolved_sites <- apply(
  
  alignment_matrix,
  
  2,
  
  function(site) {
    all(site %in% canonical_bases)
  }
  
)

resolved_matrix <- alignment_matrix[
  ,
  fully_resolved_sites,
  drop = FALSE
]

number_of_resolved_sites <- ncol(
  resolved_matrix
)


# Count fully conserved sites

if (number_of_resolved_sites > 0) {
  
  conserved_flags <- apply(
    
    resolved_matrix,
    
    2,
    
    function(site) {
      length(unique(site)) == 1
    }
    
  )
  
} else {
  
  conserved_flags <- logical(0)
  
}


conserved_sites <- sum(
  conserved_flags
)

variable_sites <- sum(
  !conserved_flags
)

incomplete_sites <- alignment_length -
  number_of_resolved_sites


# ------------------------------------------------------------
# 9. Generate alignment summary
# ------------------------------------------------------------

alignment_summary <- data.frame(
  
  number_of_sequences = number_of_sequences,
  
  alignment_length_bp = alignment_length,
  
  fully_resolved_sites = number_of_resolved_sites,
  
  conserved_sites = conserved_sites,
  
  variable_sites = variable_sites,
  
  columns_with_gaps_or_ambiguity = incomplete_sites
  
)

cat("\nALIGNMENT SUMMARY\n")

print(alignment_summary)


write.csv(
  
  alignment_summary,
  
  file.path(
    table_directory,
    "alignment_summary.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 10. Calculate pairwise sequence similarity
# ------------------------------------------------------------

cat("\nCalculating pairwise sequence similarities...\n")

if (number_of_sequences < 2) {
  stop("At least two sequences are required.")
}

pairwise_results <- data.frame()

distance_matrix <- matrix(
  
  0,
  
  nrow = number_of_sequences,
  
  ncol = number_of_sequences,
  
  dimnames = list(
    rownames(alignment_matrix),
    rownames(alignment_matrix)
  )
  
)

identity_matrix <- matrix(
  
  100,
  
  nrow = number_of_sequences,
  
  ncol = number_of_sequences,
  
  dimnames = list(
    rownames(alignment_matrix),
    rownames(alignment_matrix)
  )
  
)


# Compare every possible sequence pair

for (i in seq_len(number_of_sequences - 1)) {
  
  for (j in seq.int(i + 1, number_of_sequences)) {
    
    sequence_a <- alignment_matrix[i, ]
    
    sequence_b <- alignment_matrix[j, ]
    
    
    # Only compare canonical nucleotide positions
    
    valid_positions <-
      
      sequence_a %in% canonical_bases &
      
      sequence_b %in% canonical_bases
    
    
    comparable_sites <- sum(
      valid_positions
    )
    
    
    if (comparable_sites == 0) {
      stop("No comparable nucleotide positions found.")
    }
    
    
    matching_sites <- sum(
      
      sequence_a[valid_positions] ==
        
        sequence_b[valid_positions]
      
    )
    
    
    identity_percent <-
      
      (matching_sites / comparable_sites) * 100
    
    
    # Uncorrected p-distance
    
    p_distance <-
      
      1 - (matching_sites / comparable_sites)
    
    
    # Store results in the distance matrices
    
    distance_matrix[i, j] <- p_distance
    
    distance_matrix[j, i] <- p_distance
    
    identity_matrix[i, j] <- identity_percent
    
    identity_matrix[j, i] <- identity_percent
    
    
    # Save the pairwise comparison result
    
    pairwise_results <- rbind(
      
      pairwise_results,
      
      data.frame(
        
        organism_1 = rownames(alignment_matrix)[i],
        
        organism_2 = rownames(alignment_matrix)[j],
        
        comparable_sites = comparable_sites,
        
        matching_sites = matching_sites,
        
        percentage_identity = round(
          identity_percent,
          4
        ),
        
        p_distance = round(
          p_distance,
          6
        ),
        
        alignment_coverage_percent = round(
          comparable_sites / alignment_length * 100,
          2
        )
        
      )
      
    )
    
  }
  
}


# Display pairwise similarity results

cat("\nPAIRWISE SEQUENCE SIMILARITY\n")

print(pairwise_results)


# Save similarity table

write.csv(
  
  pairwise_results,
  
  file.path(
    table_directory,
    "pairwise_sequence_similarity.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 11. Export genetic distance matrix
# ------------------------------------------------------------

cat("\nGENETIC DISTANCE MATRIX\n")

print(
  round(distance_matrix, 5)
)

write.csv(
  
  as.data.frame(distance_matrix),
  
  file.path(
    table_directory,
    "genetic_distance_matrix.csv"
  ),
  
  row.names = TRUE
  
)


# ------------------------------------------------------------
# 12. Generate pairwise sequence identity chart
# ------------------------------------------------------------

pairwise_results$comparison <- paste(
  
  pairwise_results$organism_1,
  
  "vs",
  
  pairwise_results$organism_2
  
)


similarity_plot <- ggplot(
  
  pairwise_results,
  
  aes(
    
    x = reorder(
      comparison,
      percentage_identity
    ),
    
    y = percentage_identity
    
  )
  
) +
  
  geom_col() +
  
  geom_text(
    
    aes(
      label = paste0(
        round(percentage_identity, 2),
        "%"
      )
    ),
    
    hjust = -0.15
    
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    
    limits = c(0, 110),
    
    breaks = seq(0, 100, 20)
    
  ) +
  
  labs(
    
    title = "Pairwise 16S rRNA Sequence Similarity",
    
    x = "Sequence Comparison",
    
    y = "Sequence Identity (%)"
    
  ) +
  
  theme_minimal(base_size = 11)


print(similarity_plot)


ggsave(
  
  filename = file.path(
    
    figure_directory,
    
    "04_pairwise_sequence_similarity.png"
    
  ),
  
  plot = similarity_plot,
  
  width = 11,
  
  height = 7,
  
  dpi = 300
  
)


# ------------------------------------------------------------
# 13. Construct Neighbor-Joining phylogenetic tree
# ------------------------------------------------------------

cat("\n====================================\n")

cat("PHYLOGENETIC ANALYSIS\n")

cat("====================================\n")


# Convert the distance matrix into a distance object

genetic_distances <- stats::as.dist(
  distance_matrix
)


# Construct Neighbor-Joining tree

nj_tree <- ape::nj(
  genetic_distances
)


# Save the Newick tree

tree_file <- file.path(
  
  table_directory,
  
  "neighbor_joining_tree.nwk"
  
)

ape::write.tree(
  
  nj_tree,
  
  file = tree_file
  
)


cat("\nNeighbor-Joining tree constructed successfully!\n")


# ------------------------------------------------------------
# 14. Plot and save the phylogenetic tree
# ------------------------------------------------------------

tree_figure <- file.path(
  
  figure_directory,
  
  "05_neighbor_joining_tree.png"
  
)

png(
  
  filename = tree_figure,
  
  width = 1800,
  
  height = 1200,
  
  res = 200
  
)


plot(
  
  nj_tree,
  
  type = "phylogram",
  
  use.edge.length = TRUE,
  
  cex = 0.85,
  
  main = "Neighbor-Joining Tree of Selected Bacteria"
  
)

dev.off()


cat("Phylogenetic tree figure saved!\n")


# ------------------------------------------------------------
# 15. Final analysis summary
# ------------------------------------------------------------

cat("\n====================================\n")

cat("MAIN BIOINFORMATICS ANALYSIS COMPLETED\n")

cat("====================================\n")


cat(
  "\nNumber of sequences:",
  number_of_sequences,
  "\n"
)

cat(
  "Alignment length:",
  alignment_length,
  "bp\n"
)

cat(
  "Fully resolved sites:",
  number_of_resolved_sites,
  "\n"
)

cat(
  "Conserved sites:",
  conserved_sites,
  "\n"
)

cat(
  "Variable sites:",
  variable_sites,
  "\n"
)


cat("\nGenerated output files:\n")

cat("\nAligned sequence file:\n")

cat(aligned_fasta_file, "\n")


cat("\nResult tables:\n")

cat("alignment_summary.csv\n")

cat("pairwise_sequence_similarity.csv\n")

cat("genetic_distance_matrix.csv\n")

cat("neighbor_joining_tree.nwk\n")


cat("\nGenerated figures:\n")

cat("04_pairwise_sequence_similarity.png\n")

cat("05_neighbor_joining_tree.png\n")


cat("\nAll analysis tasks completed successfully!\n")