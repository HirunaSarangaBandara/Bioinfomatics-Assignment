# ============================================================
# Script: Exploratory Data Analysis (EDA)
# ============================================================


# ------------------------------------------------------------
# 1. Install and load required packages
# ------------------------------------------------------------

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Please install Biostrings before running this script.")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(Biostrings)
library(ggplot2)


# ------------------------------------------------------------
# 2. Define input and output paths
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


# Create output directories if they do not exist

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
# 3. Check input files
# ------------------------------------------------------------

if (!file.exists(fasta_file)) {
  stop("Processed FASTA file not found.")
}

if (!file.exists(metadata_file)) {
  stop("Processed metadata file not found.")
}


# ------------------------------------------------------------
# 4. Import processed sequences and metadata
# ------------------------------------------------------------

dna_sequences <- readDNAStringSet(
  fasta_file
)

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

cat("Data imported successfully!\n")

cat(
  "Number of sequences:",
  length(dna_sequences),
  "\n"
)


# ------------------------------------------------------------
# 5. Match sequences with metadata
# ------------------------------------------------------------

fasta_headers <- names(dna_sequences)

if (is.null(fasta_headers) ||
    any(is.na(fasta_headers))) {
  stop("FASTA sequence headers are missing.")
}

# Extract accession numbers from sequence headers

fasta_accessions <- trimws(
  sub(
    ".*\\|[[:space:]]*",
    "",
    fasta_headers
  )
)

matched_index <- match(
  metadata$accession,
  fasta_accessions
)

if (anyNA(matched_index)) {
  stop("Unable to match all sequence accessions.")
}

if (length(matched_index) != length(dna_sequences)) {
  stop("Number of sequences and metadata records differ.")
}

# Arrange sequences in metadata order

dna_sequences <- dna_sequences[matched_index]

if (any(width(dna_sequences) == 0)) {
  stop("An empty DNA sequence was detected.")
}


# ------------------------------------------------------------
# 6. Calculate sequence lengths
# ------------------------------------------------------------

sequence_lengths <- as.integer(
  width(dna_sequences)
)

cat("\nSequence lengths calculated.\n")

print(sequence_lengths)


# ------------------------------------------------------------
# 7. Extract nucleotide characters
# ------------------------------------------------------------

sequence_characters <- strsplit(
  as.character(dna_sequences),
  split = ""
)


# ------------------------------------------------------------
# 8. Count A, T, G and C nucleotides
# ------------------------------------------------------------

count_base <- function(sequence, base) {
  
  as.integer(
    sum(sequence == base)
  )
  
}

A_count <- vapply(
  sequence_characters,
  count_base,
  integer(1),
  base = "A"
)

T_count <- vapply(
  sequence_characters,
  count_base,
  integer(1),
  base = "T"
)

G_count <- vapply(
  sequence_characters,
  count_base,
  integer(1),
  base = "G"
)

C_count <- vapply(
  sequence_characters,
  count_base,
  integer(1),
  base = "C"
)


# ------------------------------------------------------------
# 9. Calculate total canonical nucleotides
# ------------------------------------------------------------

canonical_total <-
  A_count +
  T_count +
  G_count +
  C_count


# Count ambiguous nucleotide characters

ambiguous_count <-
  sequence_lengths -
  canonical_total


# ------------------------------------------------------------
# 10. Calculate GC content
# ------------------------------------------------------------

gc_content <- ifelse(
  canonical_total > 0,
  ((G_count + C_count) / canonical_total) * 100,
  NA_real_
)


# ------------------------------------------------------------
# 11. Calculate nucleotide composition percentages
# ------------------------------------------------------------

composition_counts <- data.frame(
  
  A = A_count,
  
  T = T_count,
  
  G = G_count,
  
  C = C_count
  
)


composition_percentages <- data.frame(
  
  A = ifelse(
    canonical_total > 0,
    A_count / canonical_total * 100,
    NA_real_
  ),
  
  T = ifelse(
    canonical_total > 0,
    T_count / canonical_total * 100,
    NA_real_
  ),
  
  G = ifelse(
    canonical_total > 0,
    G_count / canonical_total * 100,
    NA_real_
  ),
  
  C = ifelse(
    canonical_total > 0,
    C_count / canonical_total * 100,
    NA_real_
  )
  
)


# ------------------------------------------------------------
# 12. Create the sequence summary table
# ------------------------------------------------------------

summary_table <- data.frame(
  
  organism = metadata$organism,
  
  strain = metadata$strain,
  
  accession = metadata$accession,
  
  sequence_length_bp = sequence_lengths,
  
  A_count = A_count,
  
  T_count = T_count,
  
  G_count = G_count,
  
  C_count = C_count,
  
  ambiguous_count = ambiguous_count,
  
  gc_content_percent = round(gc_content, 2),
  
  stringsAsFactors = FALSE
  
)


# Display summary

cat("\n====================================\n")
cat("SEQUENCE SUMMARY\n")
cat("====================================\n")

print(summary_table)


# ------------------------------------------------------------
# 13. Save sequence summary table
# ------------------------------------------------------------

write.csv(
  
  summary_table,
  
  file.path(
    table_directory,
    "sequence_summary.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 14. Create nucleotide composition table
# ------------------------------------------------------------

composition_long <- data.frame(
  
  organism = rep(
    metadata$organism,
    each = 4
  ),
  
  nucleotide = rep(
    c("A", "T", "G", "C"),
    times = nrow(metadata)
  ),
  
  count = as.vector(
    t(as.matrix(composition_counts))
  ),
  
  percentage = as.vector(
    t(as.matrix(composition_percentages))
  )
  
)


# Save nucleotide composition table

write.csv(
  
  composition_long,
  
  file.path(
    table_directory,
    "nucleotide_composition.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 15. Figure 1 - Sequence Length Comparison
# ------------------------------------------------------------

length_plot <- ggplot(
  
  summary_table,
  
  aes(
    x = reorder(organism, sequence_length_bp),
    y = sequence_length_bp
  )
  
) +
  
  geom_col() +
  
  geom_text(
    aes(label = sequence_length_bp),
    hjust = -0.2
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18))
  ) +
  
  labs(
    
    title = "16S rRNA Sequence Length Comparison",
    
    subtitle = "Length of selected bacterial sequences",
    
    x = "Bacterial Species",
    
    y = "Sequence Length (bp)"
    
  ) +
  
  theme_minimal(base_size = 12)


# Display plot

print(length_plot)


# Save figure

ggsave(
  
  filename = file.path(
    figure_directory,
    "01_sequence_length_comparison.png"
  ),
  
  plot = length_plot,
  
  width = 10,
  
  height = 6,
  
  dpi = 300
  
)


# ------------------------------------------------------------
# 16. Figure 2 - GC Content Comparison
# ------------------------------------------------------------

gc_plot <- ggplot(
  
  summary_table,
  
  aes(
    x = reorder(organism, gc_content_percent),
    y = gc_content_percent
  )
  
) +
  
  geom_col() +
  
  geom_text(
    aes(
      label = paste0(
        round(gc_content_percent, 2),
        "%"
      )
    ),
    vjust = -0.5
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.08))
  ) +
  
  labs(
    
    title = "GC Content of Selected Bacterial Sequences",
    
    subtitle = "Percentage of G and C nucleotides",
    
    x = "Bacterial Species",
    
    y = "GC Content (%)"
    
  ) +
  
  theme_minimal(base_size = 12)


# Display plot

print(gc_plot)


# Save figure

ggsave(
  
  filename = file.path(
    figure_directory,
    "02_gc_content_comparison.png"
  ),
  
  plot = gc_plot,
  
  width = 10,
  
  height = 6,
  
  dpi = 300
  
)


# ------------------------------------------------------------
# 17. Figure 3 - Nucleotide Composition
# ------------------------------------------------------------

composition_plot <- ggplot(
  
  composition_long,
  
  aes(
    x = organism,
    y = percentage,
    fill = nucleotide
  )
  
) +
  
  geom_col() +
  
  coord_flip() +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  
  labs(
    
    title = "Nucleotide Composition of 16S rRNA Sequences",
    
    subtitle = "Percentage of A, T, G and C nucleotides",
    
    x = "Bacterial Species",
    
    y = "Nucleotide Percentage (%)",
    
    fill = "Nucleotide"
    
  ) +
  
  theme_minimal(base_size = 12)


# Display plot

print(composition_plot)


# Save figure

ggsave(
  
  filename = file.path(
    figure_directory,
    "03_nucleotide_composition.png"
  ),
  
  plot = composition_plot,
  
  width = 10,
  
  height = 6,
  
  dpi = 300
  
)


# ------------------------------------------------------------
# 18. Export final analysis summary
# ------------------------------------------------------------

cat("\n====================================\n")
cat("EDA COMPLETED SUCCESSFULLY\n")
cat("====================================\n")

cat("\nNumber of sequences analyzed:")

cat(length(dna_sequences), "\n")

cat("\nAverage sequence length:")

cat(
  round(mean(sequence_lengths), 2),
  "bp\n"
)

cat("\nAverage GC content:")

cat(
  round(mean(gc_content, na.rm = TRUE), 2),
  "%\n"
)

cat("\nGenerated figures:\n")

cat("01_sequence_length_comparison.png\n")

cat("02_gc_content_comparison.png\n")

cat("03_nucleotide_composition.png\n")

cat("\nGenerated tables:\n")

cat("sequence_summary.csv\n")

cat("nucleotide_composition.csv\n")

cat("\nAll EDA tasks completed!\n")