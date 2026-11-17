# ============================================================
# Script: Download and Import DNA Sequences
# Data source: NCBI GenBank
# ============================================================


# ------------------------------------------------------------
# 1. Install required packages if they are not available
# ------------------------------------------------------------

if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install(
    "Biostrings",
    ask = FALSE,
    update = FALSE
  )
}


# ------------------------------------------------------------
# 2. Load the required packages
# ------------------------------------------------------------

library(rentrez)
library(Biostrings)


# ------------------------------------------------------------
# 3. Check the project working directory
# ------------------------------------------------------------

if (!file.exists("README.md")) {
  stop(
    "README.md not found. Open the RStudio project ",
    "from the repository root directory."
  )
}


# ------------------------------------------------------------
# 4. Create the dataset information table
# ------------------------------------------------------------

dataset <- data.frame(
  
  organism = c(
    "Lactiplantibacillus plantarum",
    "Lactobacillus acidophilus",
    "Limosilactobacillus reuteri",
    "Streptococcus thermophilus"
  ),
  
  strain = c(
    "CIP 103151",
    "NBRC 13951",
    "DSM 20016",
    "ATCC 19258"
  ),
  
  accession = c(
    "NR_104573.1",
    "AB680529.1",
    "X76328.1",
    "AY188354.1"
  ),
  
  record_type = c(
    "Partial 16S rRNA (RefSeq)",
    "Partial 16S rRNA",
    "16S rRNA gene",
    "Complete 16S rRNA"
  ),
  
  stringsAsFactors = FALSE
)

# Create NCBI record links

dataset$ncbi_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/nuccore/",
  dataset$accession
)

# Display the dataset information

print(dataset)


# ------------------------------------------------------------
# 5. Download sequences from NCBI
# ------------------------------------------------------------

fasta_results <- vector(
  mode = "list",
  length = nrow(dataset)
)

for (i in seq_len(nrow(dataset))) {
  
  cat(
    "\nDownloading:",
    dataset$organism[i],
    "\n"
  )
  
  fasta_results[[i]] <- rentrez::entrez_fetch(
    
    db = "nuccore",
    
    id = dataset$accession[i],
    
    rettype = "fasta",
    
    retmode = "text"
    
  )
  
  # Avoid sending requests too frequently
  
  if (i < nrow(dataset)) {
    Sys.sleep(0.4)
  }
}


# Combine the downloaded FASTA records

fasta_text <- paste(
  unlist(fasta_results, use.names = FALSE),
  collapse = "\n"
)


# Check whether FASTA records were downloaded

if (!grepl(">", fasta_text, fixed = TRUE)) {
  
  stop(
    "No valid FASTA records were downloaded. ",
    "Please check your internet connection."
  )
  
}


# ------------------------------------------------------------
# 6. Save the downloaded sequences
# ------------------------------------------------------------

fasta_file <- file.path(
  "data",
  "raw",
  "selected_16s_sequences.fasta"
)

writeLines(
  fasta_text,
  con = fasta_file
)

cat("\nFASTA file saved successfully!\n")


# ------------------------------------------------------------
# 7. Import DNA sequences using Biostrings
# ------------------------------------------------------------

dna_sequences <- Biostrings::readDNAStringSet(
  fasta_file
)


# Check the number of downloaded sequences

if (length(dna_sequences) != nrow(dataset)) {
  
  stop(
    "Expected ",
    nrow(dataset),
    " sequences, but imported ",
    length(dna_sequences),
    ". Check the downloaded FASTA file."
  )
  
}


# ------------------------------------------------------------
# 8. Extract accession numbers from FASTA headers
# ------------------------------------------------------------

header_accessions <- sub(
  "[[:space:]].*$",
  "",
  names(dna_sequences)
)

sequence_index <- match(
  dataset$accession,
  header_accessions
)

if (anyNA(sequence_index)) {
  
  stop(
    "Some accession numbers could not be matched ",
    "with the imported FASTA records."
  )
  
}


# ------------------------------------------------------------
# 9. Create sequence metadata
# ------------------------------------------------------------

metadata <- dataset

metadata$sequence_length_bp <- as.integer(
  Biostrings::width(dna_sequences)[sequence_index]
)

metadata$fasta_header <- names(
  dna_sequences
)[sequence_index]


# ------------------------------------------------------------
# 10. Export metadata to CSV
# ------------------------------------------------------------

metadata_file <- file.path(
  "data",
  "raw",
  "sequence_metadata.csv"
)

write.csv(
  metadata,
  metadata_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Display final results
# ------------------------------------------------------------

cat("\n====================================\n")
cat("DATA IMPORT COMPLETED\n")
cat("====================================\n")

cat(
  "Number of sequences:",
  length(dna_sequences),
  "\n"
)

cat(
  "Total sequence bases:",
  sum(Biostrings::width(dna_sequences)),
  "\n"
)

cat("\nSequence metadata:\n")

print(
  metadata[, c(
    "organism",
    "accession",
    "sequence_length_bp"
  )]
)

cat("\nGenerated files:\n")

cat(fasta_file, "\n")
cat(metadata_file, "\n")

cat("\nAll tasks completed successfully!\n")