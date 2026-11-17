
# ============================================================
# Script: Data Cleaning and Preprocessing
# ============================================================


# ------------------------------------------------------------
# 1. Load required package
# ------------------------------------------------------------

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Biostrings package is not installed.")
}

library(Biostrings)


# ------------------------------------------------------------
# 2. Check required input files
# ------------------------------------------------------------

fasta_file <- file.path(
  "data",
  "raw",
  "selected_16s_sequences.fasta"
)

metadata_file <- file.path(
  "data",
  "raw",
  "sequence_metadata.csv"
)

if (!file.exists(fasta_file)) {
  stop("Raw FASTA file was not found.")
}

if (!file.exists(metadata_file)) {
  stop("Sequence metadata file was not found.")
}


# ------------------------------------------------------------
# 3. Import raw data
# ------------------------------------------------------------

raw_sequences <- readDNAStringSet(
  fasta_file
)

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

cat("Raw dataset imported successfully!\n")

cat(
  "Number of sequences:",
  length(raw_sequences),
  "\n"
)


# ------------------------------------------------------------
# 4. Validate dataset structure
# ------------------------------------------------------------

required_columns <- c(
  "organism",
  "strain",
  "accession"
)

missing_columns <- setdiff(
  required_columns,
  names(metadata)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (length(raw_sequences) != nrow(metadata)) {
  stop(
    "Number of FASTA records does not match metadata rows."
  )
}

if (anyDuplicated(metadata$accession)) {
  stop("Duplicate accession numbers detected in metadata.")
}


# ------------------------------------------------------------
# 5. Match FASTA sequences with their accession numbers
# ------------------------------------------------------------

headers <- names(raw_sequences)

matched_index <- vapply(
  metadata$accession,
  function(acc) {
    
    matches <- which(
      grepl(acc, headers, fixed = TRUE)
    )
    
    if (length(matches) == 1) {
      return(as.integer(matches))
    }
    
    return(NA_integer_)
  },
  integer(1)
)

if (anyNA(matched_index) ||
    anyDuplicated(matched_index)) {
  
  stop(
    "Could not match every metadata record ",
    "with exactly one FASTA sequence."
  )
}

raw_sequences <- raw_sequences[matched_index]


# ------------------------------------------------------------
# 6. Standardize sequence format
# ------------------------------------------------------------

sequence_text <- toupper(
  as.character(raw_sequences)
)

clean_sequences <- DNAStringSet(
  sequence_text
)

# Give each sequence a consistent descriptive name

names(clean_sequences) <- paste(
  metadata$organism,
  metadata$accession,
  sep = " | "
)

cat("Sequence formatting completed!\n")


# ------------------------------------------------------------
# 7. Calculate sequence lengths
# ------------------------------------------------------------

sequence_lengths <- as.integer(
  width(clean_sequences)
)


# ------------------------------------------------------------
# 8. Check nucleotide composition and ambiguity
# ------------------------------------------------------------

sequence_characters <- strsplit(
  as.character(clean_sequences),
  split = ""
)

valid_iupac <- c(
  "A", "C", "G", "T",
  "M", "R", "W", "S", "Y",
  "K", "V", "H", "D", "B", "N"
)

ambiguous_count <- vapply(
  sequence_characters,
  function(x) {
    sum(!x %in% c("A", "C", "G", "T"))
  },
  integer(1)
)

ambiguous_symbols <- vapply(
  sequence_characters,
  function(x) {
    
    symbols <- sort(
      unique(
        x[!x %in% c("A", "C", "G", "T")]
      )
    )
    
    if (length(symbols) == 0) {
      return("None")
    }
    
    paste(symbols, collapse = ", ")
  },
  character(1)
)

unexpected_symbols <- vapply(
  sequence_characters,
  function(x) {
    
    symbols <- sort(
      unique(x[!x %in% valid_iupac])
    )
    
    if (length(symbols) == 0) {
      return("None")
    }
    
    paste(symbols, collapse = ", ")
  },
  character(1)
)


# ------------------------------------------------------------
# 9. Detect duplicate sequences
# ------------------------------------------------------------

sequence_values <- as.character(
  clean_sequences
)

duplicate_sequences <-
  duplicated(sequence_values) |
  duplicated(sequence_values, fromLast = TRUE)


# ------------------------------------------------------------
# 10. Prepare quality control summary
# ------------------------------------------------------------

qc_status <- vapply(
  seq_along(clean_sequences),
  function(i) {
    
    issues <- character()
    
    if (ambiguous_count[i] > 0) {
      issues <- c(
        issues,
        "Ambiguous nucleotide codes present"
      )
    }
    
    if (duplicate_sequences[i]) {
      issues <- c(
        issues,
        "Exact duplicate sequence"
      )
    }
    
    if (unexpected_symbols[i] != "None") {
      issues <- c(
        issues,
        "Unexpected symbols detected"
      )
    }
    
    if (length(issues) == 0) {
      return("Pass")
    }
    
    paste(
      issues,
      collapse = "; "
    )
  },
  character(1)
)


# ------------------------------------------------------------
# 11. Create QC report
# ------------------------------------------------------------

qc_report <- data.frame(
  
  organism = metadata$organism,
  
  strain = metadata$strain,
  
  accession = metadata$accession,
  
  original_sequence_length = sequence_lengths,
  
  ambiguous_nucleotide_count = ambiguous_count,
  
  ambiguous_symbols = ambiguous_symbols,
  
  unexpected_symbols = unexpected_symbols,
  
  exact_duplicate = duplicate_sequences,
  
  qc_status = qc_status,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 12. Display quality control results
# ------------------------------------------------------------

cat("\n====================================\n")
cat("SEQUENCE QUALITY CONTROL RESULTS\n")
cat("====================================\n")

print(qc_report)


# ------------------------------------------------------------
# 13. Export processed FASTA file
# ------------------------------------------------------------

processed_fasta <- file.path(
  "data",
  "processed",
  "cleaned_16s_sequences.fasta"
)

writeXStringSet(
  clean_sequences,
  filepath = processed_fasta,
  format = "fasta",
  width = 80
)


# ------------------------------------------------------------
# 14. Export QC summary
# ------------------------------------------------------------

qc_file <- file.path(
  "data",
  "processed",
  "sequence_qc_summary.csv"
)

write.csv(
  qc_report,
  qc_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Export cleaned metadata
# ------------------------------------------------------------

cleaned_metadata <- metadata

cleaned_metadata$sequence_length_bp <-
  sequence_lengths

cleaned_metadata$ambiguous_nucleotide_count <-
  ambiguous_count

cleaned_metadata$qc_status <-
  qc_status

metadata_output <- file.path(
  "data",
  "processed",
  "cleaned_sequence_metadata.csv"
)

write.csv(
  cleaned_metadata,
  metadata_output,
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Final summary
# ------------------------------------------------------------

cat("\n====================================\n")
cat("DATA CLEANING COMPLETED\n")
cat("====================================\n")

cat(
  "Total sequences:",
  length(clean_sequences),
  "\n"
)

cat(
  "Exact duplicate sequences:",
  sum(duplicate_sequences) / 2,
  "\n"
)

cat(
  "Sequences containing ambiguous codes:",
  sum(ambiguous_count > 0),
  "\n"
)

cat("\nGenerated files:\n")

cat(processed_fasta, "\n")
cat(qc_file, "\n")
cat(metadata_output, "\n")

cat("\nAll sequences were preserved.\n")
