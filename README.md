# Comparative 16S rRNA Sequence Analysis of Selected Lactic Acid Bacteria

**Student Registration Number:** 2021/ICT/81

**Project Type:** Bioinformatics – ICA3

**Programming Language:** R

## 1. Project Overview

This project aims to investigate the sequence similarities and differences among selected lactic acid bacteria using publicly available 16S ribosomal RNA (16S rRNA) gene sequences. The analysis will be conducted using the R programming language and appropriate bioinformatics packages.

The project will examine nucleotide composition, sequence length, GC content, sequence similarity and evolutionary relationships among the selected bacterial species.

## 2. Problem Identification

Bacterial species can share similar genetic sequences while also containing nucleotide variations that help distinguish them. The 16S rRNA gene is commonly used in bacterial identification and comparative sequence analysis because it contains both conserved and variable regions.

However, understanding the genetic similarities and differences among selected bacterial species requires appropriate sequence processing, comparison and interpretation.

This project will address this problem by applying bioinformatics techniques to publicly available 16S rRNA sequences from selected lactic acid bacteria.

## 3. Research Question

What similarities and differences can be identified among the 16S rRNA gene sequences of selected lactic acid bacteria through nucleotide composition analysis, sequence comparison and phylogenetic analysis?

## 4. Research Aim

To investigate the genetic similarities and differences among selected lactic acid bacteria using comparative 16S rRNA sequence analysis in R.

## 5. Research Objectives

1. Obtain publicly available 16S rRNA gene sequences of selected lactic acid bacteria.
2. Import and preprocess the collected DNA sequences using R.
3. Analyze sequence length, nucleotide composition and GC content.
4. Compare the similarities and differences among the selected sequences.
5. Construct and interpret a phylogenetic tree to explore their evolutionary relationships.

## 6. Dataset and Data Source

The dataset for this research consists of publicly available 16S rRNA gene sequences obtained from the NCBI GenBank database.

Four bacterial species were selected for comparative sequence analysis. The selected records represent different lactic acid bacterial groups and have identifiable accession numbers.

| Organism | Strain | Accession Number |
|---|---|---|
| Lactiplantibacillus plantarum | CIP 103151 | NR_104573.1 |
| Lactobacillus acidophilus | NBRC 13951 | AB680529.1 |
| Limosilactobacillus reuteri | DSM 20016 | X76328.1 |
| Streptococcus thermophilus | ATCC 19258 | AY188354.1 |

### Data Collection Method

The nucleotide records were retrieved from NCBI using the Entrez programming interface through the `rentrez` R package.

The sequences were imported using the `Biostrings` package and stored in FASTA format.

A separate CSV file was generated to document the organism names, strain information, accession numbers and sequence lengths.

### Dataset Source

NCBI Nucleotide Database: https://www.ncbi.nlm.nih.gov/nuccore/

Individual sequence records:

1. https://www.ncbi.nlm.nih.gov/nuccore/NR_104573.1
2. https://www.ncbi.nlm.nih.gov/nuccore/AB680529.1
3. https://www.ncbi.nlm.nih.gov/nuccore/X76328.1
4. https://www.ncbi.nlm.nih.gov/nuccore/AY188354.1

### Dataset Limitations

The dataset contains one representative sequence for each selected bacterial species. The selected records also differ in their reported sequence completeness.

Therefore, comparisons will be interpreted at the level of the available sequence records, and the study will not make population-level conclusions about the bacterial species.

## 7. Data Cleaning and Preprocessing

The downloaded 16S rRNA sequences were imported into R using the Biostrings package. Before performing the main analysis, a preprocessing stage was implemented to examine the structure and quality of the collected sequences.

The FASTA records were matched with their corresponding accession numbers in the metadata table. The sequences were then standardized using uppercase nucleotide characters and consistent sequence headers.

The following quality-control checks were performed:

1. **Sequence validation:** Checked whether the FASTA records could be successfully imported and matched with the corresponding metadata.
2. **Sequence length checking:** Calculated the length of each sequence in base pairs.
3. **Nucleotide ambiguity checking:** Identified ambiguous nucleotide codes such as N, R and Y.
4. **Duplicate detection:** Checked whether any sequences were completely identical.
5. **Sequence standardization:** Standardized sequence formatting and naming.

The original raw sequences were preserved, while the standardized sequences were exported into a separate processed-data directory.

A quality-control summary and an updated metadata file were also generated in CSV format.

Sequences containing ambiguous nucleotide codes were retained because these characters can represent valid IUPAC nucleotide ambiguity codes. No records were automatically removed during this preprocessing stage.

The processed sequences will be used for exploratory data analysis and subsequent sequence comparison.

## 8. Exploratory Data Analysis

Exploratory Data Analysis (EDA) was conducted to examine the characteristics of the four selected bacterial 16S rRNA sequences.

The analysis was performed using R, with the ggplot2 package used to produce graphical representations of the results.

### 8.1 Sequence Length Analysis

The length of each DNA sequence was calculated in base pairs. A bar chart was generated to compare the sequence lengths of the selected bacterial records.

This analysis provides an initial understanding of the differences in sequence coverage among the collected records.

### 8.2 Nucleotide Composition Analysis

The numbers of adenine (A), thymine (T), guanine (G) and cytosine (C) nucleotides were calculated for each sequence.

The nucleotide counts were converted into percentages based on the total number of unambiguous nucleotides (A, T, G and C).

A stacked bar chart was created to compare the nucleotide composition of the selected records.

### 8.3 GC Content Analysis

GC content was calculated using the following formula:

GC Content (%) = ((G + C) / (A + T + G + C)) × 100

The calculated percentages were presented using a comparative bar chart.

Ambiguous nucleotide characters were excluded from the GC content denominator to avoid treating unknown nucleotide identities as G, C, A or T.

### 8.4 Visualization and Export

Three visualizations were generated:

1. 16S rRNA sequence length comparison.
2. GC content comparison.
3. Nucleotide composition comparison.

Two summary tables were also generated and exported as CSV files in the results/tables directory.

### 8.5 Interpretation and Limitations

The EDA results provide an overview of sequence length, nucleotide composition and GC content among the selected bacterial records.

The sequence length comparisons must be interpreted carefully because the selected NCBI records have different levels of sequence completeness.

Furthermore, only one representative sequence was selected for each bacterial species. Therefore, the descriptive comparisons cannot establish the full range of genetic variation within each species.

The actual numerical observations and their biological interpretations will be discussed based on the generated tables and figures.

## 9. Scope and Limitations

This study will focus on a selected group of lactic acid bacteria and their publicly available 16S rRNA gene sequences. It will not examine complete bacterial genomes or experimental laboratory samples.

The findings will depend on the quality, availability and comparability of the selected sequence records. In addition, 16S rRNA sequence analysis may not distinguish every closely related bacterial species with certainty.

## 10. Reproducibility

All data processing and analysis will be performed using R. Source code, dataset information, generated figures and analytical results will be organized within this GitHub repository.

The final project will include the required instructions, package information and report so that the analysis can be reproduced.

---

**Project Status:** Initial research problem and project scope defined. Dataset collection and analysis will be completed in subsequent development stages.