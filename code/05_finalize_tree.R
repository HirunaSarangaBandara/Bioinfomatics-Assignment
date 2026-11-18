# ============================================================
# Purpose:
# Rebuild the Neighbor-Joining tree directly from the saved
# pairwise distance results, so that the tree and the reported
# pairwise table use exactly the same distances.
# ============================================================

if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}

library(ape)

pairwise_file <- file.path(
  "results", "tables", "pairwise_sequence_similarity.csv"
)
figure_dir <- file.path("results", "figures")
table_dir <- file.path("results", "tables")

if (!file.exists(pairwise_file)) {
  stop("pairwise_sequence_similarity.csv was not found. Run code/04_sequence_analysis.R first.")
}

pairwise <- read.csv(pairwise_file, stringsAsFactors = FALSE)
required <- c("organism_1", "organism_2", "p_distance")
missing <- setdiff(required, names(pairwise))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# Use a stable order for taxa in the distance matrix.
taxa <- c(
  "Lactiplantibacillus plantarum",
  "Lactobacillus acidophilus",
  "Limosilactobacillus reuteri",
  "Streptococcus thermophilus"
)

if (!setequal(unique(c(pairwise$organism_1, pairwise$organism_2)), taxa)) {
  stop("The pairwise file does not contain the four expected organisms.")
}

D <- matrix(
  0,
  nrow = length(taxa),
  ncol = length(taxa),
  dimnames = list(taxa, taxa)
)

for (i in seq_len(nrow(pairwise))) {
  a <- pairwise$organism_1[i]
  b <- pairwise$organism_2[i]
  d <- pairwise$p_distance[i]
  
  if (!is.finite(d) || d < 0 || d > 1) {
    stop("An invalid p-distance was found in the pairwise results.")
  }
  
  D[a, b] <- d
  D[b, a] <- d
}

if (anyNA(D) || any(D[upper.tri(D)] == 0)) {
  stop("The distance matrix is incomplete or contains an unexpected zero pairwise distance.")
}

# Four-point diagnostic for the four-taxon distance data.
# The three sums correspond to the three possible unrooted splits.
four_point <- data.frame(
  split = c(
    "(plantarum, acidophilus) | (reuteri, thermophilus)",
    "(plantarum, reuteri) | (acidophilus, thermophilus)",
    "(plantarum, thermophilus) | (acidophilus, reuteri)"
  ),
  distance_sum = c(
    D[1, 2] + D[3, 4],
    D[1, 3] + D[2, 4],
    D[1, 4] + D[2, 3]
  )
)

write.csv(
  four_point,
  file.path(table_dir, "four_point_tree_diagnostic.csv"),
  row.names = FALSE
)

cat("Four-point distance diagnostic:\n")
print(four_point)

# Rebuild the Neighbor-Joining tree from the same p-distances
# used in pairwise_sequence_similarity.csv.
tree <- ape::nj(stats::as.dist(D))

# Save the Newick tree.
ape::write.tree(
  tree,
  file = file.path(table_dir, "neighbor_joining_tree.nwk")
)

# Save a fresh tree image.
png(
  filename = file.path(figure_dir, "05_neighbor_joining_tree.png"),
  width = 1800,
  height = 1200,
  res = 200
)

plot(
  tree,
  type = "phylogram",
  use.edge.length = TRUE,
  cex = 0.85,
  main = "Neighbor-Joining Tree from Pairwise p-Distances"
)

ape::add.scale.bar()
dev.off()

cat("\nFinal tree rebuilt from the saved pairwise distance table.\n")
cat("Updated files:\n")
cat("- results/figures/05_neighbor_joining_tree.png\n")
cat("- results/tables/neighbor_joining_tree.nwk\n")
cat("- results/tables/four_point_tree_diagnostic.csv\n")