# Historical import-only prototype. For the complete reproducible analysis,
# source scripts/01_qc_clustering_annotation.R from the project root.
library(Seurat)
library(Matrix)
counts <- read.csv(
  "data/GSE123818_Root_single_cell_wt_datamatrix.csv.gz",
  row.names = 1,
  check.names = FALSE
)
dim(counts)