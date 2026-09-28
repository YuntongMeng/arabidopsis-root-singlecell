#!/usr/bin/env Rscript

# ============================================================
# Arabidopsis root single-cell RNA-seq
# Recovery of Denyer et al. (2019) WT replicate identities
#
# Main purpose:
#   Recover rep1/rep2 labels for the 4,727 cells in the published
#   GSE123818 WT expression matrix without re-running the raw pipeline.
#
# Strategy:
#   Use Shahan et al. (2022) Data S3 as a cell-level reference, validate
#   the barcode-suffix crosswalk on direct overlaps, and retain separate
#   evidence labels for direct matches and suffix-based inference.
#
# Output:
#   data/reference/Denyer2019_WT_barcode_to_replicate.tsv
# ============================================================


# ------------------------------------------------------------
# 1. Check the required package
# ------------------------------------------------------------

# Only readxl is needed because this script reads two columns from the
# Shahan supplementary workbook. All other processing uses base R.

required_packages <- c("readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

# ------------------------------------------------------------
# 2. Define and validate project paths
# ------------------------------------------------------------

# The script is intended to run from the project root. Keeping the input and
# output paths explicit makes the provenance of the recovered mapping clear.

project_dir <- normalizePath(".", mustWork = TRUE)
matrix_path <- file.path(
  project_dir,
  "data",
  "GSE123818_Root_single_cell_wt_datamatrix.csv.gz"
)
metadata_path <- file.path(
  project_dir,
  "data",
  "reference",
  "Shahan2022_DataS3_WT_atlas_metadata.xltx"
)
output_path <- file.path(
  project_dir,
  "data",
  "reference",
  "Denyer2019_WT_barcode_to_replicate.tsv"
)

stopifnot(file.exists(matrix_path), file.exists(metadata_path))


# ------------------------------------------------------------
# 3. Read the relevant Shahan atlas metadata
# ------------------------------------------------------------

# Only the first two columns are required. The remaining Data S3 columns are
# deliberately skipped to keep this recovery step lightweight.
atlas_meta <- readxl::read_excel(
  metadata_path,
  sheet = "atlas meta data",
  skip = 1,
  col_types = c("text", "text", rep("skip", 58)),
  .name_repair = "minimal"
)

required_columns <- c("cell_barcode_id", "orig.ident")
if (!all(required_columns %in% names(atlas_meta))) {
  stop(
    "Data S3 is missing required columns: ",
    paste(setdiff(required_columns, names(atlas_meta)), collapse = ", ")
  )
}

denyer_meta <- atlas_meta[
  atlas_meta$orig.ident %in% c("dc1", "dc2"),
  required_columns,
  drop = FALSE
]

if (nrow(denyer_meta) == 0L) {
  stop("Data S3 contains no cells labelled dc1 or dc2.")
}

metadata_suffix <- sub("^.*_", "", denyer_meta$cell_barcode_id)
expected_suffix <- ifelse(denyer_meta$orig.ident == "dc1", "15", "16")
if (!all(metadata_suffix == expected_suffix)) {
  stop("Unexpected Data S3 barcode suffix for dc1/dc2; refusing to infer replicate identity.")
}


# ------------------------------------------------------------
# 4. Read and validate the Denyer matrix barcodes
# ------------------------------------------------------------

# Only the compressed matrix header is read here. This is sufficient to
# recover all cell-barcode columns without loading the expression matrix.

matrix_connection <- gzfile(matrix_path, open = "rt")
matrix_header <- tryCatch(
  readLines(matrix_connection, n = 1L),
  finally = close(matrix_connection)
)

matrix_fields <- strsplit(matrix_header, ",", fixed = TRUE)[[1]]
matrix_fields <- sub("^\"|\"$", "", matrix_fields)
matrix_barcodes <- matrix_fields[-1L]

if (length(matrix_barcodes) != 4727L) {
  stop(
    "Expected 4,727 cell columns in the GSE123818 WT matrix; found ",
    length(matrix_barcodes), "."
  )
}
if (anyDuplicated(matrix_barcodes)) {
  stop("The GSE123818 WT matrix contains duplicated full cell barcodes.")
}
if (!all(grepl("^[ACGTN]+-[12]$", matrix_barcodes))) {
  stop("Unexpected matrix barcode format; expected a nucleotide barcode ending in -1 or -2.")
}

normalized_matrix_barcode <- sub("-[12]$", "", matrix_barcodes)
normalized_metadata_barcode <- sub("_[0-9]+$", "", denyer_meta$cell_barcode_id)


# ------------------------------------------------------------
# 5. Match barcodes and validate the suffix crosswalk
# ------------------------------------------------------------

# Build a sample-aware key. Base-only matching is unsafe because independent
# 10x libraries can reuse the same nucleotide barcode.
metadata_matrix_key <- paste0(
  normalized_metadata_barcode,
  ifelse(denyer_meta$orig.ident == "dc1", "-1", "-2")
)

if (anyDuplicated(metadata_matrix_key)) {
  stop("Data S3 contains duplicated sample-aware dc1/dc2 barcode keys.")
}

metadata_index <- match(matrix_barcodes, metadata_matrix_key)
direct_match <- !is.na(metadata_index)

# Validate the suffix-to-sample crosswalk empirically before assigning cells
# absent from the final WT atlas metadata (typically cells removed by atlas QC).
direct_sample <- denyer_meta$orig.ident[metadata_index]
discordant_direct_matches <- sum(
  direct_match &
    ((grepl("-1$", matrix_barcodes) & direct_sample != "dc1") |
       (grepl("-2$", matrix_barcodes) & direct_sample != "dc2")),
  na.rm = TRUE
)

if (discordant_direct_matches != 0L) {
  stop(
    "Matrix suffixes conflict with Data S3 sample labels in ",
    discordant_direct_matches,
    " directly matched cells; refusing to infer the remaining cells."
  )
}


# ------------------------------------------------------------
# 6. Assign replicate labels and record evidence strength
# ------------------------------------------------------------

# Direct Data S3 overlaps and suffix-inferred cells are deliberately kept as
# separate match statuses. The final 4,727 assignments must not all be
# described as direct cross-validation.

source_sample <- ifelse(grepl("-1$", matrix_barcodes), "dc1", "dc2")
replicate <- ifelse(source_sample == "dc1", "rep1", "rep2")
source_gsm <- ifelse(
  source_sample == "dc1",
  "GSM3511858",
  "GSM3511859"
)

match_status <- ifelse(
  direct_match,
  "direct_DataS3_match",
  "suffix_inferred_not_in_final_WT_atlas_metadata"
)

evidence <- ifelse(
  direct_match,
  paste0(
    "Exact sample-aware match after ",
    ifelse(source_sample == "dc1", "_15 to -1", "_16 to -2"),
    " normalization"
  ),
  paste0(
    "Matrix suffix ",
    ifelse(source_sample == "dc1", "-1", "-2"),
    " assigned after the suffix-to-sample crosswalk was validated by all ",
    sum(direct_match),
    " Data S3 overlaps with zero discordance; cell is absent from final WT atlas metadata"
  )
)

provenance <- paste(
  "Denyer et al. 2019 GSE123818 processed WT matrix;",
  "Shahan et al. 2022 Data S3 (SuppData3_complete.xltx from Elsevier",
  "1-s2.0-S1534580722000338-mmc3.zip);",
  "GSE152766 crosswalk dc1=GSM3511858=rep1, dc2=GSM3511859=rep2"
)

mapping <- data.frame(
  original_barcode = matrix_barcodes,
  normalized_barcode = normalized_matrix_barcode,
  replicate = replicate,
  source_sample = source_sample,
  source_gsm = source_gsm,
  shahan_cell_barcode_id = ifelse(
    direct_match,
    denyer_meta$cell_barcode_id[metadata_index],
    NA_character_
  ),
  match_status = match_status,
  evidence = evidence,
  source_provenance = provenance,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 7. Write the mapping and print a validation summary
# ------------------------------------------------------------

# The TSV is a reusable reference input for the robustness workflow. The
# console summary makes direct matches, inferred assignments, and conflicts
# visible immediately after the script finishes.

write.table(
  mapping,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

base_barcode_counts <- table(mapping$normalized_barcode)
summary_lines <- c(
  sprintf("Matrix cells: %d", nrow(mapping)),
  sprintf(
    "Direct Data S3 matches: %d (%.2f%%)",
    sum(direct_match),
    100 * mean(direct_match)
  ),
  sprintf("Direct matches by sample: dc1=%d, dc2=%d",
          sum(direct_match & source_sample == "dc1"),
          sum(direct_match & source_sample == "dc2")),
  sprintf("Absent from final WT atlas metadata: %d",
          sum(!direct_match)),
  sprintf("Final assignments: rep1=%d, rep2=%d",
          sum(replicate == "rep1"),
          sum(replicate == "rep2")),
  sprintf("Duplicated full barcodes: %d", sum(duplicated(matrix_barcodes))),
  sprintf("Duplicated base barcodes across libraries: %d IDs (%d cells)",
          sum(base_barcode_counts > 1L),
          sum(base_barcode_counts[base_barcode_counts > 1L])),
  sprintf("Ambiguous sample-aware mappings: %d", 0L),
  sprintf("Discordant direct matches: %d", discordant_direct_matches),
  sprintf("Output: %s", output_path)
)

cat(paste(summary_lines, collapse = "\n"), "\n")
