#!/usr/bin/env Rscript

# ============================================================
# Robustness workflow for the Arabidopsis root single-cell analysis
#
# Main purpose:
#   Gather the completed robustness checks into one reproducible workflow
#   while leaving the main clustering and annotation analysis unchanged.
#
# Analyses:
#   1. Validate the recovered Denyer WT replicate identities.
#   2. Add replicate metadata to the unchanged QC/annotation object and assess
#      replicate mixing at UMAP, computational-cluster, and broad-population
#      levels.
#   3. Test whether best Denyer C0-C14 signature matches are stable when both
#      sides use equally sized top-50, top-100, or top-200 marker sets.
#   4. Quantify sensitivity to Denyer protoplasting-induced genes, including
#      the seven previously selected positive-interaction GO terms.
#
# This script deliberately reuses scripts 01 and 02. It does not change their
# QC, clustering, marker-calling, annotation, or strong-candidate parameters.
# Run from the project root:
#   Rscript scripts/03_robustness_checks.R
# ============================================================


# ------------------------------------------------------------
# 1. Check packages, inputs, and output directories
# ------------------------------------------------------------

# Fail early when a required package or project input is missing. This avoids
# partially written robustness outputs from an incomplete environment.

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "ggplot2",
  "readxl", "AnnotationDbi", "org.At.tair.db"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

required_project_files <- c(
  "scripts/01_qc_clustering_annotation.R",
  "scripts/02_interaction_gene_cell_context.R",
  "data/GSE123818_Root_single_cell_wt_datamatrix.csv.gz",
  "data/FULL_interaction_GO_enrichment_positive.tsv",
  "data/reference/Denyer2019_WT_barcode_to_replicate.tsv",
  "data/reference/Denyer2019_TableS1_protoplasting_induced_FCgt2_q_lt_0.05.tsv",
  "docs/reference/Denyer2019_TableS2_cluster_DEGs_markers_and_identity.zip"
)
missing_project_files <- required_project_files[!file.exists(required_project_files)]
if (length(missing_project_files) > 0L) {
  stop(
    "Run this script from the project root. Missing file(s): ",
    paste(missing_project_files, collapse = ", ")
  )
}

set.seed(42)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/reproducibility", recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

save_png_pdf <- function(plot, basename, width, height) {
  ggplot2::ggsave(
    filename = file.path("results/figures", paste0(basename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path("results/figures", paste0(basename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    bg = "white"
  )
}

# ------------------------------------------------------------
# 2. Validate the recovered replicate-identity mapping
# ------------------------------------------------------------

# Confirm the expected 4,727-cell mapping, barcode uniqueness, sample/GSM
# crosswalk, and evidence status. Direct Data S3 matches remain distinct from
# cells assigned by the empirically validated suffix convention.

replicate_map <- utils::read.delim(
  "data/reference/Denyer2019_WT_barcode_to_replicate.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

required_mapping_columns <- c(
  "original_barcode", "normalized_barcode", "replicate", "source_sample",
  "source_gsm", "shahan_cell_barcode_id", "match_status", "evidence",
  "source_provenance"
)
if (!all(required_mapping_columns %in% names(replicate_map))) {
  stop(
    "Replicate mapping is missing column(s): ",
    paste(setdiff(required_mapping_columns, names(replicate_map)), collapse = ", ")
  )
}

expected_status_levels <- c(
  "direct_DataS3_match",
  "suffix_inferred_not_in_final_WT_atlas_metadata"
)
expected_replicate <- ifelse(replicate_map$source_sample == "dc1", "rep1", "rep2")
expected_gsm <- ifelse(
  replicate_map$source_sample == "dc1", "GSM3511858", "GSM3511859"
)
expected_sample_from_suffix <- ifelse(
  grepl("-1$", replicate_map$original_barcode), "dc1", "dc2"
)
sample_aware_key <- paste(
  replicate_map$normalized_barcode,
  replicate_map$source_sample,
  sep = "::"
)
direct_rows <- replicate_map$match_status == "direct_DataS3_match"

stopifnot(
  nrow(replicate_map) == 4727L,
  !anyNA(replicate_map$original_barcode),
  !anyDuplicated(replicate_map$original_barcode),
  !anyDuplicated(sample_aware_key),
  all(grepl("^[ACGTN]+-[12]$", replicate_map$original_barcode)),
  all(replicate_map$normalized_barcode ==
        sub("-[12]$", "", replicate_map$original_barcode)),
  all(replicate_map$source_sample %in% c("dc1", "dc2")),
  all(replicate_map$replicate == expected_replicate),
  all(replicate_map$source_gsm == expected_gsm),
  all(replicate_map$source_sample == expected_sample_from_suffix),
  sum(replicate_map$replicate == "rep1") == 2367L,
  sum(replicate_map$replicate == "rep2") == 2360L,
  setequal(unique(replicate_map$match_status), expected_status_levels),
  sum(replicate_map$match_status == "direct_DataS3_match") == 4458L,
  sum(replicate_map$match_status ==
        "suffix_inferred_not_in_final_WT_atlas_metadata") == 269L,
  all(!is.na(replicate_map$shahan_cell_barcode_id[
    replicate_map$match_status == "direct_DataS3_match"
  ])),
  all(
    replicate_map$normalized_barcode[direct_rows] ==
      sub("_[0-9]+$", "", replicate_map$shahan_cell_barcode_id[direct_rows])
  ),
  all(is.na(replicate_map$shahan_cell_barcode_id[
    replicate_map$match_status ==
      "suffix_inferred_not_in_final_WT_atlas_metadata"
  ]))
)

direct_sample_from_shahan_suffix <- ifelse(
  grepl("_15$", replicate_map$shahan_cell_barcode_id), "dc1",
  ifelse(grepl("_16$", replicate_map$shahan_cell_barcode_id), "dc2", NA_character_)
)
discordant_direct_matches <- sum(
  direct_rows &
    (is.na(direct_sample_from_shahan_suffix) |
       direct_sample_from_shahan_suffix != replicate_map$source_sample)
)
ambiguous_sample_aware_mappings <- sum(duplicated(sample_aware_key))
stopifnot(
  discordant_direct_matches == 0L,
  ambiguous_sample_aware_mappings == 0L
)

# ------------------------------------------------------------
# 3. Reconstruct the unchanged main-analysis objects when needed
# ------------------------------------------------------------

# Reuse scripts 01 and 02 rather than duplicating their analysis code. In a
# clean session they are sourced automatically; their parameters and fixed
# annotation decisions are not redefined here.

required_from_01 <- c(
  "seurat_obj", "seurat_qc", "markers", "denyer_xlsm", "cluster_overview",
  "cluster_confidence_wide", "analysis_theme"
)
sourced_script_01 <- !all(vapply(required_from_01, exists, logical(1), inherits = TRUE))
if (sourced_script_01) {
  source("scripts/01_qc_clustering_annotation.R")
}
if (!all(vapply(required_from_01, exists, logical(1), inherits = TRUE))) {
  stop("Script 01 did not create every object required by the robustness workflow.")
}

required_from_02 <- c(
  "interaction_genes", "preference_strength_checked", "strong_candidates"
)
sourced_script_02 <- !all(vapply(required_from_02, exists, logical(1), inherits = TRUE))
if (sourced_script_02) {
  source("scripts/02_interaction_gene_cell_context.R")
}
if (!all(vapply(required_from_02, exists, logical(1), inherits = TRUE))) {
  stop("Script 02 did not create every object required by the robustness workflow.")
}

stopifnot(
  inherits(seurat_obj, "Seurat"),
  inherits(seurat_qc, "Seurat"),
  ncol(seurat_obj) == 4727L,
  ncol(seurat_qc) == 4685L,
  nrow(markers) > 0L,
  nrow(preference_strength_checked) == 135L,
  nrow(strong_candidates) == 56L,
  !anyDuplicated(preference_strength_checked$TAIR),
  setequal(colnames(seurat_obj), replicate_map$original_barcode),
  "cell_identity" %in% colnames(seurat_qc@meta.data),
  !anyNA(seurat_qc$cell_identity),
  length(unique(seurat_qc$cell_identity)) == 11L
)

# ------------------------------------------------------------
# 4. Add replicate metadata to the QC-passed Seurat object
# ------------------------------------------------------------

# Join by the complete original barcode, including its library suffix. Using
# the base 10x sequence alone would be unsafe because barcodes can be reused
# across independent libraries.

replicate_index <- match(colnames(seurat_qc), replicate_map$original_barcode)
stopifnot(
  length(replicate_index) == 4685L,
  !anyNA(replicate_index),
  !anyDuplicated(replicate_index),
  identical(
    replicate_map$original_barcode[replicate_index],
    colnames(seurat_qc)
  )
)

seurat_qc$original_barcode <- replicate_map$original_barcode[replicate_index]
seurat_qc$replicate <- factor(
  replicate_map$replicate[replicate_index],
  levels = c("rep1", "rep2")
)
seurat_qc$source_sample <- factor(
  replicate_map$source_sample[replicate_index],
  levels = c("dc1", "dc2")
)
seurat_qc$source_gsm <- replicate_map$source_gsm[replicate_index]
seurat_qc$replicate_match_status <- factor(
  replicate_map$match_status[replicate_index],
  levels = expected_status_levels
)

stopifnot(
  !anyNA(seurat_qc$replicate),
  !anyNA(seurat_qc$source_sample),
  !anyNA(seurat_qc$source_gsm),
  !anyNA(seurat_qc$replicate_match_status),
  sum(table(seurat_qc$replicate)) == 4685L,
  all(table(seurat_qc$replicate) > 0L)
)

# ------------------------------------------------------------
# 5. Assess replicate mixing across the reconstructed atlas
# ------------------------------------------------------------

# Visualize replicate labels on the existing UMAP and summarize rep1/rep2
# composition within computational clusters and broad populations. This is a
# visual/compositional sanity check, not proof that all batch effects are absent.

replicate_umap <- Seurat::DimPlot(
  seurat_qc,
  reduction = "umap",
  group.by = "replicate",
  shuffle = TRUE,
  seed = 42
) +
  analysis_theme +
  ggplot2::scale_color_manual(values = c("rep1" = "#0072B2", "rep2" = "#D55E00")) +
  ggplot2::labs(
    title = "Replicate composition across the reconstructed atlas",
    subtitle = "Visual sanity check; broad mixing does not rule out all batch effects",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Replicate"
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

save_png_pdf(replicate_umap, "12_replicate_umap", width = 8, height = 6)

cell_replicate_metadata <- dplyr::tibble(
  original_barcode = colnames(seurat_qc),
  cluster = as.character(SeuratObject::Idents(seurat_qc)),
  broad_population = as.character(seurat_qc$cell_identity),
  replicate = as.character(seurat_qc$replicate)
)
stopifnot(
  nrow(cell_replicate_metadata) == 4685L,
  !anyNA(cell_replicate_metadata),
  !anyDuplicated(cell_replicate_metadata$original_barcode),
  setequal(unique(cell_replicate_metadata$cluster), as.character(0:13)),
  setequal(unique(cell_replicate_metadata$replicate), c("rep1", "rep2"))
)

replicate_totals <- cell_replicate_metadata |>
  dplyr::count(replicate, name = "replicate_total")

cluster_replicate_composition <- cell_replicate_metadata |>
  dplyr::count(cluster, replicate, name = "n_cells") |>
  tidyr::complete(
    cluster = as.character(0:13),
    replicate = c("rep1", "rep2"),
    fill = list(n_cells = 0L)
  ) |>
  dplyr::group_by(cluster) |>
  dplyr::mutate(
    cluster_total = sum(n_cells),
    fraction_within_cluster = n_cells / cluster_total
  ) |>
  dplyr::ungroup() |>
  dplyr::left_join(replicate_totals, by = "replicate", relationship = "many-to-one") |>
  dplyr::mutate(fraction_of_replicate = n_cells / replicate_total) |>
  dplyr::arrange(as.integer(cluster), replicate)

population_levels <- c(
  "Mature", "Mature-like", "Meristem", "Meristem-like", "Atrichoblast",
  "Trichoblast", "Stele", "Xylem", "Endodermis", "Cortex", "QC/Columella"
)
stopifnot(setequal(unique(cell_replicate_metadata$broad_population), population_levels))

population_replicate_composition <- cell_replicate_metadata |>
  dplyr::count(broad_population, replicate, name = "n_cells") |>
  tidyr::complete(
    broad_population = population_levels,
    replicate = c("rep1", "rep2"),
    fill = list(n_cells = 0L)
  ) |>
  dplyr::group_by(broad_population) |>
  dplyr::mutate(
    population_total = sum(n_cells),
    fraction_within_population = n_cells / population_total
  ) |>
  dplyr::ungroup() |>
  dplyr::left_join(replicate_totals, by = "replicate", relationship = "many-to-one") |>
  dplyr::mutate(fraction_of_replicate = n_cells / replicate_total) |>
  dplyr::mutate(
    broad_population = factor(broad_population, levels = population_levels)
  ) |>
  dplyr::arrange(broad_population, replicate) |>
  dplyr::mutate(broad_population = as.character(broad_population))

stopifnot(
  nrow(cluster_replicate_composition) == 14L * 2L,
  nrow(population_replicate_composition) == 11L * 2L,
  sum(cluster_replicate_composition$n_cells) == 4685L,
  sum(population_replicate_composition$n_cells) == 4685L,
  all(cluster_replicate_composition$n_cells > 0L),
  all(population_replicate_composition$n_cells > 0L),
  all(abs(
    cluster_replicate_composition |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(x = sum(fraction_within_cluster), .groups = "drop") |>
      dplyr::pull(x) - 1
  ) < 1e-12),
  all(abs(
    population_replicate_composition |>
      dplyr::group_by(broad_population) |>
      dplyr::summarise(x = sum(fraction_within_population), .groups = "drop") |>
      dplyr::pull(x) - 1
  ) < 1e-12)
)

write_tsv(
  cluster_replicate_composition,
  "results/tables/robustness_cluster_replicate_composition.tsv"
)
write_tsv(
  population_replicate_composition,
  "results/tables/robustness_population_replicate_composition.tsv"
)

# ------------------------------------------------------------
# 6. Test sensitivity to marker-set size
# ------------------------------------------------------------

# Compare equally sized reconstructed and Denyer marker sets at top 50, 100,
# and 200. This tests whether the existing best reference match depends on the
# original top-100 choice; it does not redefine the biological annotation.

cutoffs <- c(50L, 100L, 200L)
reconstructed_clusters <- as.character(0:13)
denyer_clusters <- as.character(0:14)

marker_pool <- markers |>
  dplyr::transmute(
    cluster = as.character(cluster),
    gene = as.character(gene),
    avg_log2FC = as.numeric(avg_log2FC)
  ) |>
  dplyr::filter(!is.na(gene), !is.na(avg_log2FC)) |>
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC), gene) |>
  dplyr::distinct(cluster, gene, .keep_all = TRUE)

marker_pool_counts <- marker_pool |>
  dplyr::count(cluster, name = "n_available")
stopifnot(
  setequal(marker_pool_counts$cluster, reconstructed_clusters),
  all(marker_pool_counts$n_available >= max(cutoffs)),
  !anyDuplicated(marker_pool[c("cluster", "gene")])
)

denyer_identity <- cluster_overview |>
  dplyr::filter(as.character(Cluster) %in% denyer_clusters) |>
  dplyr::transmute(
    denyer_cluster = as.character(Cluster),
    denyer_identity = as.character(Identity)
  )
stopifnot(
  nrow(denyer_identity) == 15L,
  !anyDuplicated(denyer_identity$denyer_cluster),
  !anyNA(denyer_identity$denyer_identity)
)

denyer_signature_pool <- lapply(0:14, function(i) {
  signature <- readxl::read_excel(
    denyer_xlsm,
    sheet = paste0("C", i),
    .name_repair = "minimal"
  )
  if (ncol(signature) != 5L) {
    stop("Unexpected number of columns in Denyer sheet C", i, ".")
  }
  names(signature) <- c("gene", "avg_logFC", "p_value", "padj", "cluster")
  signature |>
    dplyr::transmute(
      denyer_cluster = as.character(i),
      gene = as.character(gene),
      avg_logFC = as.numeric(avg_logFC),
      padj = as.numeric(padj)
    ) |>
    dplyr::filter(
      !is.na(gene),
      !is.na(avg_logFC),
      !is.na(padj),
      avg_logFC > 0,
      padj < 0.05
    ) |>
    dplyr::arrange(dplyr::desc(avg_logFC), gene) |>
    dplyr::distinct(gene, .keep_all = TRUE)
}) |>
  dplyr::bind_rows()

denyer_pool_counts <- denyer_signature_pool |>
  dplyr::count(denyer_cluster, name = "n_available")
stopifnot(
  setequal(denyer_pool_counts$denyer_cluster, denyer_clusters),
  all(denyer_pool_counts$n_available >= max(cutoffs)),
  !anyDuplicated(denyer_signature_pool[c("denyer_cluster", "gene")])
)

score_one_cutoff <- function(cutoff) {
  our_sets <- marker_pool |>
    dplyr::group_by(cluster) |>
    dplyr::slice_head(n = cutoff) |>
    dplyr::summarise(genes = list(gene), .groups = "drop")
  denyer_sets <- denyer_signature_pool |>
    dplyr::group_by(denyer_cluster) |>
    dplyr::slice_head(n = cutoff) |>
    dplyr::summarise(genes = list(gene), .groups = "drop")

  stopifnot(
    all(lengths(our_sets$genes) == cutoff),
    all(lengths(denyer_sets$genes) == cutoff),
    all(vapply(our_sets$genes, function(x) !anyDuplicated(x), logical(1))),
    all(vapply(denyer_sets$genes, function(x) !anyDuplicated(x), logical(1)))
  )

  tidyr::crossing(
    cluster = reconstructed_clusters,
    denyer_cluster = denyer_clusters
  ) |>
    dplyr::left_join(our_sets, by = "cluster", relationship = "many-to-one") |>
    dplyr::rename(our_genes = genes) |>
    dplyr::left_join(
      denyer_sets,
      by = "denyer_cluster",
      relationship = "many-to-one"
    ) |>
    dplyr::rename(denyer_genes = genes) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      n_overlap = length(intersect(our_genes, denyer_genes)),
      jaccard = n_overlap / (2 * cutoff - n_overlap)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-our_genes, -denyer_genes) |>
    dplyr::left_join(
      denyer_identity,
      by = "denyer_cluster",
      relationship = "many-to-one"
    ) |>
    dplyr::group_by(cluster) |>
    dplyr::arrange(
      dplyr::desc(n_overlap),
      as.integer(denyer_cluster),
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(cutoff = cutoff) |>
    dplyr::select(
      cluster, cutoff, denyer_cluster, denyer_identity, n_overlap, jaccard
    )
}

best_matches_long <- lapply(cutoffs, score_one_cutoff) |>
  dplyr::bind_rows()
stopifnot(
  nrow(best_matches_long) == 14L * 3L,
  !anyDuplicated(best_matches_long[c("cluster", "cutoff")]),
  all(best_matches_long$n_overlap >= 0L),
  all(best_matches_long$jaccard >= 0 & best_matches_long$jaccard <= 1)
)

best_matches_wide <- best_matches_long |>
  dplyr::mutate(cutoff = paste0("top", cutoff)) |>
  tidyr::pivot_wider(
    names_from = cutoff,
    values_from = c(denyer_cluster, denyer_identity, n_overlap, jaccard),
    names_glue = "{.value}_{cutoff}"
  ) |>
  dplyr::mutate(
    stable_best_match =
      denyer_cluster_top50 == denyer_cluster_top100 &
      denyer_cluster_top100 == denyer_cluster_top200
  ) |>
  dplyr::arrange(as.integer(cluster))

top100_from_main <- cluster_confidence_wide |>
  dplyr::transmute(
    cluster = as.character(cluster),
    denyer_cluster_top100_main = as.character(denyer_cluster_1),
    n_overlap_top100_main = as.integer(n_overlap_1),
    jaccard_top100_main = as.numeric(jaccard_1)
  )
top100_crosscheck <- best_matches_wide |>
  dplyr::select(
    cluster, denyer_cluster_top100, n_overlap_top100, jaccard_top100
  ) |>
  dplyr::left_join(top100_from_main, by = "cluster", relationship = "one-to-one")
stopifnot(
  nrow(best_matches_wide) == 14L,
  sum(best_matches_wide$stable_best_match) == 14L,
  all(top100_crosscheck$denyer_cluster_top100 ==
        top100_crosscheck$denyer_cluster_top100_main),
  all(top100_crosscheck$n_overlap_top100 ==
        top100_crosscheck$n_overlap_top100_main),
  all(abs(top100_crosscheck$jaccard_top100 -
            top100_crosscheck$jaccard_top100_main) < 1e-12)
)

write_tsv(
  best_matches_wide,
  "results/tables/robustness_marker_cutoff_sensitivity.tsv"
)

# ------------------------------------------------------------
# 7. Quantify overlap with protoplasting-induced genes
# ------------------------------------------------------------

# Flag detected interaction genes using Denyer Table S1, then summarize the
# overlap for all detected genes, exploratory strong candidates, and positive
# versus negative interaction directions.

protoplast_ref <- utils::read.delim(
  "data/reference/Denyer2019_TableS1_protoplasting_induced_FCgt2_q_lt_0.05.tsv",
  check.names = TRUE,
  stringsAsFactors = FALSE
)
if (!"Gene.ID" %in% names(protoplast_ref)) {
  stop("Protoplasting reference does not contain the expected Gene.ID column.")
}
protoplast_gene_ids <- unique(stats::na.omit(as.character(protoplast_ref$Gene.ID)))
stopifnot(
  nrow(protoplast_ref) == 3545L,
  !anyNA(protoplast_ref$Gene.ID),
  all(nzchar(protoplast_ref$Gene.ID)),
  !anyDuplicated(protoplast_ref$Gene.ID),
  length(protoplast_gene_ids) == 3545L,
  all(nzchar(protoplast_gene_ids))
)

interaction_protoplast_check <- preference_strength_checked |>
  dplyr::mutate(
    TAIR = as.character(TAIR),
    protoplasting_induced = TAIR %in% protoplast_gene_ids
  )
stopifnot(
  nrow(interaction_protoplast_check) == 135L,
  !anyDuplicated(interaction_protoplast_check$TAIR),
  sum(interaction_protoplast_check$strong_candidate) == 56L,
  sum(interaction_protoplast_check$protoplasting_induced) == 40L
)

summarise_protoplasting <- function(x, analysis_set, interaction_direction) {
  dplyr::tibble(
    analysis_set = analysis_set,
    interaction_direction = interaction_direction,
    n_genes = nrow(x),
    n_protoplasting_induced = sum(x$protoplasting_induced),
    fraction_protoplasting_induced = mean(x$protoplasting_induced)
  )
}

protoplast_interaction_summary <- dplyr::bind_rows(
  summarise_protoplasting(
    interaction_protoplast_check,
    "detected_interaction_genes",
    "All"
  ),
  summarise_protoplasting(
    dplyr::filter(interaction_protoplast_check, strong_candidate),
    "strong_candidates",
    "All"
  ),
  summarise_protoplasting(
    dplyr::filter(interaction_protoplast_check, interaction_direction == "Negative"),
    "detected_interaction_genes",
    "Negative"
  ),
  summarise_protoplasting(
    dplyr::filter(interaction_protoplast_check, interaction_direction == "Positive"),
    "detected_interaction_genes",
    "Positive"
  )
)

expected_protoplast_interaction_counts <- dplyr::tibble(
  analysis_set = c(
    "detected_interaction_genes", "strong_candidates",
    "detected_interaction_genes", "detected_interaction_genes"
  ),
  interaction_direction = c("All", "All", "Negative", "Positive"),
  n_genes = c(135L, 56L, 80L, 55L),
  n_protoplasting_induced = c(40L, 20L, 13L, 27L)
)
stopifnot(identical(
  protoplast_interaction_summary |>
    dplyr::select(-fraction_protoplasting_induced),
  expected_protoplast_interaction_counts
))

write_tsv(
  protoplast_interaction_summary,
  "results/tables/robustness_protoplasting_interaction_overlap.tsv"
)

# ------------------------------------------------------------
# 8. Repeat the GO cellular-context summary after exclusion
# ------------------------------------------------------------

# Reconstruct the previously examined GO-gene context from GOALL so that
# inherited Biological Process annotations are included. Restrict to positive
# interaction genes, matching the original positive-interaction GO analysis.
# The exclusion analysis asks whether the Mature/Stele context disappears when
# genes previously reported as protoplasting-induced are removed.
go_terms <- c(
  "cold acclimation",
  "response to cold",
  "response to acid chemical",
  "response to water",
  "response to water deprivation",
  "jasmonic acid metabolic process",
  "long-chain fatty acid metabolic process"
)
expected_go_sizes <- c(3L, 7L, 6L, 6L, 6L, 4L, 4L)
names(expected_go_sizes) <- go_terms
expected_go_flagged <- c(3L, 5L, 4L, 4L, 4L, 2L, 2L)
names(expected_go_flagged) <- go_terms

positive_context <- interaction_protoplast_check |>
  dplyr::filter(interaction_direction == "Positive")
go_annotations <- AnnotationDbi::select(
  org.At.tair.db::org.At.tair.db,
  keys = positive_context$TAIR,
  keytype = "TAIR",
  columns = c("GOALL", "ONTOLOGYALL")
)
go_annotations <- go_annotations |>
  dplyr::filter(ONTOLOGYALL == "BP", !is.na(GOALL)) |>
  dplyr::mutate(GO_term = AnnotationDbi::Term(GOALL)) |>
  dplyr::filter(GO_term %in% go_terms) |>
  dplyr::select(TAIR, GO_ID = GOALL, GO_term) |>
  dplyr::distinct(TAIR, GO_term, .keep_all = TRUE)

go_gene_context_sensitivity <- positive_context |>
  dplyr::inner_join(
    go_annotations,
    by = "TAIR",
    relationship = "one-to-many"
  ) |>
  dplyr::mutate(GO_term = factor(GO_term, levels = go_terms)) |>
  dplyr::arrange(GO_term, TAIR) |>
  dplyr::mutate(GO_term = as.character(GO_term))

go_protoplast_summary <- go_gene_context_sensitivity |>
  dplyr::group_by(GO_term) |>
  dplyr::summarise(
    n_genes = dplyr::n(),
    n_protoplasting_induced = sum(protoplasting_induced),
    fraction_protoplasting_induced = mean(protoplasting_induced),
    .groups = "drop"
  ) |>
  dplyr::mutate(GO_term = factor(GO_term, levels = go_terms)) |>
  dplyr::arrange(GO_term) |>
  dplyr::mutate(GO_term = as.character(GO_term))

stopifnot(
  nrow(go_gene_context_sensitivity) == sum(expected_go_sizes),
  dplyr::n_distinct(go_gene_context_sensitivity$TAIR) == 13L,
  identical(go_protoplast_summary$GO_term, go_terms),
  identical(go_protoplast_summary$n_genes, unname(expected_go_sizes)),
  identical(
    go_protoplast_summary$n_protoplasting_induced,
    unname(expected_go_flagged)
  )
)

go_context_no_protoplast <- go_gene_context_sensitivity |>
  dplyr::filter(!protoplasting_induced) |>
  dplyr::count(GO_term, preferred_population, name = "n_genes") |>
  dplyr::group_by(GO_term) |>
  dplyr::mutate(
    total_remaining = sum(n_genes),
    gene_fraction = n_genes / total_remaining
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(GO_term = factor(GO_term, levels = go_terms)) |>
  dplyr::arrange(GO_term, preferred_population) |>
  dplyr::mutate(GO_term = as.character(GO_term))

remaining_go_summary <- go_context_no_protoplast |>
  dplyr::distinct(GO_term, total_remaining)
stopifnot(
  nrow(go_context_no_protoplast) == 12L,
  !"cold acclimation" %in% go_context_no_protoplast$GO_term,
  nrow(remaining_go_summary) == 6L,
  all(remaining_go_summary$total_remaining == 2L),
  all(go_context_no_protoplast$n_genes == 1L),
  all(abs(go_context_no_protoplast$gene_fraction - 0.5) < 1e-12),
  setequal(
    dplyr::filter(
      go_context_no_protoplast,
      GO_term != "response to cold"
    )$preferred_population,
    c("Mature", "Stele")
  ),
  setequal(
    dplyr::filter(
      go_context_no_protoplast,
      GO_term == "response to cold"
    )$preferred_population,
    c("Cortex", "Mature")
  )
)

write_tsv(
  go_protoplast_summary,
  "results/tables/robustness_protoplasting_GO_overlap.tsv"
)
write_tsv(
  go_context_no_protoplast,
  "results/tables/robustness_GO_context_without_protoplasting_genes.tsv"
)

# ------------------------------------------------------------
# 9. Write reproducibility records and interpretation boundaries
# ------------------------------------------------------------

# Save the key expected counts together with explicit limits on interpretation.
# Also record package versions, session information, and the complete list of
# required outputs so a partial run cannot be mistaken for a successful one.

n_clusters_both_replicates <- cluster_replicate_composition |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(both_replicates = all(n_cells > 0L), .groups = "drop") |>
  dplyr::summarise(n = sum(both_replicates)) |>
  dplyr::pull(n)
n_populations_both_replicates <- population_replicate_composition |>
  dplyr::group_by(broad_population) |>
  dplyr::summarise(both_replicates = all(n_cells > 0L), .groups = "drop") |>
  dplyr::summarise(n = sum(both_replicates)) |>
  dplyr::pull(n)
stopifnot(
  n_clusters_both_replicates == 14L,
  n_populations_both_replicates == 11L
)

robustness_validation <- c(
  "Arabidopsis root single-cell robustness workflow validation",
  "",
  "Replicate identity recovery",
  sprintf("- Mapping rows: %d (expected 4727)", nrow(replicate_map)),
  sprintf(
    "- Direct Data S3 cross-validation: %d/4727 (%.2f%%)",
    sum(direct_rows),
    100 * mean(direct_rows)
  ),
  sprintf(
    "- Suffix-inferred cells absent from final WT atlas metadata: %d/4727 (%.2f%%)",
    sum(!direct_rows),
    100 * mean(!direct_rows)
  ),
  sprintf("- Discordant direct matches: %d", discordant_direct_matches),
  sprintf("- Ambiguous sample-aware mappings: %d", ambiguous_sample_aware_mappings),
  "- Evidence boundary: 4458 assignments are direct sample-aware cross-validation; 269 are suffix inference after the crosswalk was validated with zero discordance. The 4727 assignments must not all be described as direct matches.",
  "",
  "QC-passed object and replicate/batch sanity check",
  sprintf("- QC-passed cells with replicate labels: %d/%d", sum(!is.na(seurat_qc$replicate)), ncol(seurat_qc)),
  sprintf("- Replicate totals after QC: rep1=%d; rep2=%d", sum(seurat_qc$replicate == "rep1"), sum(seurat_qc$replicate == "rep2")),
  sprintf("- Computational clusters represented by both replicates: %d/14", n_clusters_both_replicates),
  sprintf("- Broad populations represented by both replicates: %d/11", n_populations_both_replicates),
  "- Interpretation: no obvious replicate-driven segregation; broad mixing is observed.",
  "- Boundary: this visual/compositional sanity check does not demonstrate complete absence of batch effects and is not formal replicate-level differential testing.",
  "",
  "Marker-cutoff sensitivity",
  "- Comparison: reconstructed top50 vs Denyer top50; top100 vs top100; top200 vs top200.",
  "- Jaccard formula at cutoff k: overlap / (2*k - overlap).",
  sprintf("- Stable best Denyer match across all three cutoffs: %d/14", sum(best_matches_wide$stable_best_match)),
  "- Boundary: this is a sensitivity check on the existing annotation evidence; it does not redefine clustering or annotation.",
  "",
  "Protoplasting sensitivity",
  sprintf("- Detected interaction genes flagged: %d/%d (%.1f%%)", 40L, 135L, 100 * 40 / 135),
  sprintf("- Strong candidates flagged: %d/%d (%.1f%%)", 20L, 56L, 100 * 20 / 56),
  sprintf("- Negative interaction genes flagged: %d/%d (%.1f%%)", 13L, 80L, 100 * 13 / 80),
  sprintf("- Positive interaction genes flagged: %d/%d (%.1f%%)", 27L, 55L, 100 * 27 / 55),
  "- cold acclimation: 3/3 genes flagged; this term has no genes remaining after exclusion.",
  "- Other GO terms: response to cold 5/7; response to acid chemical 4/6; response to water 4/6; response to water deprivation 4/6; jasmonic acid metabolic process 2/4; long-chain fatty acid metabolic process 2/4 flagged.",
  "- Interpretation: protoplasting is an important confounder, especially for positive interaction/stress-related genes. After exclusion, five terms retain a Mature/Stele pattern and response to cold retains Mature/Cortex, so the Mature/Stele cellular context is not completely explained by flagged genes.",
  "- Small-n boundary: each surviving GO term contains only 2 genes, and the GO terms share genes; this is supportive sensitivity evidence, not independent or strong statistical confirmation.",
  "",
  "Execution",
  sprintf("- Script 01 sourced in this session: %s", sourced_script_01),
  sprintf("- Script 02 sourced in this session: %s", sourced_script_02),
  "- Main-analysis QC, clustering, annotation, and candidate thresholds were reused unchanged."
)
writeLines(
  robustness_validation,
  "results/reproducibility/robustness_validation.txt"
)

writeLines(
  sub("[[:space:]]+$", "", capture.output(sessionInfo())),
  "results/reproducibility/robustness_sessionInfo.txt"
)
robustness_package_versions <- dplyr::tibble(
  package = required_packages,
  version = vapply(
    required_packages,
    function(package) as.character(utils::packageVersion(package)),
    character(1)
  )
)
write_tsv(
  robustness_package_versions,
  "results/reproducibility/robustness_package_versions.tsv"
)

expected_outputs <- c(
  "results/tables/robustness_cluster_replicate_composition.tsv",
  "results/tables/robustness_population_replicate_composition.tsv",
  "results/tables/robustness_marker_cutoff_sensitivity.tsv",
  "results/tables/robustness_protoplasting_interaction_overlap.tsv",
  "results/tables/robustness_protoplasting_GO_overlap.tsv",
  "results/tables/robustness_GO_context_without_protoplasting_genes.tsv",
  "results/figures/12_replicate_umap.png",
  "results/figures/12_replicate_umap.pdf",
  "results/reproducibility/robustness_validation.txt",
  "results/reproducibility/robustness_sessionInfo.txt",
  "results/reproducibility/robustness_package_versions.tsv"
)
stopifnot(
  all(file.exists(expected_outputs)),
  all(file.info(expected_outputs)$size > 0L)
)

message(
  "Robustness workflow completed: 4685/4685 QC-passed cells labelled; ",
  "14/14 marker matches stable; protoplasting overlap reproduced; ",
  "all expected outputs written."
)
