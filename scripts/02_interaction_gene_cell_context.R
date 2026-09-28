# ============================================================
# Arabidopsis root single-cell RNA-seq
# Cellular context of Project 1 genotype × drought interaction genes
#
# Main purpose:
#   Use the untreated WT root atlas from Project 2 to describe where the
#   Project 1 interaction genes are normally detected across broad root
#   cell populations.
#
# Important interpretation:
#   The atlas provides cellular context only. It cannot demonstrate that a
#   BRL3 × drought interaction occurs within a particular cell population.
#
# Run directly with Rscript, or source after script 01. If the annotated
# Seurat object is absent, script 01 is sourced automatically.
# ============================================================


# ------------------------------------------------------------
# 1. Reconstruct or validate the annotated single-cell object
# ------------------------------------------------------------

# Reuse the unchanged object produced by script 01. In a clean R session,
# source the complete QC, clustering, and annotation workflow first.

required_from_01 <- c("seurat_qc", "analysis_theme", "save_figure", "save_table")
if (!all(vapply(required_from_01, exists, logical(1), inherits = TRUE))) {
  source("scripts/01_qc_clustering_annotation.R")
}
stopifnot(inherits(seurat_qc, "Seurat"),
          "cell_identity" %in% colnames(seurat_qc@meta.data),
          !anyNA(seurat_qc$cell_identity),
          length(unique(seurat_qc$cell_identity)) == 11L)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Read and classify Project 1 interaction genes
# ------------------------------------------------------------

# Retain the Project 1 significance definition and classify direction from
# the interaction log2 fold change. No new differential-expression testing
# is performed in the single-cell dataset.

interaction_genes <- read.delim(
  "data/FULL_interaction_DEGs_padj0.05_log2FC1.tsv",
  check.names = FALSE, stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    TAIR = as.character(TAIR),
    interaction_direction = dplyr::if_else(
      log2FoldChange > 0, "Positive", "Negative"
    )
  )
stopifnot(nrow(interaction_genes) == 137L,
          !anyDuplicated(interaction_genes$TAIR),
          !anyNA(interaction_genes$interaction_direction))


# ------------------------------------------------------------
# 3. Verify presence and detection in the single-cell matrix
# ------------------------------------------------------------

# Separate genes that are absent from the atlas feature list from genes that
# are present but have zero counts in every retained cell. Only detected genes
# proceed to population-level expression summaries.

raw_counts <- SeuratObject::LayerData(seurat_qc, assay = "RNA", layer = "counts")
normalized_expression <- SeuratObject::LayerData(
  seurat_qc, assay = "RNA", layer = "data"
)
present_genes <- intersect(interaction_genes$TAIR, rownames(raw_counts))
detected_present <- Matrix::rowSums(
  raw_counts[present_genes, , drop = FALSE] > 0
) > 0
detected_genes <- names(detected_present)[detected_present]
gene_input_status <- interaction_genes %>%
  dplyr::transmute(
    TAIR, SYMBOL, interaction_direction,
    present_in_atlas = TAIR %in% present_genes,
    detected_in_atlas = TAIR %in% detected_genes
  )
stopifnot(
  length(present_genes) == 136L,
  length(detected_genes) == 135L,
  identical(gene_input_status$TAIR[!gene_input_status$present_in_atlas], "AT5G07985"),
  identical(
    gene_input_status$TAIR[
      gene_input_status$present_in_atlas & !gene_input_status$detected_in_atlas
    ], "AT1G66950"
  )
)
save_table(gene_input_status, "interaction_gene_input_status")


# ------------------------------------------------------------
# 4. Summarize expression across the 11 broad populations
# ------------------------------------------------------------

# For each detected gene, calculate the fraction of cells with non-zero
# counts and the mean normalized expression in every broad population.
# Per-gene z-scores describe relative expression across populations.

population_levels <- c(
  "Mature", "Mature-like", "Meristem", "Meristem-like", "Atrichoblast",
  "Trichoblast", "Stele", "Xylem", "Endodermis", "Cortex", "QC/Columella"
)
cell_population <- factor(
  as.character(seurat_qc$cell_identity), levels = population_levels
)
stopifnot(!anyNA(cell_population), all(table(cell_population) > 0L))

population_summary <- lapply(population_levels, function(population) {
  population_cells <- which(cell_population == population)
  data.frame(
    TAIR = detected_genes,
    population = population,
    n_cells = length(population_cells),
    detection_rate = as.numeric(Matrix::rowMeans(
      raw_counts[detected_genes, population_cells, drop = FALSE] > 0
    )),
    average_normalized_expression = as.numeric(Matrix::rowMeans(
      normalized_expression[detected_genes, population_cells, drop = FALSE]
    )),
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(
    interaction_genes %>%
      dplyr::select(TAIR, SYMBOL, log2FoldChange, interaction_direction),
    by = "TAIR", relationship = "many-to-one"
  ) %>%
  dplyr::group_by(TAIR) %>%
  dplyr::mutate(
    expression_sd_across_populations = stats::sd(average_normalized_expression),
    expression_z = dplyr::if_else(
      expression_sd_across_populations > 0,
      (average_normalized_expression - mean(average_normalized_expression)) /
        expression_sd_across_populations,
      0
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    TAIR, SYMBOL, interaction_direction, log2FoldChange, population, n_cells,
    detection_rate, average_normalized_expression, expression_z
  )
stopifnot(nrow(population_summary) == 135L * 11L,
          nrow(population_summary) == 1485L,
          !anyNA(population_summary$expression_z),
          all(dplyr::count(population_summary, TAIR)$n == 11L))
save_table(population_summary, "population_summary")


# ------------------------------------------------------------
# 5. Build the 135 × 11 display matrix
# ------------------------------------------------------------

# Construct an explicit numeric matrix and cluster genes only to determine
# heatmap row order. This display clustering does not define candidates or
# change the biological population annotations.

# Convert to a base data.frame before assigning row names. This avoids tibble
# row-name and list-column failure modes.
expression_wide <- population_summary %>%
  dplyr::select(TAIR, population, expression_z) %>%
  tidyr::pivot_wider(names_from = population, values_from = expression_z) %>%
  as.data.frame(check.names = FALSE)
expression_matrix <- as.matrix(expression_wide[, population_levels, drop = FALSE])
rownames(expression_matrix) <- expression_wide$TAIR
storage.mode(expression_matrix) <- "double"
stopifnot(identical(dim(expression_matrix), c(135L, 11L)),
          length(expression_matrix) == 1485L,
          !is.list(expression_matrix), !anyNA(expression_matrix),
          identical(colnames(expression_matrix), population_levels))

# Euclidean/ward.D2 controls only heatmap gene order, not candidate status.
# Cluster Positive and Negative interaction genes separately so the two
# directions remain visually distinct rather than forming one mixed block.
interaction_direction_levels <- c("Positive", "Negative")
gene_direction <- population_summary %>%
  dplyr::distinct(TAIR, interaction_direction)
stopifnot(
  !anyDuplicated(gene_direction$TAIR),
  identical(
    as.integer(table(factor(
      gene_direction$interaction_direction,
      levels = interaction_direction_levels
    ))),
    c(55L, 80L)
  )
)
gene_clustering_by_direction <- lapply(
  interaction_direction_levels,
  function(direction) {
    direction_genes <- gene_direction %>%
      dplyr::filter(interaction_direction == direction) %>%
      dplyr::pull(TAIR)
    stats::hclust(
      stats::dist(
        expression_matrix[direction_genes, , drop = FALSE],
        method = "euclidean"
      ),
      method = "ward.D2"
    )
  }
)
names(gene_clustering_by_direction) <- interaction_direction_levels
clustered_gene_order <- unlist(
  lapply(interaction_direction_levels, function(direction) {
    clustering <- gene_clustering_by_direction[[direction]]
    clustering$labels[clustering$order]
  }),
  use.names = FALSE
)
stopifnot(length(clustered_gene_order) == 135L,
          !anyDuplicated(clustered_gene_order),
          setequal(clustered_gene_order, rownames(expression_matrix)))
heatmap_data <- population_summary %>%
  dplyr::mutate(
    TAIR = factor(TAIR, levels = rev(clustered_gene_order)),
    population = factor(population, levels = population_levels),
    interaction_direction = factor(
      interaction_direction,
      levels = interaction_direction_levels
    )
  )
clustered_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(x = population, y = TAIR, fill = expression_z)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(interaction_direction),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  analysis_theme +
  ggplot2::labs(
    title = "Cell-population context of interaction genes",
    subtitle = paste0(
      "Per-gene relative expression z-scores: blue = higher, red = lower; ",
      "genes clustered within interaction direction"
    ),
    x = "Broad cell population", y = NULL,
    fill = "Relative\nexpression\nz-score"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold")
  )
if (interactive()) print(clustered_heatmap)
save_figure(clustered_heatmap, "09_interaction_gene_clustered_heatmap",
            width = 10, height = 11)


# ------------------------------------------------------------
# 6. Identify population-preferential exploratory candidates
# ------------------------------------------------------------

# Rank populations separately for each gene, then combine the top1-versus-top2
# z-score difference with detection support. The fixed exploratory threshold
# is delta-z >= 1 and detection in at least 10% of the preferred population.

ranked_population <- population_summary %>%
  dplyr::mutate(population = factor(population, levels = population_levels)) %>%
  dplyr::group_by(TAIR) %>%
  dplyr::arrange(dplyr::desc(expression_z), population, .by_group = TRUE) %>%
  dplyr::mutate(preference_rank = dplyr::row_number()) %>%
  dplyr::ungroup()
preferred_population <- ranked_population %>%
  dplyr::filter(preference_rank == 1L) %>%
  dplyr::transmute(
    TAIR, SYMBOL, interaction_direction, log2FoldChange,
    preferred_population = as.character(population),
    preferred_average_normalized_expression = average_normalized_expression,
    preferred_z = expression_z,
    preferred_population_detection_rate = detection_rate
  )
preference_strength_checked <- ranked_population %>%
  dplyr::filter(preference_rank <= 2L) %>%
  dplyr::select(
    TAIR, SYMBOL, interaction_direction, preference_rank, population,
    average_normalized_expression, expression_z, detection_rate
  ) %>%
  tidyr::pivot_wider(
    names_from = preference_rank,
    values_from = c(
      population, average_normalized_expression, expression_z, detection_rate
    ),
    names_glue = "{.value}_rank{preference_rank}"
  ) %>%
  dplyr::mutate(delta_z = expression_z_rank1 - expression_z_rank2) %>%
  dplyr::left_join(
    preferred_population %>%
      dplyr::select(TAIR, preferred_population, preferred_population_detection_rate),
    by = "TAIR", relationship = "one-to-one"
  ) %>%
  dplyr::mutate(
    strong_candidate =
      delta_z >= 1 & preferred_population_detection_rate >= 0.10
  )
stopifnot(nrow(preferred_population) == 135L,
          nrow(preference_strength_checked) == 135L,
          all(preference_strength_checked$population_rank1 ==
                preference_strength_checked$preferred_population))

strong_candidates <- preference_strength_checked %>%
  dplyr::filter(strong_candidate) %>%
  dplyr::arrange(
    interaction_direction, preferred_population, dplyr::desc(delta_z)
  )
strong_candidate_counts <- strong_candidates %>%
  dplyr::count(
    interaction_direction, preferred_population, name = "n_strong_candidates"
  ) %>%
  dplyr::arrange(interaction_direction, preferred_population)
expected_strong_candidate_counts <- dplyr::tribble(
  ~interaction_direction, ~preferred_population, ~n_strong_candidates,
  "Negative", "Atrichoblast", 2L,
  "Negative", "Cortex", 2L,
  "Negative", "Endodermis", 7L,
  "Negative", "Mature", 2L,
  "Negative", "QC/Columella", 1L,
  "Negative", "Stele", 9L,
  "Negative", "Trichoblast", 7L,
  "Positive", "Mature", 20L,
  "Positive", "Mature-like", 1L,
  "Positive", "QC/Columella", 1L,
  "Positive", "Stele", 4L
) %>%
  dplyr::arrange(interaction_direction, preferred_population)
stopifnot(nrow(strong_candidates) == 56L,
          identical(strong_candidate_counts, expected_strong_candidate_counts))
save_table(preferred_population, "preferred_population")
save_table(preference_strength_checked, "preference_strength_checked")
save_table(strong_candidates, "strong_candidates")
save_table(strong_candidate_counts, "strong_candidate_counts")


# ------------------------------------------------------------
# 7. Visualize preference strength and detection support
# ------------------------------------------------------------

# Show all detected interaction genes while highlighting only those that pass
# both exploratory criteria. Color represents the preferred broad population;
# point shape retains the Project 1 interaction direction.

population_colors <- c(
  "Mature" = "#D55E00", "Mature-like" = "#E69F00",
  "Meristem" = "#F0E442", "Meristem-like" = "#CC79A7",
  "Atrichoblast" = "#56B4E9", "Trichoblast" = "#0072B2",
  "Stele" = "#009E73", "Xylem" = "#1B9E77",
  "Endodermis" = "#7570B3", "Cortex" = "#E7298A",
  "QC/Columella" = "#A6761D"
)
scatter_data <- preference_strength_checked %>%
  dplyr::mutate(
    preferred_population = factor(preferred_population, levels = population_levels),
    interaction_direction = factor(
      interaction_direction, levels = c("Negative", "Positive")
    )
  )
highlight_scatter <- ggplot2::ggplot(
  scatter_data,
  ggplot2::aes(
    x = preferred_population_detection_rate, y = delta_z,
    shape = interaction_direction
  )
) +
  ggplot2::geom_vline(
    xintercept = 0.10, linetype = "dashed", linewidth = 0.5, color = "grey35"
  ) +
  ggplot2::geom_hline(
    yintercept = 1, linetype = "dashed", linewidth = 0.5, color = "grey35"
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(scatter_data, !strong_candidate),
    color = "grey25", alpha = 0.32, size = 2.4, stroke = 0.7
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(scatter_data, strong_candidate),
    ggplot2::aes(color = preferred_population),
    alpha = 0.9, size = 2.8, stroke = 0.8
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::label_percent(accuracy = 1), limits = c(0, NA),
    expand = ggplot2::expansion(mult = c(0.01, 0.05))
  ) +
  ggplot2::scale_color_manual(values = population_colors, drop = TRUE) +
  ggplot2::scale_shape_manual(values = c("Negative" = 16, "Positive" = 17)) +
  analysis_theme +
  ggplot2::labs(
    title = "Interaction genes with population-preferential expression",
    subtitle = paste0(
      "Colored points: delta-z >= 1 and detection >= 10% in the preferred population (n = ",
      nrow(strong_candidates), ")"
    ),
    x = "Detection rate in preferred population",
    y = "Preference strength (top1 - top2 z-score)",
    color = "Preferred population", shape = "Interaction direction"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.title = ggplot2::element_text(face = "bold")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      order = 1, override.aes = list(alpha = 1, size = 3)
    ),
    shape = ggplot2::guide_legend(
      order = 2, override.aes = list(alpha = 1, color = "grey20")
    )
  )
if (interactive()) print(highlight_scatter)
save_figure(highlight_scatter, "10_interaction_gene_highlight_scatter",
            width = 10, height = 7)


# ------------------------------------------------------------
# 8. Summarize enriched pathways across preferred populations
# ------------------------------------------------------------

# Reconstruct the seven previously selected positive-interaction GO terms from
# inherited Biological Process annotations. Each tile shows the fraction of a
# pathway's genes whose preferred atlas population is the indicated column.
# This is a descriptive cellular-context summary, not an enrichment test.

pathway_go_terms <- c(
  "cold acclimation",
  "response to cold",
  "response to acid chemical",
  "response to water",
  "response to water deprivation",
  "jasmonic acid metabolic process",
  "long-chain fatty acid metabolic process"
)
expected_pathway_gene_counts <- c(3L, 7L, 6L, 6L, 6L, 4L, 4L)
names(expected_pathway_gene_counts) <- pathway_go_terms

positive_preference_context <- preference_strength_checked %>%
  dplyr::filter(interaction_direction == "Positive")
positive_go_annotations <- AnnotationDbi::select(
  org.At.tair.db::org.At.tair.db,
  keys = positive_preference_context$TAIR,
  keytype = "TAIR",
  columns = c("GOALL", "ONTOLOGYALL")
) %>%
  dplyr::filter(ONTOLOGYALL == "BP", !is.na(GOALL)) %>%
  dplyr::mutate(GO_term = AnnotationDbi::Term(GOALL)) %>%
  dplyr::filter(GO_term %in% pathway_go_terms) %>%
  dplyr::transmute(TAIR, GO_ID = GOALL, GO_term) %>%
  dplyr::distinct(TAIR, GO_term, .keep_all = TRUE)

go_gene_context <- positive_preference_context %>%
  dplyr::inner_join(
    positive_go_annotations,
    by = "TAIR",
    relationship = "one-to-many"
  ) %>%
  dplyr::mutate(GO_term = factor(GO_term, levels = pathway_go_terms)) %>%
  dplyr::arrange(GO_term, TAIR) %>%
  dplyr::mutate(GO_term = as.character(GO_term))

observed_pathway_gene_counts <- go_gene_context %>%
  dplyr::count(GO_term, name = "n_genes") %>%
  dplyr::mutate(GO_term = factor(GO_term, levels = pathway_go_terms)) %>%
  dplyr::arrange(GO_term)
stopifnot(
  nrow(go_gene_context) == sum(expected_pathway_gene_counts),
  dplyr::n_distinct(go_gene_context$TAIR) == 13L,
  identical(as.character(observed_pathway_gene_counts$GO_term), pathway_go_terms),
  identical(
    observed_pathway_gene_counts$n_genes,
    unname(expected_pathway_gene_counts)
  )
)

pathway_population_levels <- c("Cortex", "Mature", "Mature-like", "Stele")
pathway_context <- go_gene_context %>%
  dplyr::count(GO_term, preferred_population, name = "n_genes") %>%
  dplyr::group_by(GO_term) %>%
  dplyr::mutate(
    pathway_total = sum(n_genes),
    gene_fraction = n_genes / pathway_total,
    fraction_label = scales::percent(gene_fraction, accuracy = 0.01)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    GO_term = factor(GO_term, levels = pathway_go_terms),
    preferred_population = factor(
      preferred_population,
      levels = pathway_population_levels
    )
  ) %>%
  dplyr::arrange(GO_term, preferred_population)
stopifnot(
  nrow(pathway_context) == 20L,
  sum(pathway_context$n_genes) == 36L,
  !anyNA(pathway_context$preferred_population),
  all(pathway_context$pathway_total ==
        unname(expected_pathway_gene_counts[as.character(pathway_context$GO_term)]))
)

pathway_heatmap <- ggplot2::ggplot(
  pathway_context,
  ggplot2::aes(
    x = preferred_population,
    y = GO_term,
    fill = gene_fraction
  )
) +
  ggplot2::geom_tile(color = "black", linewidth = 0.6) +
  ggplot2::geom_text(
    ggplot2::aes(label = fraction_label),
    size = 4
  ) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "#2166AC",
    labels = scales::label_percent(accuracy = 0.01)
  ) +
  analysis_theme +
  ggplot2::labs(
    title = "Cell-population context of enriched interaction pathways",
    x = "Preferred cell population",
    y = NULL,
    fill = "Gene fraction"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )
if (interactive()) print(pathway_heatmap)
save_figure(pathway_heatmap, "11_pathway_cell_context_heatmap",
            width = 12, height = 7)
save_table(go_gene_context, "pathway_GO_gene_context")
save_table(pathway_context, "pathway_cell_context")


# ------------------------------------------------------------
# 9. Write result tables and the reproducibility record
# ------------------------------------------------------------

# Record the expected input, detection, matrix, and candidate counts so that
# future reruns fail visibly if the upstream object or input files change.

# Untreated WT data provide cellular context only; they cannot prove that a
# BRL3 x drought effect occurs in the preferred population.
validation_lines <- c(
  "Project 2 interaction-gene cellular contextualization validation",
  paste("Interaction genes:", nrow(interaction_genes)),
  paste("Present in atlas:", length(present_genes)),
  paste("Detected in atlas:", length(detected_genes)),
  "Absent: AT5G07985",
  "Present but all-zero: AT1G66950",
  paste("Population-summary rows:", nrow(population_summary)),
  paste("Expression matrix:", paste(dim(expression_matrix), collapse = " x ")),
  paste("Strong candidates:", nrow(strong_candidates)),
  "Strong cutoff: delta_z >= 1 and preferred-population detection_rate >= 0.10",
  "Heatmap gene order only: Euclidean distance + ward.D2 clustering within interaction direction",
  paste("Pathway GO-term x gene rows:", nrow(go_gene_context)),
  paste("Pathway unique genes:", dplyr::n_distinct(go_gene_context$TAIR)),
  paste("Pathway population-summary rows:", nrow(pathway_context)),
  "Pathway context is descriptive and is not a pathway-enrichment statistic"
)
writeLines(validation_lines,
           "results/reproducibility/interaction_gene_validation.txt")
print(strong_candidate_counts, n = Inf)
message(
  "Completed: 137 input; 136 present; 135 detected; 135 x 11 = 1485 summaries; ",
  "56 exploratory strong candidates; 7 pathway contexts."
)
