###############################################################################
# Project 1 interaction-gene contextualization in the Project 2 cell atlas
# Run directly with Rscript, or source after script 01. If the annotated
# seurat_qc object is absent, script 01 is sourced automatically.
###############################################################################

# 0. Reconstruct or validate the annotated single-cell object -----------------
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

# 1. Read and classify Project 1 interaction genes ----------------------------
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

# 2. Verify presence and detection in the single-cell matrix ------------------
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

# 3. Summarize expression across the 11 broad populations ---------------------
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

# 4. Build an explicit 135 x 11 matrix; cluster genes for display only --------
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
gene_clustering <- stats::hclust(
  stats::dist(expression_matrix, method = "euclidean"), method = "ward.D2"
)
clustered_gene_order <- rownames(expression_matrix)[gene_clustering$order]
heatmap_data <- population_summary %>%
  dplyr::mutate(
    TAIR = factor(TAIR, levels = rev(clustered_gene_order)),
    population = factor(population, levels = population_levels)
  )
clustered_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(x = population, y = TAIR, fill = expression_z)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
  ) +
  analysis_theme +
  ggplot2::labs(
    title = "Cell-population context of interaction genes",
    subtitle = "Per-gene z-scores; genes ordered by Euclidean / ward.D2 clustering",
    x = "Broad cell population", y = "135 detected interaction genes",
    fill = "Expression\nz-score"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )
print(clustered_heatmap)
save_figure(clustered_heatmap, "09_interaction_gene_clustered_heatmap",
            width = 10, height = 11)

# 5. Preferred population, top1-vs-top2 delta-z, and detection support --------
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

# 6. Final highlight scatter --------------------------------------------------
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
  ggplot2::scale_shape_manual(values = c("Negative" = 17, "Positive" = 16)) +
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
print(highlight_scatter)
save_figure(highlight_scatter, "10_interaction_gene_highlight_scatter",
            width = 10, height = 7)

# 7. Final reproducibility record ---------------------------------------------
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
  "Gene order only: Euclidean distance + ward.D2 hierarchical clustering"
)
writeLines(validation_lines,
           "results/reproducibility/interaction_gene_validation.txt")
print(strong_candidate_counts, n = Inf)
message(
  "Completed: 137 input; 136 present; 135 detected; 135 x 11 = 1485 summaries; ",
  "56 exploratory strong candidates."
)
