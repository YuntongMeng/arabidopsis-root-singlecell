###############################################################################
# Arabidopsis root single-cell RNA-seq
# QC, dimensionality reduction, clustering, and cell-population annotation
#
# Dataset:
#   Denyer et al. (2019), Arabidopsis root single-cell RNA-seq
#
# Main purpose:
#   1. Perform an independent clustering of the published WT root scRNA-seq data.
#   2. Identify marker genes for the reconstructed clusters.
#   3. Compare reconstructed cluster signatures with Denyer et al. clusters.
#   4. Assign broad biological identities based on this comparison.
#
# Important interpretation:
#   Computational clusters and biological cell identities are NOT the same thing.
#   Clusters are generated from expression similarity.
#   Biological identities are assigned afterward using marker/signature evidence.
#
# Project question:
#   This atlas will later be used to contextualize BRL3 × drought interaction genes
#   identified in the independent bulk RNA-seq project.
###############################################################################


# =============================================================================
# 0. Load packages
# =============================================================================

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# Arabidopsis gene annotation
library(org.At.tair.db)
library(AnnotationDbi)


# Run from the project root, including when using RStudio Source.
set.seed(42)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/reproducibility", recursive = TRUE, showWarnings = FALSE)
analysis_theme <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.title = element_blank())
save_figure <- function(plot, name, width = 9, height = 6) {
  for (ext in c("png", "pdf")) {
    ggsave(file.path("results/figures", paste0(name, ".", ext)),
           plot = plot, width = width, height = height, dpi = 300, bg = "white")
  }
}
save_table <- function(x, name) {
  write.csv(x, file.path("results/tables", paste0(name, ".csv")), row.names = FALSE)
}

# =============================================================================
# 1. Load the published single-cell expression matrix
# =============================================================================

# Rows = genes
# Columns = individual cell barcodes
#
# The matrix contains 27,629 genes and 4,727 cells.

counts <- read.csv(
  "data/GSE123818_Root_single_cell_wt_datamatrix.csv.gz",
  row.names = 1,
  check.names = FALSE
)

dim(counts)

# Expected:
# 27629 genes × 4727 cells


# =============================================================================
# 2. Create Seurat object
# =============================================================================

seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "ArabidopsisRoot"
)

head(seurat_obj@meta.data)

# Important QC variables:
#
# nCount_RNA
#   Total RNA counts captured for each cell.
#
# nFeature_RNA
#   Number of genes detected with non-zero counts in each cell.


summary(seurat_obj$nCount_RNA)
summary(seurat_obj$nFeature_RNA)


# =============================================================================
# 3. Initial QC visualization
# =============================================================================

qc_violin <- VlnPlot(
  seurat_obj,
  features = c("nCount_RNA", "nFeature_RNA"),
  ncol = 2,
  pt.size = 0
) &
  analysis_theme
print(qc_violin)
save_figure(qc_violin, "01_qc_violin")

# Relationship between sequencing depth and detected gene number.
#
# Cells with unusually high values for BOTH metrics may represent
# potential doublets/multiplets, although these values alone cannot prove this.

qc_depth_features <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
) +
  analysis_theme +
  labs(title = "Sequencing depth and detected genes")
print(qc_depth_features)
save_figure(qc_depth_features, "02_qc_depth_features")


# =============================================================================
# 4. Mitochondrial and chloroplast transcript percentages
# =============================================================================

# TAIR locus prefixes:
#
# ATMG = mitochondrial genes
# ATCG = chloroplast genes

grep("^ATMG", rownames(seurat_obj), value = TRUE)[1:20]
grep("^ATCG", rownames(seurat_obj), value = TRUE)[1:20]


seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^ATMG"
)

seurat_obj[["percent.cp"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^ATCG"
)


summary(seurat_obj$percent.mt)
summary(seurat_obj$percent.cp)

# Note:
# PercentageFeatureSet() already returns percentage units.
# For example, 4.35 means 4.35%, not 0.0435.


qc_mitochondrial <- FeatureScatter(
  seurat_obj,
  feature1 = "nFeature_RNA",
  feature2 = "percent.mt"
) +
  analysis_theme +
  labs(title = "Mitochondrial transcript percentage")
print(qc_mitochondrial)
save_figure(qc_mitochondrial, "03_qc_mitochondrial")


qc_chloroplast <- FeatureScatter(
  seurat_obj,
  feature1 = "nFeature_RNA",
  feature2 = "percent.cp"
) +
  analysis_theme +
  labs(title = "Chloroplast transcript percentage")
print(qc_chloroplast)
save_figure(qc_chloroplast, "04_qc_chloroplast")


# In this dataset, mitochondrial/chloroplast percentages are generally low.
# Therefore, no hard mitochondrial or chloroplast cutoff was applied.


# =============================================================================
# 5. Examine high-count / high-feature tails
# =============================================================================

quantile(
  seurat_obj$nCount_RNA,
  probs = c(0.90, 0.95, 0.975, 0.99, 0.995, 0.999)
)

quantile(
  seurat_obj$nFeature_RNA,
  probs = c(0.90, 0.95, 0.975, 0.99, 0.995, 0.999)
)


# Conservative upper thresholds were chosen based on the extreme tails:
#
# nCount_RNA   < 220,000
# nFeature_RNA < 11,000
#
# No lower cutoff was applied because the published matrix was already processed
# and the low-count cells did not show an obvious low-quality population.


sum(seurat_obj$nCount_RNA >= 220000)
sum(seurat_obj$nFeature_RNA >= 11000)

sum(
  seurat_obj$nCount_RNA >= 220000 |
    seurat_obj$nFeature_RNA >= 11000
)


# =============================================================================
# 6. QC filtering
# =============================================================================

seurat_qc <- subset(
  seurat_obj,
  subset =
    nCount_RNA < 220000 &
    nFeature_RNA < 11000
)

dim(seurat_qc)

# Expected:
# 27,629 genes × 4,685 cells
#
# 42 extreme high-count/high-feature cells were removed.
#
# Important:
# These cells are described as potential multiplet/outlier candidates,
# NOT confirmed doublets.


summary(seurat_qc$nCount_RNA)
summary(seurat_qc$nFeature_RNA)


# =============================================================================
# 7. Normalize expression
# =============================================================================

seurat_qc <- NormalizeData(
  seurat_qc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Concept:
#
# First correct for differences in total RNA counts between cells,
# then log-transform the normalized expression values.


# =============================================================================
# 8. Identify highly variable genes
# =============================================================================

seurat_qc <- FindVariableFeatures(
  seurat_qc,
  selection.method = "vst",
  nfeatures = 2000
)

head(
  VariableFeatures(seurat_qc),
  20
)

# These 2,000 genes show relatively high biological variability across cells
# and are therefore informative for distinguishing different cell populations.


# =============================================================================
# 9. Scale variable genes
# =============================================================================

seurat_qc <- ScaleData(
  seurat_qc,
  features = VariableFeatures(seurat_qc)
)


# =============================================================================
# 10. Principal component analysis
# =============================================================================

seurat_qc <- RunPCA(
  seurat_qc,
  seed.use = 42,
  features = VariableFeatures(seurat_qc)
)


pca_elbow <- ElbowPlot(
  seurat_qc,
  ndims = 50
) +
  analysis_theme +
  labs(
    title = "PCA elbow plot"
  )
print(pca_elbow)
save_figure(pca_elbow, "05_pca_elbow")


# The elbow plot does not contain a single sharp breakpoint.
# The curve begins to flatten approximately around PCs 20–30.
#
# Therefore, 25 PCs were selected as a reasonable working choice.
#
# This value should be interpreted as a practical dimensionality choice,
# not as a uniquely correct mathematical answer.


# =============================================================================
# 11. Construct cell-cell neighbor graph
# =============================================================================

seurat_qc <- FindNeighbors(
  seurat_qc,
  dims = 1:25
)

# Each cell is represented using its coordinates in PCs 1–25.
# FindNeighbors identifies transcriptionally similar cells and builds
# a cell-cell graph.


# =============================================================================
# 12. Unsupervised clustering
# =============================================================================

seurat_qc <- FindClusters(
  seurat_qc,
  random.seed = 0,
  resolution = 0.5
)

table(Idents(seurat_qc))

length(levels(Idents(seurat_qc)))

# Result:
# 14 computational clusters (0–13)
#
# Important:
# A computational cluster is NOT automatically equivalent to a cell type.


# =============================================================================
# 13. UMAP visualization
# =============================================================================

seurat_qc <- RunUMAP(
  seurat_qc,
  seed.use = 42,
  dims = 1:25
)


umap_clusters <- DimPlot(
  seurat_qc,
  reduction = "umap",
  label = FALSE
) +
  analysis_theme +
  labs(
    title = "UMAP of Arabidopsis root single cells",
    subtitle = "Clusters identified at resolution = 0.5"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11)
  )


umap_clusters <- LabelClusters(
  plot = umap_clusters,
  id = "ident",
  repel = TRUE,
  size = 3,
  fontface = "bold",
  box = FALSE,
  seed = 42,
  min.segment.length = 0
)

print(umap_clusters)
save_figure(umap_clusters, "06_cluster_umap")

# Final figure formatting can be adjusted later.
# In particular, labels should remain readable without covering small clusters.


# =============================================================================
# 14. Identify marker genes for each reconstructed cluster
# =============================================================================

markers <- FindAllMarkers(
  seurat_qc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# FindAllMarkers tests which genes are enriched in each cluster relative
# to cells outside that cluster.
#
# only.pos = TRUE
#   Only positively enriched genes are retained.
#
# min.pct = 0.25
#   A gene must be detected in at least 25% of cells in the relevant comparison.
#
# logfc.threshold = 0.25
#   Removes very small expression differences.


# =============================================================================
# 15. Inspect top marker genes
# =============================================================================

top_markers <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  )


# Convert TAIR locus IDs into common Arabidopsis gene symbols where available.

top_markers$SYMBOL <- AnnotationDbi::mapIds(
  org.At.tair.db,
  keys = top_markers$gene,
  keytype = "TAIR",
  column = "SYMBOL",
  multiVals = "first"
)


top_markers %>%
  dplyr::select(
    cluster,
    gene,
    SYMBOL,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  ) %>%
  print(n = 140)


# Some TAIR loci have no commonly assigned SYMBOL.
# NA therefore does NOT indicate that the gene is missing from the genome.


# Human-readable marker summary.

marker_summary <- top_markers %>%
  dplyr::mutate(
    label = ifelse(
      is.na(SYMBOL),
      gene,
      paste0(SYMBOL, " (", gene, ")")
    )
  ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    top10_markers = paste(label, collapse = "; "),
    .groups = "drop"
  )

marker_summary


# =============================================================================
# 16. Load Denyer et al. supplementary cluster information
# =============================================================================

denyer_zip <-
  "docs/reference/Denyer2019_TableS2_cluster_DEGs_markers_and_identity.zip"

denyer_extract_dir <-
  "docs/reference/Denyer2019_TableS2_extracted"


dir.create(
  denyer_extract_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


unzip(
  denyer_zip,
  exdir = denyer_extract_dir
)


denyer_xlsm <- file.path(
  denyer_extract_dir,
  "devcel_4493_mmc3.xlsm"
)


excel_sheets(denyer_xlsm)


# Relevant sheets include:
#
# Cluster Overview
# C0 ... C14
# Marker Genes
#
# The original source also contains subcluster analyses and specificity scores.


# =============================================================================
# 17. Published Denyer cluster identities
# =============================================================================

cluster_overview <- read_excel(
  denyer_xlsm,
  sheet = "Cluster Overview"
)


cluster_overview %>%
  dplyr::select(
    Cluster,
    `Total cell number`,
    Identity
  ) %>%
  dplyr::filter(Cluster %in% as.character(0:14))


# Published broad identities include:
#
# Mature
# Meristem
# Atrichoblast
# Trichoblast
# Stele
# Cortex
# QC/Columella
# Xylem
# Endodermis
#
# Note:
# Some spelling in the source workbook is retained as originally published
# (e.g. "Meristerm", "QC/Colummela").


# =============================================================================
# 18. First annotation check using published cell-type marker lists
# =============================================================================

denyer_markers <- read_excel(
  denyer_xlsm,
  sheet = "Marker Genes"
)


# Convert the wide marker table into long format.

denyer_marker_long <- denyer_markers %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "identity",
    values_to = "gene"
  ) %>%
  dplyr::filter(!is.na(gene))


# Use the top 100 reconstructed markers per cluster.
#
# Top 10 is useful for human inspection but too small for systematic
# marker-overlap analysis.

our_top100 <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = 100,
    with_ties = FALSE
  ) %>%
  dplyr::select(
    cluster,
    gene
  )


# Count shared genes between our cluster marker sets
# and Denyer cell-type reference marker sets.

marker_overlap <- our_top100 %>%
  dplyr::inner_join(
    denyer_marker_long,
    by = "gene",
    relationship = "many-to-many"
  ) %>%
  dplyr::count(
    cluster,
    identity,
    name = "n_overlap"
  )


# Number of reference markers available for each identity.

denyer_marker_counts <- denyer_marker_long %>%
  dplyr::count(
    identity,
    name = "n_reference_markers"
  )


marker_overlap_norm <- marker_overlap %>%
  dplyr::left_join(
    denyer_marker_counts,
    by = "identity"
  ) %>%
  dplyr::mutate(
    recovery = n_overlap / n_reference_markers
  )


marker_overlap_norm %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = recovery,
    n = 3,
    with_ties = TRUE
  ) %>%
  dplyr::arrange(
    cluster,
    dplyr::desc(recovery)
  )


# This classic-marker approach provided useful preliminary identities,
# but many clusters had only a small number of overlaps.
#
# Therefore, a more complete cluster-signature comparison was performed next.


# =============================================================================
# 19. Read complete Denyer C0–C14 cluster signatures
# =============================================================================

# For a fair comparison:
#
# Our cluster      = top 100 marker genes
# Denyer cluster   = top 100 positively enriched genes
#
# Both sides therefore use marker sets of equal size.


denyer_top100 <- lapply(
  0:14,
  function(i) {

    x <- read_excel(
      denyer_xlsm,
      sheet = paste0("C", i)
    )

    names(x) <- c(
      "gene",
      "avg_logFC",
      "p_value",
      "padj",
      "cluster"
    )

    x %>%
      dplyr::filter(
        !is.na(gene),
        avg_logFC > 0,
        padj < 0.05
      ) %>%
      dplyr::arrange(
        dplyr::desc(avg_logFC)
      ) %>%
      dplyr::slice_head(
        n = 100
      ) %>%
      dplyr::transmute(
        denyer_cluster = as.character(i),
        gene
      )
  }
) %>%
  dplyr::bind_rows()


denyer_top100 %>%
  dplyr::count(denyer_cluster)


# =============================================================================
# 20. Compare our 14 clusters with Denyer's 15 clusters
# =============================================================================

stopifnot(all(table(our_top100$cluster) == 100L),
          all(table(denyer_top100$denyer_cluster) == 100L),
          !anyDuplicated(our_top100), !anyDuplicated(denyer_top100))

cluster_overlap <- our_top100 %>%
  dplyr::inner_join(
    denyer_top100,
    by = "gene",
    relationship = "many-to-many"
  ) %>%
  dplyr::count(
    cluster,
    denyer_cluster,
    name = "n_overlap"
  )


# Generate every possible comparison:
#
# 14 reconstructed clusters × 15 Denyer clusters = 210 comparisons.

cluster_scores <- tidyr::crossing(
  cluster = unique(our_top100$cluster),
  denyer_cluster = as.character(0:14)
) %>%
  dplyr::left_join(
    cluster_overlap,
    by = c(
      "cluster",
      "denyer_cluster"
    )
  ) %>%
  dplyr::mutate(
    n_overlap = tidyr::replace_na(
      n_overlap,
      0L
    ),

    # Because both marker sets contain 100 genes:
    #
    # Jaccard =
    # shared genes / total unique genes
    #
    # = overlap / (100 + 100 - overlap)

    jaccard =
      n_overlap /
      (200 - n_overlap)
  )


# Attach Denyer's published biological identity.

denyer_identity <- cluster_overview %>%
  dplyr::filter(
    Cluster %in% as.character(0:14)
  ) %>%
  dplyr::transmute(
    denyer_cluster = Cluster,
    denyer_identity = Identity
  )


cluster_scores <- cluster_scores %>%
  dplyr::left_join(
    denyer_identity,
    by = "denyer_cluster"
  )


# =============================================================================
# 21. Inspect top Denyer matches for each reconstructed cluster
# =============================================================================

cluster_scores %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = n_overlap,
    n = 3,
    with_ties = FALSE
  ) %>%
  dplyr::arrange(
    cluster,
    dplyr::desc(n_overlap)
  ) %>%
  dplyr::select(
    cluster,
    denyer_cluster,
    denyer_identity,
    n_overlap,
    jaccard
  ) %>%
  print(n = 42)


# Strong examples observed:
#
# Our C10 -> Denyer C11 -> QC/Columella
#   51 shared top-100 genes
#
# Our C9 -> Denyer C13 -> Endodermis
#   38 shared genes
#
# Our C1 -> Denyer C10 -> Trichoblast
#   37 shared genes
#
# Our C8 -> Denyer C12 -> Xylem
#   34 shared genes
#
# Our C11 -> Denyer C9 -> Cortex
#   33 shared genes


# =============================================================================
# 22. Visualize complete 14 × 15 similarity matrix
# =============================================================================

heatmap_df <- cluster_scores %>%
  dplyr::mutate(
    cluster = factor(
      cluster,
      levels = as.character(0:13)
    ),

    denyer_cluster = factor(
      denyer_cluster,
      levels = as.character(0:14)
    )
  )


similarity_heatmap <- ggplot(
  heatmap_df,
  aes(
    x = denyer_cluster,
    y = cluster,
    fill = jaccard
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = n_overlap, colour = jaccard < 0.18),
    size = 3
  ) +
  scale_colour_manual(values = c(`FALSE` = "#111111", `TRUE` = "white"),
                      guide = "none") +
  scale_fill_viridis_c() +
  analysis_theme +
  labs(
    title =
      "Cluster similarity between reconstructed and Denyer clusters",
    subtitle =
      "Tile color = Jaccard similarity; number = shared top-100 marker genes",
    x = "Denyer cluster",
    y = "Our reconstructed cluster",
    fill = "Jaccard"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(similarity_heatmap)
save_figure(similarity_heatmap, "07_similarity_heatmap", width = 11, height = 8)

# Interpretation:
#
# Many reconstructed clusters show a clear strongest match to a single
# published Denyer cluster.
#
# Differences are particularly apparent among several "Mature" populations,
# which may explain part of the difference between our 14-cluster solution
# and the published 15-cluster solution.
#
# The analysis does NOT attempt to tune clustering parameters merely
# to reproduce exactly 15 clusters.


# =============================================================================
# 23. Quantify best vs second-best cluster matches
# =============================================================================

cluster_confidence <- cluster_scores %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(
    dplyr::desc(jaccard),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 2) %>%
  dplyr::mutate(
    rank = dplyr::row_number()
  ) %>%
  dplyr::select(
    cluster,
    rank,
    denyer_cluster,
    denyer_identity,
    n_overlap,
    jaccard
  )


cluster_confidence_wide <- cluster_confidence %>%
  tidyr::pivot_wider(
    names_from = rank,
    values_from = c(
      denyer_cluster,
      denyer_identity,
      n_overlap,
      jaccard
    ),
    names_glue = "{.value}_{rank}"
  ) %>%
  dplyr::mutate(
    jaccard_gap =
      jaccard_1 - jaccard_2,

    overlap_gap =
      n_overlap_1 - n_overlap_2
  )


cluster_confidence_wide %>%
  dplyr::arrange(
    dplyr::desc(jaccard_gap)
  ) %>%
  print(n = 14)


# Important distinction:
#
# Cluster-level confidence:
#   How clearly does our cluster match one specific Denyer cluster?
#
# Identity-level confidence:
#   Even if two Denyer clusters are close matches, do they share
#   the same broad biological identity?
#
# Example:
#   Denyer C5 and C10 are both Trichoblast.
#   Therefore, ambiguity between C5 and C10 does not necessarily imply
#   ambiguity in the biological identity "Trichoblast".


# =============================================================================
# 24. Manual biological annotation
# =============================================================================

# IMPORTANT:
#
# The following mapping is NOT produced automatically by Seurat.
#
# It is a manual biological interpretation based on:
#
#   1. reconstructed cluster marker genes
#   2. Denyer classic cell-type markers
#   3. 14 × 15 cluster-signature similarity
#   4. best vs second-best matching patterns
#
# Low-confidence matches are deliberately labeled "-like" rather than
# assigned an overly strong cell identity.


stopifnot(identical(levels(Idents(seurat_qc)), as.character(0:13)))

annotation_table <- dplyr::tibble(
  cluster = factor(
    as.character(0:13),
    levels = as.character(0:13)
  ),

  identity = c(
    "Mature-like",       # C0
    "Trichoblast",       # C1
    "Atrichoblast",      # C2
    "Meristem-like",     # C3
    "Meristem",          # C4
    "Meristem-like",     # C5
    "Stele",             # C6
    "Meristem",          # C7
    "Xylem",             # C8
    "Endodermis",        # C9
    "QC/Columella",      # C10
    "Cortex",            # C11
    "Mature-like",       # C12
    "Mature"             # C13
  ),

  confidence = c(
    "Low",
    "High",
    "High",
    "Moderate",
    "High",
    "Moderate-Low",
    "Moderate",
    "High",
    "High",
    "High",
    "High",
    "High",
    "Low",
    "High"
  )
)


annotation_table


# =============================================================================
# 25. Build a final annotation evidence table
# =============================================================================

final_annotation_table <- cluster_confidence_wide %>%
  dplyr::left_join(
    annotation_table,
    by = "cluster"
  ) %>%
  dplyr::select(
    cluster,

    denyer_cluster_1,
    denyer_identity_1,
    n_overlap_1,
    jaccard_1,

    denyer_cluster_2,
    denyer_identity_2,
    n_overlap_2,
    jaccard_2,

    overlap_gap,
    jaccard_gap,

    identity,
    confidence
  ) %>%
  dplyr::arrange(
    as.numeric(
      as.character(cluster)
    )
  )


final_annotation_table


# This table explicitly separates:
#
# Computational evidence:
#   Denyer matches, overlaps, Jaccard similarity
#
# from:
#
# Biological interpretation:
#   identity and confidence


# =============================================================================
# 26. Add broad identity to Seurat metadata
# =============================================================================

cluster_to_identity <- setNames(
  annotation_table$identity,
  annotation_table$cluster
)


# IMPORTANT:
#
# unname() is required here.
#
# Otherwise Seurat interprets the names "0", "1", ... "13"
# as if they were cell barcodes and reports:
#
# "No cell overlap between new meta data and Seurat object"


seurat_qc$cell_identity <- unname(
  cluster_to_identity[
    as.character(
      Idents(seurat_qc)
    )
  ]
)


# Verify cluster-to-identity assignment.

table(
  Cluster = Idents(seurat_qc),
  Identity = seurat_qc$cell_identity
)


# =============================================================================
# 27. Final annotated UMAP
# =============================================================================

annotated_umap <- DimPlot(
  seurat_qc,
  reduction = "umap",
  group.by = "cell_identity",
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  label.box = FALSE
) +
  analysis_theme +
  labs(
    title =
      "Annotated Arabidopsis root cell populations",
    subtitle =
      "Cell identities inferred by comparison with Denyer et al. cluster signatures"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

print(annotated_umap)
save_figure(annotated_umap, "08_annotated_umap", width = 11, height = 7)

# NOTE:
#
# There are still 14 computational clusters.
#
# However, several clusters share the same broad biological identity:
#
#   C0 + C12 -> Mature-like
#   C3 + C5  -> Meristem-like
#   C4 + C7  -> Meristem
#
# Therefore, the annotated UMAP contains 11 broad identity labels.
#
# No clusters were deleted or merged computationally.
# Only the biological annotation layer was grouped.


###############################################################################
# END OF FIRST-PASS CLUSTERING AND ANNOTATION
#
# Current result:
#
#   4,685 QC-passed cells
#   14 reconstructed computational clusters
#   11 broad biological identity labels
#
# This independently reconstructed cluster structure shows substantial
# correspondence with the published Denyer root atlas.
#
# These identities will next be used to examine where BRL3 × drought
# interaction genes from the independent bulk RNA-seq project are detected
# across normal Arabidopsis root cell populations.
#
# TODO / later robustness analyses:
#
#   - Compare clustering with 20 / 25 / 30 PCs
#   - Explore alternative clustering resolutions
#   - Quantitatively investigate elbow-point selection
#   - Add protoplasting-induced gene sensitivity analysis
###############################################################################

# 28. Persist small evidence tables and reproducibility records.
# Fail rather than apply the fixed manual mapping to an unexpected clustering.
stopifnot(ncol(seurat_qc) == 4685L,
          identical(levels(Idents(seurat_qc)), as.character(0:13)),
          nrow(cluster_scores) == 210L,
          !anyNA(seurat_qc$cell_identity),
          length(unique(seurat_qc$cell_identity)) == 11L)
save_table(data.frame(input_cells = ncol(seurat_obj),
                      retained_cells = ncol(seurat_qc),
                      removed_cells = ncol(seurat_obj) - ncol(seurat_qc),
                      features = nrow(seurat_qc)), "qc_summary")
save_table(as.data.frame(table(cluster = Idents(seurat_qc))), "cluster_sizes")
save_table(top_markers, "top10_markers")
save_table(marker_summary, "marker_summary")
save_table(marker_overlap_norm, "classic_marker_overlap")
save_table(cluster_scores, "denyer_similarity_14x15")
save_table(cluster_confidence_wide, "best_second_evidence")
save_table(annotation_table, "manual_annotation")
save_table(final_annotation_table, "final_annotation_evidence")
save_table(as.data.frame(table(cluster = Idents(seurat_qc),
                              identity = seurat_qc$cell_identity)), "cluster_identity_counts")
writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())),
           "results/reproducibility/sessionInfo.txt")
packages <- c("Seurat", "SeuratObject", "Matrix", "dplyr", "tidyr", "ggplot2",
              "ggrepel", "patchwork", "readxl", "org.At.tair.db", "AnnotationDbi", "uwot")
write.csv(data.frame(package = packages,
                     version = vapply(packages, function(x) as.character(packageVersion(x)), "")),
          "results/reproducibility/package_versions.csv", row.names = FALSE)
print(final_annotation_table, n = 14)
message("Completed: 4,685 cells; 14 computational clusters; 11 provisional broad identities.")
