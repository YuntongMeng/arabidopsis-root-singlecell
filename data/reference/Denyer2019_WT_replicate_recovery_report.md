# Denyer 2019 WT replicate identity recovery

## Result

The 4,727 columns in `GSE123818_Root_single_cell_wt_datamatrix.csv.gz`
can be assigned to the two Denyer WT replicates without re-running the raw
scRNA-seq pipeline.

- Matrix suffix `-1` maps to `dc1`, `GSM3511858`, `rep1`.
- Matrix suffix `-2` maps to `dc2`, `GSM3511859`, `rep2`.
- 4,458/4,727 cells (94.31%) have a direct, sample-aware barcode match in
  Shahan 2022 Data S3.
- Direct matches: 2,224 `dc1` cells and 2,234 `dc2` cells.
- 269 cells are absent from the final WT atlas metadata: 143 `-1` cells and
  126 `-2` cells. These are retained in the GSE123818 processed matrix but
  were not included in the final Shahan WT atlas metadata.
- Across all 4,458 direct overlaps, the matrix suffix and Data S3 sample label
  agree without exception. This validates applying the same suffix crosswalk
  to the 269 cells absent from the final atlas metadata.
- Final assignments: 2,367 `rep1` cells and 2,360 `rep2` cells.
- There are no duplicated full matrix barcodes and no ambiguous sample-aware
  mappings. Seven 16-nt base barcode sequences occur once in each library
  (14 cells total), which is expected because 10x barcode sequences can be
  reused across independent libraries. The `-1`/`-2` suffix must therefore be
  retained for joins.

## Source inspection

The publisher copy of Shahan 2022 Data S3 is stored as
`Shahan2022_DataS3_WT_atlas_metadata.xltx`. It contains two sheets:

- SHA-256: `fc2ce681fa252d0ce08907acc3b55656b03c15573b7dd6253ffefab1f58142a1`
- The 62 MiB workbook is a public regeneration input and is intentionally not
  tracked in Git; the derived barcode-to-replicate table is tracked.

1. `atlas meta data`: 110,427 cells and 60 metadata columns.
2. `Legend for column name`: definitions for the metadata columns.

The fields used for recovery are `cell_barcode_id` and `orig.ident`. Data S3
contains 3,378 `dc1` cells, all ending in `_15`, and 3,283 `dc2` cells, all
ending in `_16`. The merged GSE123818 matrix contains 2,367 barcodes ending in
`-1` and 2,360 ending in `-2`. A sample-aware normalization
(`dc1: _15 -> -1`; `dc2: _16 -> -2`) gives the direct overlaps above.

Sources:

- Shahan et al. 2022 article and Data S3 description:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC9014886/
- Publisher Data S3 archive used here:
  https://ars.els-cdn.com/content/image/1-s2.0-S1534580722000338-mmc3.zip
- GSE152766 sample crosswalk (`dc1`/`dc2`):
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152766
- Denyer GSE123818 processed matrix and replicate accessions:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123818

## Implemented replicate sanity check

The mapping is consumed by `scripts/03_robustness_checks.R`. It is joined to
the QC-passed Seurat object by exact `original_barcode` without changing the
clustering, annotation, or main-analysis parameters. The formal workflow saves
the replicate UMAP and cluster/population composition tables under `results/`.

```r
rep_map <- read.delim(
  "data/reference/Denyer2019_WT_barcode_to_replicate.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

idx <- match(colnames(seurat_qc), rep_map$original_barcode)
stopifnot(!anyNA(idx), !anyDuplicated(rep_map$original_barcode))

seurat_qc$replicate <- rep_map$replicate[idx]
seurat_qc$source_sample <- rep_map$source_sample[idx]
seurat_qc$replicate_match_status <- rep_map$match_status[idx]

# Visual mixing check used by the formal workflow
Seurat::DimPlot(
  seurat_qc,
  reduction = "umap",
  group.by = "replicate"
)

# Cluster composition
cluster_counts <- with(
  seurat_qc@meta.data,
  table(seurat_clusters, replicate)
)
cluster_proportions <- prop.table(cluster_counts, margin = 1)

# Broad-population composition
population_counts <- with(
  seurat_qc@meta.data,
  table(broad_population, replicate)
)
population_proportions <- prop.table(population_counts, margin = 1)
```

Interpretation should focus on whether both replicates contribute across the
UMAP, clusters, and broad populations, while recognizing that cell counts are
not independent biological replicates for formal differential testing.
