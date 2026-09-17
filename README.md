# Arabidopsis root single-cell contextualization

This project asks where genes with a **BRL3 genotype × drought interaction** in
bulk Arabidopsis root RNA-seq are normally detected across root cell
populations.

## Connection to Project 1

Project 1 identified 137 genes whose drought response differs between BRL3
overexpression and wild type. This project carries that candidate set into the
independent Denyer et al. (2019) wild-type root single-cell atlas (GSE123818) to
add cell-population context to the bulk signal.

## Prepared inputs

- GSE123818 wild-type count matrix: 27,629 TAIR features × 4,727 cells
  (downloaded locally, intentionally not tracked in Git)
- Project 1 interaction-gene table: 137 genes
- Original Denyer et al. Supplementary Tables S1 and S2
- A 3,545-gene protoplasting-induced flag set derived from Table S1 using
  log2FC > 1 and q < 0.05
- The original 15-cluster identity/marker framework from Table S2

The candidate overlap has been verified: 136/137 genes are present in the
matrix and 135/137 are detected in at least one cell. `AT5G07985` is absent;
`AT1G66950` is present but all-zero.

Exact provenance, checksums, download instructions, and the annotation
limitation are recorded in [docs/data_sources.md](docs/data_sources.md).

## Completed first-pass analysis (2026-09-16)

The reproducible workflow now runs **QC → normalization → 2,000 variable genes
→ PCA → clustering → markers → Denyer 14×15 comparison → manual annotation**.
QC retains 4,685 of 4,727 cells (42 extreme high-count/high-feature outliers
removed; these are potential multiplets, not confirmed doublets). Filtering uses
`nCount_RNA < 220000` and `nFeature_RNA < 11000`; no hard organelle-percentage
cutoff is applied. LogNormalize uses scale factor 10,000, neighbors/UMAP use
PCs 1–25, and clustering uses resolution 0.5 with `set.seed(42)`, PCA/UMAP seed 42 and clustering seed 0.
**Exactly 15 clusters were not forced**: the independent result contains 14
computational clusters (0–13).

Positive cluster markers (`min.pct = 0.25`, `logfc.threshold = 0.25`) are ranked
by log2 fold change. Each cluster's top 100 markers is compared with Denyer
C0–C14's top 100 positive significant DEGs (adjusted p-value < 0.05).
All 210 overlaps and Jaccard scores, plus best/second-best matches and gaps,
are calculated before interpretation and saved as evidence tables. The script verifies that both sets have 100
unique genes. Equal scores are ordered by the reference-cluster input order;
a zero gap is evidence of a tie, not a decisive match.

The **11 broad identities** are provisional manual interpretations, including
conservative Mature-like/Meristem-like labels. Manual confidence categories
are qualitative judgments, not calibrated probabilities or algorithmic scores.
C0/C12 share Mature-like, C3/C5 share Meristem-like, and C4/C7 share Meristem;
the original 14 cluster assignments remain unchanged. The fixed manual mapping
is guarded against unexpected cluster IDs and should be reviewed if analysis
parameters or package versions change.

### Run and outputs

From the project root, with the matrix restored as described in
[docs/data_sources.md](docs/data_sources.md):

```r
source("scripts/01_qc_clustering_annotation.R")
```

Or run `Rscript scripts/01_qc_clustering_annotation.R` from a terminal.
Required packages: Seurat, Matrix, dplyr, tidyr, ggplot2, readxl,
org.At.tair.db, AnnotationDbi, and Seurat's plotting/UMAP dependencies.
[Package versions](results/reproducibility/package_versions.csv) and
[sessionInfo](results/reproducibility/sessionInfo.txt) record the verified environment.
The earlier `01_qc_preprocessing.R` is only an import prototype; use the complete
workflow above for analysis.

- [Figures](results/figures): eight final figures, each in PNG and PDF, covering
  QC violin/scatter plots, PCA elbow, cluster UMAP, the complete similarity
  heatmap, and annotated UMAP. Every final plot is assigned and explicitly
  printed in Source mode, with a shared classic 12-point white theme and bold
  titles. UMAP labels use small repelled text without boxes.
- [Summary tables](results/tables): QC counts, cluster sizes, top markers,
  classic-marker overlaps, all 210 similarities, best/second evidence,
  manual annotation, final evidence, and cluster-to-identity counts.
- Large expression matrices, serialized objects, extracted reference workbook,
  and RStudio session files remain ignored; final figures and small tables are tracked.

Candidate-gene population summaries, PC/resolution sensitivity checks,
protoplasting sensitivity analysis, and quantitative elbow exploration remain
future work.

## Interpretation boundary

This is an untreated, normal wild-type atlas. It can show which root cell
populations normally express the bulk interaction candidates, but it cannot
demonstrate that BRL3 changes the drought response within a particular cell
type. That would require genotype- and drought-resolved single-cell data.

## Citation

Denyer T, Ma X, Klesen S, Scacchi E, Nieselt K, Timmermans MCP. 2019.
*Spatiotemporal Developmental Trajectories in the Arabidopsis Root Revealed
Using High-Throughput Single-Cell RNA Sequencing.* Developmental Cell 48:
840–852.e5. https://doi.org/10.1016/j.devcel.2019.02.022
