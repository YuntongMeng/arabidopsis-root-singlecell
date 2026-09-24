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

## Completed Project 2 analysis

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

- [Figures](results/figures): ten final figures, each in PNG and PDF, covering
  QC violin/scatter plots, PCA elbow, cluster UMAP, the complete similarity
  heatmap, annotated UMAP, interaction-gene heatmap, and highlight scatter.
  Every final plot is assigned and explicitly printed in Source mode, with a
  shared classic 12-point white theme and bold titles. UMAP labels use small
  repelled text without boxes.
- [Summary tables](results/tables): QC counts, cluster sizes, top markers,
  classic-marker overlaps, all 210 similarities, best/second evidence,
  manual annotation, final evidence, cluster-to-identity counts, population
  expression summaries, preferred populations, preference strength, and
  exploratory strong candidates/counts.
- Large expression matrices, serialized objects, extracted reference workbook,
  and RStudio session files remain ignored; final figures and small tables are tracked.

Project 2 has now progressed from cluster annotation to **Project 1 →
single-cell contextualization**. Run the complete workflow from the project
root with:

```r
source("scripts/01_qc_clustering_annotation.R")
source("scripts/02_interaction_gene_cell_context.R")
```

The second script can also be run directly with
`Rscript scripts/02_interaction_gene_cell_context.R`; it reconstructs the 01
analysis when `seurat_qc` is not already available and therefore does not rely
on an undocumented interactive session.

Of the 137 Project 1 interaction genes, 136 are present among atlas features
and 135 have at least one non-zero count. `AT5G07985` is absent from the
expression matrix, while `AT1G66950` is a matrix feature but is all-zero in
these cells. The analysis calculates detection rate and mean log-normalized
expression for every detected gene in each of the 11 broad populations,
producing an explicitly checked 135 × 11 matrix (1,485 combinations).

The clustered heatmap displays a **per-gene z-score across populations**: zero
is that gene's across-population mean and positive/negative values indicate
relative enrichment/depletion for that gene. It does not compare absolute
expression between different genes. Euclidean distance with ward.D2 linkage is
used only to order heatmap genes; it does not define candidate status.

For each gene, the preferred population is the population with the highest
expression z-score. `delta_z` is the difference between its highest and
second-highest population z-scores, and preferred-population detection rate is
the fraction of that population's cells with a non-zero raw count. The final
highlight scatter uses detection rate on the x-axis and `delta_z` on the
y-axis; shape indicates Positive/Negative Project 1 interaction direction.
Points meeting the exploratory cutoff (`delta_z >= 1` and detection rate
`>= 0.10`) are colored by preferred population, while all other points are
semi-transparent gray. This yields 56 exploratory strong candidates. The
clustered heatmap and highlight scatter are figures 09 and 10; population,
preference, candidate, and count tables are in `results/tables/`.

### Next checkpoint

- pathway-level cellular contextualization
- replicate/batch sanity check
- top-50/100/200 marker-overlap robustness
- protoplasting sensitivity
- add `scripts/00_prepare_inputs.R` and/or `scripts/00_verify_inputs.R`
- final README polish

## Interpretation boundary

This is an untreated, normal wild-type atlas. It can show which root cell
populations normally express the bulk interaction candidates and provide
cellular context, but it cannot demonstrate that a BRL3 × drought effect
occurs within the preferred cell type. That would require genotype- and
drought-resolved single-cell data.

## Citation

Denyer T, Ma X, Klesen S, Scacchi E, Nieselt K, Timmermans MCP. 2019.
*Spatiotemporal Developmental Trajectories in the Arabidopsis Root Revealed
Using High-Throughput Single-Cell RNA Sequencing.* Developmental Cell 48:
840–852.e5. https://doi.org/10.1016/j.devcel.2019.02.022
