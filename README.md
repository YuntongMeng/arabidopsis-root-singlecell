# Arabidopsis root single-cell contextualization

This project asks where genes with a **BRL3 genotype × drought interaction**
from Project 1 are normally detected across Arabidopsis root cell populations.
It uses the untreated wild-type root single-cell atlas from Denyer et al.
(2019; GSE123818) to provide cellular context for the bulk RNA-seq signal.

## Analysis overview

The final workflow contains four scripts, which are intended to be run in
numerical order from the project root:

| Script | Purpose | Main validated result |
|---|---|---|
| `scripts/00_recover_replicate_metadata.R` | Recovers Denyer WT replicate labels and records direct versus suffix-inferred evidence | 4,727/4,727 cells assigned; 4,458 direct matches; 269 suffix-inferred; 0 ambiguous |
| `scripts/01_qc_clustering_annotation.R` | Performs QC, normalization, dimensionality reduction, clustering, marker comparison, and provisional broad annotation | 4,685 cells; 14 computational clusters; 11 broad identities |
| `scripts/02_interaction_gene_cell_context.R` | Adds population-level context for Project 1 interaction genes and summarizes seven selected GO pathways | 135 detected genes; 1,485 gene-population summaries; 56 exploratory strong candidates; 7 pathway contexts |
| `scripts/03_robustness_checks.R` | Checks replicate mixing, marker-cutoff stability, and protoplasting sensitivity | 4,685/4,685 cells labelled; 14/14 marker matches stable; expected protoplasting overlaps reproduced |

Scripts 01 and 02 keep the previously selected clustering, annotation, and
candidate thresholds unchanged. Script 03 is a sensitivity workflow and does
not redefine the main analysis.

## Inputs and provenance

The main inputs are:

- GSE123818 wild-type count matrix: 27,629 TAIR features × 4,727 cells.
- Project 1 interaction-gene table: 137 genes.
- Denyer et al. Supplementary Tables S1 and S2.
- A 3,545-gene protoplasting-induced reference set defined by
  `log2FC > 1` and `q < 0.05` in Denyer Table S1.
- Shahan et al. (2022) Data S3, used only to validate and recover the original
  Denyer WT replicate identities.

The public count matrix and the large Shahan workbook are intentionally not
tracked in Git. Their exact provenance, checksums, and preparation notes are
recorded in [`docs/data_sources.md`](docs/data_sources.md) and
[`data/reference/Denyer2019_WT_replicate_recovery_report.md`](data/reference/Denyer2019_WT_replicate_recovery_report.md).

To restore the count matrix:

```bash
curl -L \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123818/suppl/GSE123818_Root_single_cell_wt_datamatrix.csv.gz \
  -o data/GSE123818_Root_single_cell_wt_datamatrix.csv.gz
```

To regenerate the replicate map with script 00, download the Shahan et al.
Data S3 archive from the publisher, extract `SuppData3_complete.xltx`, and save
it as:

```text
data/reference/Shahan2022_DataS3_WT_atlas_metadata.xltx
```

Publisher archive:
<https://ars.els-cdn.com/content/image/1-s2.0-S1534580722000338-mmc3.zip>

## Running the complete workflow

Start R in the project root and run:

```r
source("scripts/00_recover_replicate_metadata.R")
source("scripts/01_qc_clustering_annotation.R")
source("scripts/02_interaction_gene_cell_context.R")
source("scripts/03_robustness_checks.R")
```

The same scripts can be run separately with `Rscript`. Scripts 02 and 03
reconstruct their required upstream objects when run from a clean R session;
therefore, they do not depend on an undocumented interactive workspace.

Core R packages include Seurat, Matrix, dplyr, tidyr, ggplot2, readxl,
AnnotationDbi, and org.At.tair.db. Verified package versions and session
information are stored under [`results/reproducibility`](results/reproducibility).

The message
`'select()' returned 1:many mapping between keys and columns` is expected when
GO annotations are queried because one TAIR gene can map to multiple GO terms.
It is not an error; downstream dimensions and counts are checked explicitly.

## Main analysis

### QC, clustering, and annotation

QC removes 42 extreme high-count/high-feature cells, retaining 4,685 of 4,727
cells. Filtering uses `nCount_RNA < 220000` and `nFeature_RNA < 11000`; no hard
organelle-percentage cutoff is applied. The workflow uses LogNormalize with a
scale factor of 10,000, 2,000 variable genes, PCs 1–25, clustering resolution
0.5, and fixed random seeds. Exactly 15 clusters are not forced: the independent
analysis produces 14 computational clusters.

Positive cluster markers (`min.pct = 0.25`, `logfc.threshold = 0.25`) are ranked
by log2 fold change. Each reconstructed cluster's top 100 markers is compared
with the top 100 positive significant genes from each Denyer C0–C14 signature.
The complete 14 × 15 overlap and Jaccard matrix is calculated before biological
interpretation. The resulting 11 broad identities are provisional manual
summaries; the original 14 reconstructed cluster assignments remain unchanged.

### Interaction-gene cellular context

Of the 137 Project 1 interaction genes, 136 are present in the atlas feature
list and 135 are detected in at least one cell. `AT5G07985` is absent;
`AT1G66950` is present but all-zero. Expression is summarized for every detected
gene across 11 broad populations, yielding 135 × 11 = 1,485 checked
gene-population combinations.

Figure 09 displays a **per-gene z-score across populations**, separated into
Positive and Negative interaction blocks. Blue means relatively higher
expression for that gene and red means relatively lower expression. These are
relative-expression values, not enrichment/depletion statistics, and absolute
color values should not be compared between genes.

For each gene, `delta_z` is the difference between its highest and second-highest
population z-scores. The exploratory strong-candidate definition is
`delta_z >= 1` plus detection in at least 10% of cells in the preferred
population; 56 genes meet both criteria. In Figure 10, Negative interactions
are circles and Positive interactions are triangles.

Figure 11 summarizes the preferred-population distribution of genes in seven
previously selected positive-interaction GO terms. It is a descriptive pathway
cellular-context view, not a new pathway-enrichment test.

## Robustness checks

### Replicate identity and mixing

The replicate crosswalk contains 4,727 unique full barcodes:

- 4,458 assignments (94.31%) are directly cross-validated against Shahan Data S3.
- 269 assignments (5.69%) are suffix-inferred after the crosswalk is validated
  with zero discordant direct matches.
- There are 0 ambiguous sample-aware mappings and 0 duplicated full barcodes.
- All 4,685 QC-passed cells receive a replicate label, with both replicates
  represented in all 14 clusters and all 11 broad populations.

The supported interpretation is **no obvious replicate-driven segregation and
broad mixing**. This visual/compositional sanity check does not establish the
complete absence of batch effects and is not replicate-level differential
testing.

### Marker-cutoff sensitivity

The existing annotation evidence is compared fairly at three equal-size
cutoffs: reconstructed top 50 versus Denyer top 50, top 100 versus top 100, and
top 200 versus top 200. Jaccard is calculated as
`overlap / (2 × cutoff - overlap)`. All 14 reconstructed clusters retain the
same best Denyer match at all three cutoffs. This supports the stability of the
existing matching and does not redefine annotation.

### Protoplasting sensitivity

The protoplasting-induced reference flags:

- 40/135 detected interaction genes (29.6%).
- 20/56 exploratory strong candidates (35.7%).
- 13/80 Negative interaction genes (16.2%).
- 27/55 Positive interaction genes (49.1%).

Protoplasting is therefore an important confounder, particularly for positive
interaction and stress-related genes. The `cold acclimation` term is fully
flagged (3/3). After flagged genes are excluded, the other surviving GO terms
often contain only two genes. Mature/Stele context is not completely removed,
but this small-n, overlapping-gene result is sensitivity evidence rather than
strong independent confirmation.

## Outputs

- [`results/figures`](results/figures): QC, clustering, annotation,
  interaction-gene, pathway-context, and replicate UMAP figures (Figures
  01–12; PNG/PDF where applicable).
- [`results/tables`](results/tables): annotation evidence, population summaries,
  candidate tables, pathway contexts, replicate composition, marker-cutoff
  sensitivity, and protoplasting overlap tables.
- [`results/reproducibility`](results/reproducibility): validation summaries,
  package versions, and session information.

The principal robustness audit is
[`results/reproducibility/robustness_validation.txt`](results/reproducibility/robustness_validation.txt).
It records the expected counts and the interpretation boundaries for direct
cross-validation, suffix inference, replicate mixing, marker stability, and
protoplasting sensitivity.

## Interpretation boundary

This is an untreated, normal wild-type atlas. It identifies cell populations
that normally express the bulk interaction candidates and adds cellular context
to the Project 1 result. It cannot demonstrate that the BRL3 × drought effect
occurs within a preferred population; that would require genotype- and
drought-resolved single-cell data.

## Citation

Denyer T, Ma X, Klesen S, Scacchi E, Nieselt K, Timmermans MCP. 2019.
*Spatiotemporal Developmental Trajectories in the Arabidopsis Root Revealed
Using High-Throughput Single-Cell RNA Sequencing.* Developmental Cell 48:
840–852.e5. <https://doi.org/10.1016/j.devcel.2019.02.022>
