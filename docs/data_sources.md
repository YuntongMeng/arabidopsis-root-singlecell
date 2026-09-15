# Data sources and verified preparation state

Verified on 2026-09-15. Formal Seurat analysis has not started.

## Primary single-cell matrix

**Local file (not tracked by Git):**
`data/GSE123818_Root_single_cell_wt_datamatrix.csv.gz`

- Source record: NCBI GEO [GSE123818](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123818)
- Direct download: [GSE123818_Root_single_cell_wt_datamatrix.csv.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123818/suppl/GSE123818_Root_single_cell_wt_datamatrix.csv.gz)
- Size: 25,409,960 bytes (24.2 MiB)
- MD5: `96af50339415691a42514800d923ec54` (matches the GEO-hosted file)
- SHA-256: `767a78fda19e9cccbc7791618dbf45657f72d705b3a5ef8e3fbe376b0dcb273a`
- Verified dimensions: 27,629 feature rows × 4,727 cell columns
- Feature identifiers: TAIR locus IDs

This public, easily re-downloadable matrix is excluded from Git. To restore it:

```bash
curl -L \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123818/suppl/GSE123818_Root_single_cell_wt_datamatrix.csv.gz \
  -o data/GSE123818_Root_single_cell_wt_datamatrix.csv.gz
```

GEO links the matrix to Denyer et al. (2019), lists the two WT single-cell
replicates and the same processing framework described in the paper, and reports
processed data under the GSE123818 series.

## Project 1 interaction genes

**Tracked file:** `data/FULL_interaction_DEGs_padj0.05_log2FC1.tsv`

- Origin: Project 1 bulk root RNA-seq genotype × drought interaction results
- Selection represented by the filename: adjusted p-value < 0.05 and absolute
  log2 fold change > 1
- Rows: 137 genes
- SHA-256: `7ea91942d622d8514522c178c87cec65f96f0b9d8866da3a64ca5decb31b2c3e`

Independent verification against the GSE123818 matrix found:

- 137 candidates total
- 136 candidates matched a matrix feature
- 135 candidates had a non-zero count in at least one cell
- `AT5G07985`: no matching matrix feature
- `AT1G66950`: matching feature, but zero in all 4,727 cells

## Denyer et al. cell-cluster annotation framework

**Original publisher supplement (tracked):**
`docs/reference/Denyer2019_TableS2_cluster_DEGs_markers_and_identity.zip`

- Exact source URL: <https://ars.els-cdn.com/content/image/1-s2.0-S1534580719301455-mmc3.zip>
- Publisher object: Supplementary Table 2 (`devcel_4493_mmc3.xlsm` inside the
  downloaded ZIP)
- SHA-256: `1ca359904a62f1b361473892ffb0afa9db60cfbec4e024d51c4491a006aafab9`
- Verified contents: cluster overview, cluster/subcluster differentially
  expressed genes, markers used for cell-identity calling, specificity scores,
  and cluster-level expression summaries

**Convenience extract (tracked):**
`data/reference/Denyer2019_TableS2_cluster_overview.csv`

The original supplement defines 15 clusters (0–14) and gives their cell counts
and broad identities: Mature, Meristem, Atrichoblast, Stele, Trichoblast,
Cortex, QC/Columella, Xylem, and Endodermis. Spelling in the convenience extract
is preserved exactly from the source. Counts sum to 4,727 cells.

### Important annotation limitation

The original GSE123818 supplementary matrix contains barcodes but no cluster
column. The original paper supplement provides the authoritative cluster-level
identity and marker framework, but it does **not** provide a clean per-cell table
mapping each of the 4,727 barcodes to a published cluster. The author-hosted
[Plant scRNA-seq Browser tutorial](https://www.zmbp-resources.uni-tuebingen.de/timmermans/plant-single-cell-browser-root-atlas/)
confirms the Denyer atlas/cluster interpretation but likewise does not expose an
original downloadable barcode-to-cluster metadata table.

No later reprocessed atlas has been substituted and no per-cell labels have been
invented. During formal analysis, cluster assignments must either be reproduced
from the stated Denyer workflow and reconciled with Table S2 markers/identities,
or an original author-provided barcode mapping must be obtained and separately
verified.

## Denyer et al. protoplasting reference

**Original publisher supplement (tracked):**
`docs/reference/Denyer2019_TableS1_RNAseq_metrics_and_protoplasting_DEGs.xlsx`

- Exact source URL: <https://ars.els-cdn.com/content/image/1-s2.0-S1534580719301455-mmc2.xlsx>
- Publisher object: Supplementary Table 1
- SHA-256: `e17d85b2cc7c12dde8914f920e20ffcdc3ac40b063928a5c15bca29a522bcd60`
- Verified contents: sequencing metrics and 6,063 genes differentially expressed
  between protoplasted and un-protoplasted bulk root tissue, including log2FC,
  p-value, and q-value

**Induced-gene subset (tracked):**
`data/reference/Denyer2019_TableS1_protoplasting_induced_FCgt2_q_lt_0.05.tsv`

- Rule: `Log2 FC > 1` (equivalent to fold change > 2) and `q value < 0.05`
- Direction: positive log2FC, i.e. induced in protoplasted relative to
  un-protoplasted root tissue
- Result: 3,545 genes
- SHA-256: `69baaa6d95003b25b29e41906a48aaff27a95748060415bb7112ea588ffbdb84`

Table S1 contains both positive and negative differentially expressed genes;
the 3,545-row file is the explicitly induced subset relevant to the paper's
protoplasting caution. These genes will be added as a flag to candidate-level
and cell-population summaries. They will **not** be automatically excluded,
because this project is contextualizing a biologically defined bulk candidate
set and should show when an interpretation may be sensitive to protoplasting.

## Primary citation

Denyer T, Ma X, Klesen S, Scacchi E, Nieselt K, Timmermans MCP. 2019.
*Spatiotemporal Developmental Trajectories in the Arabidopsis Root Revealed
Using High-Throughput Single-Cell RNA Sequencing.* Developmental Cell 48:
840–852.e5. DOI: [10.1016/j.devcel.2019.02.022](https://doi.org/10.1016/j.devcel.2019.02.022).
PMID: [30913408](https://pubmed.ncbi.nlm.nih.gov/30913408/).
