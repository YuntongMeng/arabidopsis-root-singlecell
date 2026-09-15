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

## Planned analysis (not yet started)

The next stage will import the matrix into Seurat, perform quality control and
dimensional reduction, reproduce or reconcile the published cluster structure,
apply the published marker/identity framework, map the 137 candidates, and
summarize their detection and expression by annotated population.
Protoplasting-induced genes will be flagged for sensitivity and interpretation,
not automatically excluded.

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
