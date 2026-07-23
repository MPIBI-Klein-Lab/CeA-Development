# CeA-Development

This repository contains the scripts used to reproduce the analyses presented in:

> **Developmental trajectories organize the assembly of valence circuits in the central amygdala**

The repository is organized by manuscript figures. Each figure folder contains plotting scripts and, where applicable, dedicated analysis scripts used to generate intermediate datasets and supplementary tables.

---

# Repository Structure

```
CeA-Development/
│
├── Figure_1/
├── Figure_2/
├── Figure_3/
├── Figure_6/
├── data/
└── README.md
```

**Analysis scripts** perform computational analyses and generate intermediate datasets.

**Figure scripts** reproduce the manuscript figures directly from archived intermediate results.

---

# Contents

## Figure 1

Early developmental specification of central amygdala neurons and reconstruction of developmental trajectories.

---

### `01_Fig1.R`

**Input**

- `CeA.dev.rds`
- `Amygdala.dev.rds`
- `urd.tree.rds`

**Analysis**

- Reproduces Figure 1 and associated Extended Data figures
- Differential expression analysis at E15
- Manhattan distance analysis (coming soon)

---

### `02_Build.URD.Tree.R`

**Input**

- `CeA.dev.rds`

**Analysis**

- Step-by-step reconstruction of the URD developmental tree

---

### `03_P21.Adult.Mapping.R`

**Input**

- `P21_adult.integrated.rds`

**Analysis**

- Adult-informed annotation of P21 cells
- Gene-expression correlation between adult and P21 cell types (coming soon)
- Comparison of developmental trajectory annotations with adult-informed annotations (coming soon)

---

# Figure 2

Divergent germinal zone origins of valence-specific CeA neurons: lineage analysis and ganglionic eminences classification

---

### `01_Fig2.STICR.CeA.Integration.R`

**Input**

- `STICR.CeA.Integrated.rds`

**Analysis**

- Annotation and validation of CeA cell types in the STICR-CeA integrated dataset
- Reproduces Figure 2 and associated Extended Data figures

---

# Figure 3

Gene-family decoding identifies molecular determinants associated with appetitive and aversive developmental trajectories.

---

### `MetaNeighbor.TimeSeries.source.R`

Source file defining the functions used for gene-family based trajectory decoding.

Functions:

- `GeneSetTest()`
    - Compare functional gene sets against randomly sampled gene sets.

- `GeneNumberTest()`
    - Evaluate the effect of gene-set size on decoding performance.

- `GeneSetTestMatrix()`
    - Screen multiple gene families and calculate ΔAUROC.

- `GeneSetTestMatrix2()`
    - Screen multiple gene families and calculate both ΔAUROC and empirical P values.

---

### `01_Fig3.R`

**Input**

- `CeA.dev.rds`
- `Manually.Curated.Rinput.xlsx`
- `Table4.dAUROC.Trajectory.csv`

**Analysis**

- Reproduces Figure 3 and associated Extended Data figures
- Gene-family MetaNeighbor decoding
- Gene-family screening
- Gene-expression visualization

---

### `02_Metaneighbor.Trajectory.Decoding.R`

**Input**

- `CeA.dev.rds`
- `Manually.Curated.Rinput.xlsx`
- `HGNC.Rinput.xlsx`

**Output**

- `Table4.dAUROC.Trajectory.csv`

**Analysis**

- Time-course MetaNeighbor decoding using representative gene families
- Gene-family size analysis
- Genome-wide gene-family screening

---

# Figure 6

Cross-species integration of human and mouse central amygdala datasets and identification of conserved molecular programs.

---

### `01_Fig6.R`

**Input**

- `Human.MouseAdult.Integrated.rds`
- `Human.CeA.consensus.rds`
- `Human.MouseDev.Merged.rds`
- `Table5.dAUROC.CrossSpecies.csv`

**Analysis**

- Human–mouse adult cell-type association analysis
- Human–mouse developmental trajectory comparison using MetaNeighbor
- Reproduces Figure 6 and associated Extended Data figures

---

### `02_Metaneighbor.CrossSpecies.Decoding.R`

**Input**

- `Human.MouseDev.Merged.rds`
- `Human.CeA.consensus.rds`
- `Manually.Curated.Rinput.xlsx`
- `HGNC.Rinput.xlsx`

**Output**

- `Table5.dAUROC.CrossSpecies.csv`

**Analysis**

- Cross-species MetaNeighbor decoding
- Gene-family size analysis
- Representative gene-family analysis using empirical randomization-based null distributions
- Genome-wide screening of conserved gene-family programs

---

# System Requirements

| Software | Version |
|-----------|---------|
| R | 4.x |
| Seurat | 5.3 |
| Harmony | latest |
| URD | latest |
| MetaNeighbor | latest |
| SingleCellNet | latest |
| SCENIC | latest |

Recommended memory: **≥ 32 GB RAM**

---

# Getting Started

1. Clone this GitHub repository.

2. Download the primary datasets from GEO.

3. Download processed integration datasets from Zenodo.

4. Place all files into the `data/` directory.

5. Open the project in RStudio or VS Code.

6. Run

```R
Figure_1/01_Fig1.R
```

to verify the installation.

---

# Data Requirements

| File | Source |
|------|--------|
| `Amygdala.dev.rds` | GEO: GSE334812 (available soon) |
| `CeA.dev.rds` | GEO: GSE334812 (available soon) |
| `urd.tree.rds` | Zenodo (available soon) |
| `P21_adult.integrated.rds` | Zenodo (available soon) |
| `STICR.CeA.Integrated.rds` | Zenodo (available soon) |
| `Human.MouseDev.Merged.rds` | Zenodo (available soon) |
| `Human.CeA.consensus.rds` | Zenodo (available soon) |
| `Human.MouseAdult.Integrated.rds` | Zenodo (available soon) |

---

# Cross-reference Dataset Provenance

| Resource | Source | DOI / Accession |
|----------|--------|-----------------|
| CeA developmental scRNA-seq | This study | GEO: GSE334812 |
| Adult CeA reference | Peters *et al.*, 2023 | GEO: GSE231790 |
| STICR reference atlas | Bandler *et al.*, 2022 | GEO: GSE188528 |
| Human CeA | Human Brain Cell Atlas | CELLxGENE collection |
| Human–mouse integration | This repository | Zenodo DOI (coming soon) |
| P21–Adult integration | This repository | Zenodo DOI (coming soon) |
| CeA–STICR integration | This repository | Zenodo DOI (coming soon) |

---

# Citation

If you use this repository or the accompanying datasets, please cite:

#> He S. *et al.*  
#> **Developmental trajectories organize the assembly of valence circuits in the central amygdala.**
#> bioRxiv (2026). DOI: (coming soon)
