# CeA-Development

This repository contains scripts used to reproduce the analyses presented in

Developmental trajectories organize the assembly of valence circuits in the central amygdala

## Contents

#### Figure_1

01_Fig1.R
Input: 
    - CeA.dev.rds  
    - Amygdala.dev.rds  
    - urd.tree.rds 
Analysis: 
    Reproduction of plots in Figure 1 and associated extended data figures; 
    Defferential expression analysis at E15; 
    Manhattan distance analysis (Cao)

02_Build.URD.Tree.R
Input: 
    - CeA.dev.rds
Analysis:
    A Step-by-step build of the URD tree model

03_P21.Adult.Mapping.R
Input: 
    - P21_adult.integrated.rds
Analysis: 
    Adult-informed annotation of P21 cells; 
    gene expression correlation analysis of Adult and P21 annotated cell clusters (Chao);
    comparison between developmental trajecory annotations and adult-informed annotations of P21 cells (Chao)


 
#### Figure_2

01_Fig2.STICR.CeA.Integration.R
Input:
    - STICR.CeA.Integrated.rds
Analysis: 
    Annotation and validationof CeA cell types in the STICR-CeA integrated dataset; 
    Reproduction of plots for lineage analysis in Figure 2 and associated extended data figures.



 
#### Figure_3

MetaNeighbor.TimeSeries.source.R
Source file defined the following functions used for gene family based CeA developmental trajectory decoding:
    GeneSetTest() - Function for making random chance line and compare with functionally-defined gene set
    GeneNumberTest() - Function that test the effect of gene numbers on decoding scores (AUROC)
    GeneSetTestMatrix() - Function that test a list of gene sets, normalize the AUROC score with the chance line (∆AUROC)
    GeneSetTestMatrix2() - Function that test a list of gene sets, return both the ∆AUROC scores and p-value comparing the paired gene set and randomly selected genes


01_Fig3.R
Input:
    - CeA.dev.rds, 
    - Manually.Curated.Rinput.xlsx,
    - Table4.dAUROC.Trajecory.csv
Analysis:
Reproduction of plots for MetaNeighbor decoding, gene family screening, and gene expression in Figure 3 and associated extended data figures.


02_Metaneighbor.Trajectory.Decoding.R
Input:
    - CeA.dev.rds, 
    - Manually.Curated.Rinput.xlsx,
    - HGNC.Rinput.xlsx
Analysis/Output:
    Plots for time-course decoding of CeA trajectories using example gene families 
    ‘Table4.dAUROC.Trajecory.csv’ as one of the outputs

 
#### Figure_6

01_Fig6.R
Inputs
    - Human.MouseAdult.Integrated.rds
    - Human.CeA.consensus.rds
    - Human.MouseDev.Merged.rds 
    - Table5.dAUROC.CrossSpecies.csv
Analysis
    Human–mouse adult cell type association analysis
    Human–mouse developmental cell type association analysis using MetaNeighbor
    Reproduction of plots for Figure 6 and associated extended data figures


02_Metaneighbor.CrossSpecies.Decoding.R
Inputs
    - Human.MouseDev.Merged.rds
    - Human.CeA.consensus.rds
    - Manually.Curated.Rinput.xlsx
    - HGNC.Rinput.xlsx
Analysis/output
Cross-species decoding of CeA consensus cell types using the MetaNeighbor algorithm as a decoder, including:
    Detailed analysis of the gene family size-effects; 
    Example gene family analysis with bootstrapping for their null distributions; 
    Computational screening for conserved gene-family based genetic programs   
    ‘Table5.dAUROC.CrossSpecies.csv’ as one of the outputs



## System Requirements
R 4.x
Seurat 5.3
Harmony
URD
MetaNeighbor
SingleCellNet
SCENIC


## Getting Started
1. Clone the GitHub repository
2. Download the processed GEO objects (GSE334812, available soon)
3. Download the processed .rds files from Zenodo (DOI: ..., available soon)
4. Place all required files into data/.
5. Open the project in R/VS Code.
6. Run 'Figure_1/01_Fig1.R' 


## Data Requirements
| File                      | Source |
| ------------------------- | ------ |
| Amygdala.dev.rds          | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334812 (available soon) |
| CeA.dev.rds               | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334812 (available soon) |
| urd.tree.rds              | Zenodo, available soon |
| P21_adult.integrated.rds  | Zenodo, available soon |  
| STICR.CeA.Integrated.rds  | Zenodo, available soon |
| Human.MouseDev.Merged.rds | Zenodo, available soon |
| Human.CeA.consensus.rds   | Zenodo, available soon |
| Human.MouseAdult.Integrated.rds | Zenodo, will be available soon |


## Cross-reference dataset provenance
Resource	Source	DOI / accession
CeA developmental scRNA-seq	This study	GEO: GSE334812
Adult CeA reference	Peters et al., 2023	GEO: GSE231790
STICR dataset	Bandler et al.,2022	GEO: GSExxxxx
Human CeA	Human Brain Cell Atlas	CELLxGENE collection
Integrated datasets	This repository / Zenodo	DOI: cooming soon

