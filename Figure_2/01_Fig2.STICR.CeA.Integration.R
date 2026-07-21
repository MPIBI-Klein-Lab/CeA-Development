library(Seurat)
library(patchwork)
library('harmony')
library(RColorBrewer)
library('metap')
library("dplyr")
library(ggplot2)

################################################################################ read the STICR-CeA integrated dataset
STICR.scCeA.inte <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.2/STICR.CeA.Integrated.rds")

################################################################################ Checking the alignment of CeA-cell types between developmental dataset and STICR-CeA integrated clusters
#### Integrated clusters were assigned as corresponding CeA types when more than 90% of cells originating from the annotated CeA developmental amygdala dataset belongs to the same integrated cluster (majority-vote assignment) 

#### Checking integrated clusters
DimPlot(STICR.scCeA.inte, group.by = "integrated.annotation.numbers", reduction = "tsne", label = T) 
#### Calculating alignment
table(STICR.scCeA.inte$integrated.annotation.numbers, STICR.scCeA.inte$annotation.updated)
tab <- table(STICR.scCeA.inte$annotation.updated, STICR.scCeA.inte$integrated.annotation.numbers)
####column-wise proportions
tab_prop <- prop.table(tab, margin = 2)

# 19, 16, 10, 26, 27, 22, 25, 13 can be assgined >90% onto a specific CeA related cell type
# 28 ambiguous between CeA aversvie and Ast

################################################################################ Load annotations
Idents(STICR.scCeA.inte) <- STICR.scCeA.inte$annotation.allen.final.2025

#### Assign cluster colors
inter.colors <- c("#8DD3C7", "#C3E8BD", "#FAFDB3", "#E2E1C3", "#C3C0D6", "#D5A3B1", "#F3877F", "#D09193", "#95A8C2", "#A6B1B0", "#E1B37A", "#E9BE63", 
                  "#C6D367", "#C2DA83", "#E5D2BF", "#F5CFE2", "#E5D4DD", "#D5CDD5", "#C7A2C7", "#BD89BD", "#C5BCC1", "#CEEBC1", "#E41A1C","#377EB8")
names(inter.colors) <- levels(STICR.scCeA.inte$annotation.allen.final.2025)

################################################################################ Fig. S5 A
DimPlot(STICR.scCeA.inte, group.by = "annotation.allen.final.2025", label = F, reduction = "tsne", na.value = "gray90", 
        cols = inter.colors,
        alpha = 0.8,
        pt.size = 0.5,
        label.size = 10)

### check expression of CeA functional- and spacial- specific markers
FeaturePlot(STICR.scCeA.inte, features = c("Sox5", "Chodl", "Calcrl", "Dlk1", "Prkcd"), reduction = "tsne")
FeaturePlot(STICR.scCeA.inte, features = c("Sst", "Nts", "Il1rapl2", "Tafa1", "Pnoc"), reduction = "tsne")

################################################################################ Fig. 2 A
DimPlot(STICR.scCeA.inte, group.by = "region", label = F, reduction = "tsne", na.value = "gray90", 
        cols = c("#8FABD3", "#D45B97"),
        alpha = 0.8,
        pt.size = 0.5)

################################################################################ Fig. 2 B Violin plot comparing CeA and Bandler dataset
# Subset only CeA aversive and appetitive populations
inter.CeA <- subset(STICR.scCeA.inte, annotation.allen.final.2025 %in% c("CeA_Aversive", "CeA_Appetitive"))
This.CeA <- subset(inter.CeA, region == "CeA")
Bandler.CeA <- subset(inter.CeA, region == "forebrain")

#selected expression, SCT data layer works the best
VlnPlot(This.CeA, features = c("Sox5", "Sox2", "Penk", "Cyp26b1", "Igfbp5",  "Ccnd2", "Plekhh2", "Pax6", "Adora2a","Prkcd", 
                               "Tshz2", "Pcdh7", "Kcnb2", "Sst",  "Isl1", "Ebf1", "Sema6a", "Nefl", "Pnoc", "Nts","S100a10", "Nefm"), 
        stack = T, fill.by = "ident", flip = F, assay = "SCT", layer = "data", cols = c("#E41A1C","#377EB8")) +  theme(
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # X-axis labels larger & angled
          axis.text.y = element_text(size = 20),  # Y-axis labels larger
          axis.title.x = element_text(size = 12, face = "bold"),  # X-axis title larger & bold
          axis.title.y = element_text(size = 12, face = "bold"),
          strip.text = element_text(size = 16)
        )
VlnPlot(Bandler.CeA, features = c("Sox5", "Sox2", "Penk", "Cyp26b1", "Igfbp5",  "Ccnd2", "Plekhh2", "Pax6", "Adora2a","Prkcd", 
                                  "Tshz2", "Pcdh7", "Kcnb2", "Sst",  "Isl1", "Ebf1", "Sema6a", "Nefl", "Pnoc", "Nts","S100a10", "Nefm"), 
        stack = T, fill.by = "ident", flip = F, assay = "SCT", layer = "data", cols = c("#E41A1C","#377EB8")) +  theme(
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # X-axis labels larger & angled
          axis.text.y = element_text(size = 20),  # Y-axis labels larger
          axis.title.x = element_text(size = 12, face = "bold"),  # X-axis title larger & bold
          axis.title.y = element_text(size = 12, face = "bold"),
          strip.text = element_text(size = 16)
        )



################################################################################ Fig. S5 B
##### Find DE genes (using RNA data assay):
DefaultAssay(STICR.scCeA.inte) <- "RNA"
inte.markrs.rna.wilcox <- FindAllMarkers(STICR.scCeA.inte, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.3, test.use = "wilcox", slot = "data")
inte.markrs.rna.wilcox %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20.rna.wilcox

##### or import from saved file:
#top20.rna.wilcox <- read.csv("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/Bandler.CeA.markers.top20.csv")

#### Visualization of expression
agg.rna.data <- AggregateExpression(STICR.scCeA.inte, assays = "RNA", group.by = "ident", slot = "counts", return.seurat = TRUE) 
DoHeatmap(agg.rna.data, features = top20.rna.wilcox$gene, group.colors = inter.colors, assay = "RNA", size = 6, raster = F, draw.lines = F) + 
  scale_fill_gradientn(colors=c("white", "gray90", "gray30"))

