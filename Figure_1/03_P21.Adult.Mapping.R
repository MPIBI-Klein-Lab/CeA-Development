# This script documents the integration of P21 developmental cells with an
# adult CeA reference and the assignment of adult-informed P21 annotations.
#
# The finalized P21 annotations used in the manuscript are already stored in
# the `P21.annotated` metadata column of CeA.dev.rds. This script does not
# overwrite that archived object.

library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(sctransform)
library(harmony)
library(RColorBrewer)
library(here)

data_dir <- here("data")

########################### load P21-adult integrated dataset ########################### 
P21_adult.se <- readRDS(file.path(data_dir, "P21_adult.integrated.rds"))

### Inspect P21-adult integrated clusters
DimPlot(P21_adult.se, group.by = "inte.number", label = T)

### Inspect developmental stages
DimPlot(P21_adult.se, group.by = "stage")

### Inspect individual dataset
# Note that the Adult dataset contains single nuclei extracted from both Fasted and fed animals, 
# because these conditions do not affect cell clustering (Peters et al., 2023),
# they are treated equally here.

DimPlot(P21_adult.se, group.by = "batch")

### Inspect adult annotation
DimPlot(P21_adult.se, group.by = "adult.annotation", label = T, na.value = "gray80",
        cols = colorRampPalette(brewer.pal(8, "Set1"))(
          length(levels(as.factor(P21_adult.se$adult.annotation)))),
        alpha = 0.9,
        pt.size = 1,
        repel = T)

########################### Inspect the maximum adult-reference contribution for each integrated cluster ########################### 
tab <- table(P21_adult.se$adult.annotation, P21_adult.se$inte.number)
# column-wise proportions
tab_prop <- prop.table(tab, margin = 2)

# assign cluster name of inte.number directly from adult.annotation if there is a >0.75 column proportion match
mapping.summary <- data.frame(
  inte.number = colnames(tab_prop),
  adult.annotation = apply(
    tab_prop,
    2,
    function(x) rownames(tab_prop)[which.max(x)]
  ),
  maximum.proportion = apply(tab_prop, 2, max)
)
mapping.summary
sum(mapping.summary$maximum.proportion > 0.75)

### Adult identities were assigned directly when one adult population accounted for more than 75% of an integrated cluster. 
### Ambiguous clusters were resolved by inspecting marker-gene expression below:

DefaultAssay(P21_adult.se) <- "integrated"
Idents(P21_adult.se) <- P21_adult.se$inte.number
P21_adult.markers <- FindAllMarkers(P21_adult.se, only.pos = TRUE, min.pct = 0.02, logfc.threshold = 0.25)
P21_adult.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20.rna.v2
DoHeatmap(P21_adult.se, features = top20.rna.v2$gene) + NoLegend() + scale_fill_gradientn(colors=c("#6395c7", "white", "#e06ead"))

########################### Assign adult-informed cluster identities ########################### 
# some updates of adult annotations: old --> new names:
# CeM_Cacna1g --> CeM_Il1rapl2.Tafa1
# CeM_Gfra1.Dlk1 --> CeM_Dlk1
P21_adult.se <- RenameIdents(P21_adult.se, "0" = "CeL_Prkcd",
                               "4" = "CeC_Cdh9.Calcrl",
                               "12" = "CeM_Dlk1",
                               "8" = "IPAC_Pde1c",
                               "2" = "Ast_Synpo2",
                               "3" = "Ast_Dach2",
                               "9" = "BLA_Gbx1.Lhx8",
                               "13" = "BLA_Nkx2-1",
                               "15" = "BLA_Nkx2-1",
                               "7" = "ITC_Foxp2",
                               "14" = "CeM_Tac1.Sst",
                               "11" = "CeM_Drd2.Rai14",
                               "1" = "CeM_Il1rapl2.Tafa1",
                               "10" = "CeM_Vdr",
                               "6" = "CeL_Sst",
                               "5" = "CeL_Nts.Tac2")
DimPlot(P21_adult.se, label = T)

# Optional: save a reconstructed mapping object without overwriting archived data
# saveRDS(
#   P21_adult.se,
#   file.path(data_dir, "P21_adult.integrated.reannotated.rds")
# )

########################### Archived P21 annotations ###########################
# The finalized identities generated from this analysis were transferred to the developmental Seurat object during the original analysis. 
#They are already available in CeA.dev.rds as the metadata column `P21.annotated`:

### read in CeA developmental dataset
CeA.dev <- readRDS(file.path(data_dir, "CeA.dev.rds"))

DimPlot(CeA.dev, group.by = "P21.annotated", label = T, na.value = "gray90",
        cols = colorRampPalette(brewer.pal(8, "Set2"))(
          length(levels(CeA.dev$P21.annotated))),
        alpha = 0.9,
        pt.size = 1,
        repel = T)

########################### Correlation and river plots ###########################
# To be completed by Chao


