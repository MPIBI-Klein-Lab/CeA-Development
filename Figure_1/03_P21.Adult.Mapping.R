library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(sctransform)
library(harmony)
library(RColorBrewer)

########################### load P21-adult integrated dataset ########################### 
P21_adult.se <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.1/P21_adult.integrated.rds")

### Inspect P21-adult integrated clusters
DimPlot(P21_adult.se, group.by = "inte.number", label = T)

### Inspect developmental stages
DimPlot(P21_adult.se, group.by = "stage")

### Inspect individual dataset
# Note that the Adult dataset contains single nuclei extracted from both Fasted and fed animals, 
# since these conditions doesn't affect cell clustering (https://www.science.org/doi/10.1126/sciadv.adf6521), these conditions are treated equally here
DimPlot(P21_adult.se, group.by = "batch")

### Inspect adult annotation
DimPlot(P21_adult.se, group.by = "adult.annotation", label = T, na.value = "gray80",
        cols = colorRampPalette(brewer.pal(8, "Set1"))(
          length(levels(as.factor(P21_adult.se$adult.annotation)))),
        alpha = 0.9,
        pt.size = 1,
        repel = T)

########################### label transfer via "majority vote / maximum overlap assignment" ########################### 
tab <- table(P21_adult.se$adult.annotation, P21_adult.se$inte.number)
# column-wise proportions
tab_prop <- prop.table(tab, margin = 2)
# assign cluster name of inte.number directly from adult.annotation if there is a >0.75 column proportion match

###########################  for < 0.75 column proportions, ambiguous mapping, checking markers instead " ########################### 
DefaultAssay(P21_adult.se) <- "integrated"
Idents(P21_adult.se) <- P21_adult.se$inte.number
P21_adult.markers <- FindAllMarkers(P21_adult.se, only.pos = TRUE, min.pct = 0.02, logfc.threshold = 0.25)
P21_adult.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20.rna.v2
DoHeatmap(P21_adult.se, features = top20.rna.v2$gene) + NoLegend() + scale_fill_gradientn(colors=c("#6395c7", "white", "#e06ead"))

########################### rename clusters ########################### 
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
#saveRDS(P21_adult.se, "~/Documents/Dev_manuscript/fig.1_trajectory/P21_adult.integration.rds")

### Update metadata of CeA developmental dataset 'CeA.dev.rds':
#P21.anno <- subset(P21_adult.se, subset = stage == "P21")
#table(Idents(P21.anno))
#CeA.dev <- AddMetaData(CeA.dev, metadata = Idents(P21.anno), col.name = "P21.annotated")
#DimPlot(CeA.dev, group.by = "P21.annotated", label = T, na.value = "gray90",
#        cols = colorRampPalette(brewer.pal(8, "Set2"))(
#          length(levels(CeA.dev$P21.annotated))),
#        alpha = 0.9,
#        pt.size = 1,
#        repel = T)

########################### correlation, river plot (Chao) ########################### 


