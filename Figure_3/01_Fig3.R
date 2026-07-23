########################################################## ######################################################## 
# Figure 3 visualization script
#
# This script reproduces the figure panels from previously generated trajectory-decoding results.
#
# To reproduce the decoding analyses themselves, run:
#   02_Metaneighbor.Trajectory.Decoding.R
#
# The helper functions are implemented in:
#   MetaNeighbor.TimeSeries.source.R
##########################################################  ######################################################## 

library(Seurat)
library(RColorBrewer)
library(svMisc)
library(ggplot2)
library(cowplot)
library(dplyr)
library(rlang)
library(pheatmap)
library(here)

data_dir <- here("data")

########################################################## Fig. 3C
### read mean dAUROC value
dAUROC.table <- read.csv(file = file.path(data_dir, "Table4.dAUROC.Trajectory.csv"), row.names = 1,
                                  check.names = FALSE)

### read manually curated gene family names
library(readxl)
manual.curated <- readxl::read_excel(file.path(data_dir, "Manually.Curated.Rinput.xlsx"))

manual.curated <- lapply(manual.curated, as.character)
names(manual.curated) <- make.names(names(manual.curated))
manual.curated <- lapply(manual.curated, function(x) x[!is.na(x)])


### column metadata
col.meta <- data.frame(
  stage = rep(c("E15", "E18", "P0", "P4", "P10", "P21"), 2),
  trajectory = rep(c("Appetitive", "Aversive"), each = 6)
)
col.meta$stage <- factor(col.meta$stage, levels = c("E15", "E18", "P0", "P4", "P10", "P21"))


### row metadata
is_manual <- rownames(dAUROC.table) %in% names(manual.curated)
source_label <- ifelse(is_manual, "Manually.curated", "HGNC.database")
row.meta <- data.frame(Source = source_label, row.names = rownames(dAUROC.table))


stage_colors <- brewer.pal(6, "PRGn")
names(stage_colors) <- c("E15", "E18", "P0", "P4", "P10", "P21")

my_colour = list(
  trajectory = c(Appetitive = "#377EB8", Aversive = "#E41A1C"),
  stage = stage_colors,
  Source = c(Manually.curated = "coral3", HGNC.database =   "burlywood1"))


### P-values are not displayed across the complete heatmap to preserve figure readability.
### Statistical significance was evaluated separately for the enlarged panels using GeneSetTestMatrix2().
### Asterisks were added manually based on those calculated p-values.

pheatmap(dAUROC.table,
         cluster_rows = T,
         cluster_cols = F,
         #color = colorRampPalette(brewer.pal(n = 9, name ="BuPu"))(100), # good
         #color = colorRampPalette(brewer.pal(n = 9, name ="YlOrBr"))(100), #good
         #color = colorRampPalette(brewer.pal(n = 9, name ="OrRd"))(100), # good
         color = colorRampPalette(brewer.pal(n = 6, name ="RdPu"))(100), # good
         #color = colorRampPalette(brewer.pal(n = 6, name ="PuBu"))(100),
         #cellwidth = 10,
         #cellheight = 10,
         fontsize_col=10,
         fontsize_row=10,
         cutree_rows =7,
         gaps_col = 6,
         annotation_col = col.meta,
         annotation_row = row.meta,
         annotation_colors = my_colour,
         angle_col = "270"
)


########################################################## Prepare data for expression plots ######################################################## 
### load CeA developmental dataset
se.data.all <- readRDS(file.path(data_dir, "CeA.dev.rds"))

### SCT normalization
se.data.list <- SplitObject(se.data.all, split.by = "batch")
se.data.list <- lapply(X = se.data.list, FUN = function(x) {
  SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE)
})
se.data.sct <- merge(se.data.list[[1]], se.data.list[2:length(se.data.list)])
se.data.sct$stage <- factor(se.data.sct$stage, levels = c("E15", "E18", "P0", "P4", "P10", "P21"))

### subset only CeA cell types (tree segments) for expression plots
se.data.sct <- subset(se.data.sct, segments.simple %in% c("P_app", "CeL_appetitive", "CeM_appetitive", 
                                                          "P_ave", "CeC/L_aversive", "CeM_aversive"))

se.data.sct$segments.simple <- factor(se.data.sct$segments.simple, levels = c("P_app", "CeL_appetitive", "CeM_appetitive", 
                                                                              "P_ave", "CeC/L_aversive", "CeM_aversive"))


##################################################### HEATMAP gene expression summary (Fig. 3D related) ######################################################## 

### Aggregate SCT-normalized expression across developmental subpopulations for heatmap visualization
DefaultAssay(se.data.sct) <- "SCT"
agg.sct <- AggregateExpression(se.data.sct, assays = "SCT", group.by = c("segments.simple", "stage"), return.seurat = T)
agg.sct$segments.simple <- factor(agg.sct$segments.simple, levels = c("P-app", "CeL-appetitive", "CeM-appetitive", 
                                                                      "P-ave", "CeC/L-aversive", "CeM-aversive"))
### Gene expression heatmap
DoHeatmap(agg.sct, 
          slot = "scale.data",
          draw.lines = FALSE,
          #disp.min = -0.5,
          #disp.max = 1,
          raster = F,
          features = c("Sema3a", "Sema6a", "Sema3c", "Sema5a","Sema3e",
                       
                       "Gal", "Rln1", "Pnoc","Nts",  "Sst", "Tac2", "Crh", "Pdyn", "Nmb",  "Tac1", "Penk", "Cartpt", "Bdnf", "Cort", "Calca","Rspo1", "Adcyap1", "Cck", "Ghrh", "Adm",
                       "Syt4", "Syt13",  "Syt7","Syt6", "Syt1", "Syt10",
                       "Gnas", "Scg5", "Scg3", "Scg2",
                       
                       "Cpne8","Cpne6","Cpne5", "Cpne7", "Cpne2",  "Cpne4", "Cpne3",
                       "Camk2n1",	"Camk2b",  "Camk2a","Camk2d",	"Camk2g",		"Camk2n2",
                       "Cacna2d1", "Cacna2d3", "Lin7b","Lin7a","Lin7c", "Prkca", "Prkcg", "Prkcd", "Syndig1l", "Cacng5"),
          #group.by = "segments.simple",
          group.colors = c("#377EB8", "#377EB8", "#377EB8", "#E41A1C", "#E41A1C", "#E41A1C"),
          group.bar.height = 0.02,
          size = 6,
          #angle = 0,
          #hjust = 0.5,
          #vjust = 0.5
) + 
  scale_fill_gradientn(colors=c("white", "grey90", "gray20"))


######################################################## violin plots ######################################################## 
### Subset out precursors
se.data.sct.subset <- subset(se.data.sct, segments.simple %in% c("CeL_appetitive", "CeM_appetitive", 
                                                                 "CeC/L_aversive", "CeM_aversive"))
### Subset individual stages
P4 <- subset(se.data.sct.subset, stage == "P4")
P0 <- subset(se.data.sct.subset, stage == "P0")
P10 <- subset(se.data.sct.subset, stage == "P10")
P4.P10 <- subset(se.data.sct.subset, stage%in% c("P4", "P10"))

### Axon guidance -- Supplementary Fig.6D
VlnPlot(P4.P10, assay = "SCT", layer = "data",
        features = c("Sema3a", "Sema6a", "Sema3c", "Sema5a","Sema3e",
                     "Plxna4","Plxna2","Plxnd1", "Plxna3","Nrp2", "Nrp1"), 
        group.by = "segments.simple", stack = T, fill.by = "ident", cols =  c("#377EB8", "#377EB8",  "#E41A1C", "#E41A1C"), flip=TRUE)

### Transcriptional factors -- Supplementary Fig.6E
VlnPlot(P4, assay = "SCT", layer = "data",
        features = c("Pou3f3","Isl1", "Dlx6", "Sox4","Klf7", "Sox2", "Sox5"), 
        group.by = "segments.simple", stack = T, fill.by = "ident", cols =  c("#377EB8", "#377EB8",  "#E41A1C", "#E41A1C"), flip=TRUE)
VlnPlot(P0, assay = "SCT", layer = "data",
        features = c("Pou3f3", "Isl1", "Dlx6", "Sox4","Klf7", "Sox2", "Sox5"), 
        group.by = "segments.simple", stack = T, fill.by = "ident", cols =  c("#377EB8", "#377EB8",  "#E41A1C", "#E41A1C"), flip=TRUE)
VlnPlot(P10, assay = "SCT", layer = "data",
        features = c("Pou3f3","Isl1", "Dlx6","Sox4", "Klf7", "Sox2", "Sox5"), 
        group.by = "segments.simple", stack = T, fill.by = "ident", cols =  c("#377EB8", "#377EB8",  "#E41A1C", "#E41A1C"), flip=TRUE)



######################################################## Regulon plots (Chao) ######################################################## 











