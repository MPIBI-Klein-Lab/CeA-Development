library(Libra)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(VennDiagram)
library(Nebulosa)
library(URD)
library(Matrix)
library(pheatmap)
library(SingleCellExperiment)
library(here)

stage.colors <- brewer.pal(6, "PRGn")
data_dir <- here("data")


########################################################################### Fig. S1
### read in subpallial amygdala developmental atlas
Amygdala.dev <- readRDS(file.path(data_dir, "Amygdala.dev.rds"))
DimPlot(Amygdala.dev, group.by = "annotation.updated", label = T, label.size = 5, repel = T, cols = colorRampPalette(brewer.pal(8, "Set3"))(23))
DimPlot(Amygdala.dev, group.by = "stage", cols = stage.colors, pt.size = 0.5)
DimPlot(Amygdala.dev, group.by = "rep", pt.size = 0.5, alpha = 0.5, cols = c("#D8B365", "#5AB4AC", "#35978F"))

CeA.all.har.markers <- FindAllMarkers(object = Amygdala.dev, only.pos = TRUE, min.pct = 0.2, return.thresh = 0.001)

########################################################################### Fig.1 trajectory
### read in CeA developmental dataset
CeA.dev <- readRDS(file.path(data_dir, "CeA.dev.rds"))

########################################################################### Fig.1C Stage(without Velocity), 
###Velocity will be provide as Python script (Alisson)

DimPlot(CeA.dev, group.by = "stage", label = F, na.value = "gray90",
        cols = stage.colors,
        alpha = 0.2,
        pt.size = 1,
        repel = T) 

########################################################################### Fig.1D
col.trajectories <- c( "#E41A1C","#377EB8",  "#F4CAE4","#E6F5C9","#d9d9d9", "#B3E2CD", "#FFF2AE", "#FDCDAC")
names(col.trajectories) <- c("Aversive", "Appetitve", "Ast region", "IPAC", "Precursor Dlx1","GP", "BLA", "ITC")
DimPlot(CeA.dev, group.by = "trajectories", label = F, na.value = "gray90",
        cols = col.trajectories,
        alpha = 0.25,
        pt.size = 1,
        repel = T) 

########################################################################### Fig.1E
col.P21 <- c("#bd4146", "maroon3", "violet", 
             "#6adc88", "#7EA2DF", "#66b0d7", "#80ba8a",  
             "#285B90", "cyan4",
             "grey90", "grey90", "grey90", "grey90", "grey90", "grey90")

names(col.P21) <- c("CeC_Cdh9.Calcrl", "CeL_Prkcd", "CeM_Dlk1", 
                    "CeM_Tac1.Sst", "CeM_Drd2.Rai14", "CeM_Il1rapl2.Tafa1", "CeM_Vdr", 
                    "CeL_Sst", "CeL_Nts.Tac2", 
                    "Ast_Synpo2", "Ast_Dach2",  "IPAC_Pde1c", "ITC_Foxp2", "BLA_Gbx1.Lhx8", "BLA_Nkx2-1")

CeA.dev$P21.annotated <- factor(CeA.dev$P21.annotated, levels = c("CeC_Cdh9.Calcrl", "CeL_Prkcd", "CeM_Dlk1", 
                                                                  "CeM_Il1rapl2.Tafa1", "CeM_Vdr", "CeM_Drd2.Rai14", "CeM_Tac1.Sst",
                                                                  "CeL_Sst", "CeL_Nts.Tac2",
                                                                  "Ast_Synpo2", "Ast_Dach2",  "IPAC_Pde1c", "ITC_Foxp2", "BLA_Gbx1.Lhx8", "BLA_Nkx2-1"))
DimPlot(CeA.dev, 
        group.by = "P21.annotated", label = F, 
        cols = col.P21, na.value = "grey90", 
        pt.size = 1, alpha = 0.6, 
        repel = T, order = F,
        label.size = 0)

########################################################################### Fig.S3 A
DimPlot(CeA.dev, group.by = "seurat_clusters", label = T)

########################################################################### Fig.S3 B
FeaturePlot(CeA.dev, features = c("Syt4", "Pnoc", "Tafa1", "Nts", "Sst", 
                                  "Sox5", "Chodl", "Prkcd", "Dlk1", "Calcrl",
                                  "Gucy1a1", "Drd2", "Rarb", "Pde1c"), cols = c("#FFEDA0", "firebrick4"), ncol = 6)



########################################################################### Fig.S3 E-G P21-Adult mapping (Chao) 


########################################################################### Fig. S4 A, also in '02_Build.URD.Tree.R'
se.data.cea <- subset(CeA.dev, trajectories %in% c("Appetitve", "Aversive"))
VlnPlot(se.data.cea, features = "pseudotime", group.by = "stage", cols = stage.colors, alpha = 0.3, pt.size = 0.5) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=14))
### Pseudotime is calculated from diffusion map embedded in the URD package

###########################################################################  URD trajectory, also in '02_Build.URD.Tree.R'
### Read URD tree object
axial.tree <- readRDS(file.path(data_dir, "urd.tree.rds"))

############################################################## Fig. S4 C
plotDim(axial.tree, "visitfreq.log.1", plot.title="CeL_Prkcd", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.2", plot.title="CeC_Cdh9.Calcrl", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.3", plot.title="CeM_Dlk1", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.11", plot.title="CeL_Sst", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.12", plot.title="CeL_Nts.Tac2", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.7", plot.title="CeM_Il1rapl2.Tafa1", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.9", plot.title="CeM_Vdr", transitions.plot=10000)
### Detailed building of URD tree will be provided in another .r file

########################################################################### Fig. 1 F, also in '02_Build.URD.Tree.R'
plotTree(axial.tree, "stage", title="Tdtomato", discrete.colors = stage.colors,
         tree.alpha = 0.5, cell.alpha = 0.7, cell.size = 2.5, hide.y.ticks = F, tree.size = 2) + 
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        legend.text = element_text(size=24),  # Increase legend text size
        legend.title = element_text(size=24)) # Increase legend title size

########################################################################### Fig. 1 G, , also in '02_Build.URD.Tree.R'
plotTree(axial.tree, "trajectories", title="Tdtomato", discrete.colors = c("#377EB8","#E41A1C"),
         tree.alpha = 0.5, cell.alpha = 0.6, cell.size = 2.5, hide.y.ticks = F, tree.size = 2) + 
  theme(axis.text=element_text(size=32),
        axis.title=element_text(size=32),
        legend.text = element_text(size=32),  # Increase legend text size
        legend.title = element_text(size=32)) # Increase legend title size

########################################################################### Fig. 1 H, also in '02_Build.URD.Tree.R'
### Function for mapping gene expression from individual dataset onto the tree
AddScaledExpression <- function(
    tree.object, # URD tree object for adding the metadata
    sct.list, # List of SCT transformed seurat objects
    gene # Character of length 1. genes that are not present in the scale.data of all the objects in the list, thus will not present in the same slot of the merged/integrated object 
){
  
  newcolumn <- paste0("expression.", gene)
  tree.object@meta[, newcolumn] <- 0
  for (i in 1:length(sct.list)) {
    if (gene %in% rownames(sct.list[[i]]@assays$SCT@scale.data)) {
      cells <- intersect(colnames(tree.object@logupx.data), colnames(sct.list[[i]]))
      tree.object@meta[cells, newcolumn] <-sct.list[[i]]@assays$SCT@scale.data[gene, cells]
    }
  }
  return(tree.object)
}

### Normalization using sct for expression visualization
DefaultAssay(CeA.dev) <- "RNA"
se.data.list <- SplitObject(CeA.dev, split.by = "batch")
se.data.list <- lapply(X = se.data.list, FUN = function(x) {
  SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE)
})

axial.tree <- AddScaledExpression(axial.tree, se.data.list, "Rbp1")
axial.tree <- AddScaledExpression(axial.tree, se.data.list, "Tdtomato")
axial.tree <- AddScaledExpression(axial.tree, se.data.list, "Sox5")

### Tdtomato does not exist in all brain samples
Tdtomato.samples <- colnames(axial.tree@count.data)[axial.tree@meta$batch %in% c("E15_1_all", "E15_2_CeA", "E18_1_A", "E18_2_P", "P0_2_Htr", "P21_2_Htr", "P4_2_Htr", "P10_2_Htr")]
plotTree(axial.tree, "expression.Tdtomato", title="Tdtomato", 
         color.limits = c(-5, 7),
         continuous.colors = c("gold", "#FFEDA0", "firebrick4"),
         tree.alpha = 1, cell.alpha = 1, cell.size = 2, hide.y.ticks = F, tree.size = 2,
         cells.highlight = Tdtomato.samples,
         cells.highlight.alpha = 1,
         cells.highlight.size = 2) + 
  theme(axis.text=element_text(size=28),axis.title=element_text(size=28))

plotTree(axial.tree, "Sox5", title="Sox5", 
         #color.limits = c(-7, 10),
         continuous.colors = c("gold", "#FFEDA0", "firebrick4"),
         tree.alpha = 1, cell.alpha = 1, cell.size = 2, hide.y.ticks = F, tree.size = 2,
         #cells.highlight = Tdtomato.samples,
         #cells.highlight.alpha = 1,
         #cells.highlight.size = 1
) + 
  theme(axis.text=element_text(size=28),axis.title=element_text(size=28))

########################################################################### Fig. S4 G Define tree segments
plotTree(axial.tree, "segment",
         tree.alpha = 0.3, cell.alpha = 1, cell.size = 1) + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=14))  

segment.colors.P_ave <- c("gray90", "gray90", "gray90", "gray90", "#E41A1C", 
                          "gray90", "gray90", "gray90", "gray90", "gray90",
                          "gray90", "gray90", "gray90", "gray90")
names(segment.colors.P_ave) <- c("1", "2", "3", "13", "18",
                                 "6", "7", "8", "9", "16",
                                 "11", "12", "17", "19")

segment.colors.CeLC_ave <- c("#E41A1C", "#E41A1C", "gray90", "#E41A1C", "gray90", 
                             "gray90", "gray90", "gray90", "gray90", "gray90",
                             "gray90", "gray90", "gray90", "gray90")
names(segment.colors.CeLC_ave) <- c("1", "2", "3", "13", "18",
                                    "6", "7", "8", "9", "16",
                                    "11", "12", "17", "19")

segment.colors.CeM_ave <- c("gray90", "gray90", "#E41A1C", "gray90", "gray90", 
                            "gray90", "gray90", "gray90", "gray90", "gray90",
                            "gray90", "gray90", "gray90", "gray90")
names(segment.colors.CeM_ave) <- c("1", "2", "3", "13", "18",
                                   "6", "7", "8", "9", "16",
                                   "11", "12", "17", "19")


segment.colors.P_app <- c("gray90", "gray90", "gray90", "gray90", "gray90", 
                          "gray90", "gray90", "gray90", "gray90", "gray90",
                          "gray90", "gray90", "gray90", "#377EB8")
names(segment.colors.P_app) <- c("1", "2", "3", "13", "18",
                                 "6", "7", "8", "9", "16",
                                 "11", "12", "17", "19")

segment.colors.CeM_app <- c("gray90", "gray90", "gray90", "gray90", "gray90", 
                            "#377EB8", "#377EB8", "#377EB8", "#377EB8", "#377EB8",
                            "gray90", "gray90", "gray90", "gray90")
names(segment.colors.CeM_app) <- c("1", "2", "3", "13", "18",
                                   "6", "7", "8", "9", "16",
                                   "11", "12", "17", "19")

segment.colors.CeL_app <- c("gray90", "gray90", "gray90", "gray90", "gray90", 
                            "gray90", "gray90", "gray90", "gray90", "gray90",
                            "#377EB8", "#377EB8", "#377EB8", "gray90")
names(segment.colors.CeL_app) <- c("1", "2", "3", "13", "18",
                                   "6", "7", "8", "9", "16",
                                   "11", "12", "17", "19")

plotTree(axial.tree, "segment", discrete.colors = segment.colors.P_ave,
         tree.alpha = 0.3, cell.alpha = 1, cell.size = 3) + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=14))

########################################################################### Fig. S4 pseudobulk method for deferential expression
E15 <- subset(CeA.dev, stage == "E15")
E15.cea <- subset(E15, subset = trajectories %in% c("Aversive", "Appetitve"))
E15.cea$orig.ident <- E15$trajectories

### Note: The Libra package was built for older versions of Seurat (v3/v4) and relies on GetAssayData(..., slot = "counts")
### Bypass Seurat v5 Using SingleCellExperiment
sce_obj <- as.SingleCellExperiment(E15.cea, assay = "RNA")

### pseudobulk DE run
DE <- run_de(sce_obj, 
             cell_type_col = "stage", 
             replicate_col = "rep", 
             label_col = "orig.ident", 
             de_family = "pseudobulk", 
             de_method = "edgeR", 
             de_type = "LRT")

DE <- as.data.frame(DE)
#write.csv(DE, file = "~/Documents/backup_mac20250207/Dev_manuscript/data.submission/E15.DE.pseudobulk.csv")

#checking the top p-value genes
#DE %>% filter(-log10(p_val_adj) > 1.5) -> top.de
#top.de <- top.de %>% arrange(p_val)
#saveRDS(top.de, "~/Documents/backup_mac20250207/Dev_manuscript/DEgeneList/E15.pseudobulk.significant.rds")   


### Volcano plots
DE$diffexpressed <- "NO"
DE$diffexpressed[DE$avg_logFC > 1] <- "Highlit_up"
DE$diffexpressed[DE$avg_logFC < -1] <- "Highlit_down"

# Volcano plot colors
mycolors <- c("#E41A1C", "#377EB8", "gray90")
names(mycolors) <- c("Highlit_down", "Highlit_up", "NO")
DE$diffexpressed <- as.factor(DE$diffexpressed)

# List of genes to be labeled
DE$gene_label <- ifelse(DE$gene %in% c("Pax6", "Klf5", "Igfbp5", "Sox5", "Tdtomato", 
                                       "Pcdh7", "Htr2a", "Sox2", "Rbp1", "Cyp26b1", "Adora2a", "Rbp1", "S100a10", "Trpc5",
                                       "Cck", "Ebf1", "Syt4", "Syt6", "Nts", "Calcrl", "Ecel1", "Sema6a", "Rspo3"),
                        paste0("italic('", DE$gene, "')"), NA)

### Plot with threshold genes
ggplot(data=DE, aes(x=avg_logFC, y=-log10(p_val), col=diffexpressed)) + 
  geom_point(alpha = 0.8, size = 3) + 
  theme_minimal() +
  geom_text_repel(data = DE[!is.na(DE$gene_label), ],  # Use only selected genes
                  aes(label = gene_label),
                  size = 10,
                  color = "gray20",
                  max.overlaps = 10, force = 41,
                  parse = TRUE) +  # Enable parsing for italics
  #ylim(0, 20) +
  scale_colour_manual(values = mycolors) +
  theme_bw(30) +  # Adjust base font size
  theme(panel.grid.minor = element_blank()) 

########################################################################### Fig. 1J Manhattan distance (Chao)













