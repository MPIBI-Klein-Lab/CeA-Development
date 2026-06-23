library("Libra")
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(VennDiagram)
library(Nebulosa)
library("URD")
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library("Matrix")
library("pheatmap")



########################################################################### Fig. S1
CeA.all.har.sim <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.1/Amygdala.dev.rds")

DimPlot(CeA.all.har.sim, group.by = "annotation.updated", label = T, label.size = 5, repel = T, cols = colorRampPalette(brewer.pal(8, "Set3"))(23))
DimPlot(CeA.all.har.sim, group.by = "stage", cols = stage.colors, pt.size = 0.5)
DimPlot(CeA.all.har.sim, group.by = "rep", pt.size = 0.5, alpha = 0.5, cols = c("#D8B365", "#5AB4AC", "#35978F"))

CeA.all.har.markers <- FindAllMarkers(object = CeA.all.har.sim, only.pos = TRUE, min.pct = 0.2, return.thresh = 0.001)





########################################################################### Fig.1 trajectory




CeA.dev <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.1/CeA.dev.rds")


########################################################################### Fig.1C Stage(without Velocity), 
###Velocity will be provide as Python script (Alisson)

stage.colors <- brewer.pal(6, "PRGn")

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


########################################################################### Fig. S4 A (Pseudotime)

se.data.sub.2 <- subset(CeA.dev, trajectories %in% c("Appetitve", "Aversive"))


#Idents(se.data.sub.2) <- se.data.sub.2$annotation.final
#se.data.sub.2 <- RenameIdents(se.data.sub.2, "CeM_Dlk1" = "CeM_aversive",
 #                             "CeL/C_Prkcd" = "CeL/C_aversive")
#se.data.sub.2$annotation.final <- factor(se.data.sub.2$annotation.final, levels=c("NP_aversive", "CeM_aversive", "CeL/C_aversive", "NP_appetitve", "CeM_appetitive", "CeL_appetitive"))



VlnPlot(se.data.sub.2, features = "pseudotime", group.by = "stage", cols = stage.colors, alpha = 0.3, pt.size = 0.5) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=14))

### Pseudotime is calculated in from diffusion map embeded in the URD package

###########################################################################  URD trajectory
axial.tree <- readRDS(file="~/Documents/backup_mac20250207/Dev_manuscript/data.submission/urd.tree.k40.s45.final.final.rds")

############################################################## Fig. S4 C
plotDim(axial.tree, "visitfreq.log.1", plot.title="CeL_Prkcd", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.2", plot.title="CeC_Cdh9.Calcrl", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.3", plot.title="CeM_Dlk1", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.11", plot.title="CeL_Sst", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.12", plot.title="CeL_Nts.Tac2", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.7", plot.title="CeM_Il1rapl2.Tafa1", transitions.plot=10000)
plotDim(axial.tree, "visitfreq.log.9", plot.title="CeM_Vdr", transitions.plot=10000)
### Detailed building of URD tree will be provided in another .r file

########################################################################### Fig. 1 F (URD)


plotTree(axial.tree, "stage", title="Tdtomato", discrete.colors = stage.colors,
         tree.alpha = 0.5, cell.alpha = 0.7, cell.size = 2.5, hide.y.ticks = F, tree.size = 2) + 
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        legend.text = element_text(size=24),  # Increase legend text size
        legend.title = element_text(size=24)) # Increase legend title size

########################################################################### Fig. 1 G (URD)





plotTree(axial.tree, "trajectories", title="Tdtomato", discrete.colors = c("#377EB8","#E41A1C"),
         tree.alpha = 0.5, cell.alpha = 0.6, cell.size = 2.5, hide.y.ticks = F, tree.size = 2) + 
  theme(axis.text=element_text(size=32),
        axis.title=element_text(size=32),
        legend.text = element_text(size=32),  # Increase legend text size
        legend.title = element_text(size=32)) # Increase legend title size



########################################################################### Fig. 1 H (URD)

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



DefaultAssay(CeA.dev) <- "RNA"

# redo normalization using sct
se.data.list <- SplitObject(CeA.dev, split.by = "batch")
se.data.list <- lapply(X = se.data.list, FUN = function(x) {
  SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE)
})


axial.tree <- AddScaledExpression(axial.tree, se.data.list, "Rbp1")
axial.tree <- AddScaledExpression(axial.tree, se.data.list, "Tdtomato")
axial.tree <- AddScaledExpression(axial.tree, se.data.list, "Sox5")



Tdtomato.samples <- colnames(axial.tree@count.data)[axial.tree@meta$batch %in% c("E15_1_all", "E15_2_CeA", "E18_1_A", "E18_2_P", "P0_2_Htr", "P21_2_Htr", "P4_2_Htr", "P10_2_Htr")]
#axial.tree.Tdtomato <- urdSubset(axial.tree, cells.keep = Tdtomato.samples)      

#axial.tree.Tdtomato <- AddScaledExpression(axial.tree.Tdtomato, se.data.list, "Tdtomato")
#axial.tree.Tdtomato <- AddScaledExpression(axial.tree.Tdtomato, se.data.list, "Sox5")

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


########################################################################### Fig. S4 pseudobulk method




E15.cea <- subset(scCeA.gaba.E15, idents = c("Aversive", "Appetitve"))
E15.cea$orig.ident <- Idents(E15.cea)
DE <- run_de(E15.cea, 
             cell_type_col = "stage", replicate_col = "rep", label_col = "orig.ident", de_family = "pseudobulk", de_method = "edgeR", de_type = "LRT")
DE <- as.data.frame(DE)



#saveRDS(DE, "~/Documents/backup_mac20250207/Dev_manuscript/DEgeneList/E15.pseudobulk.full.rds")
#DE <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/DEgeneList/E15.pseudobulk.full.rds")
write.csv(DE, file = "~/Documents/backup_mac20250207/Dev_manuscript/data.submission/E15.DE.pseudobulk.csv")

#checking the top p-value genes
DE %>% filter(-log10(p_val_adj) > 1.5) -> top.de
top.de <- top.de %>% arrange(p_val)
saveRDS(top.de, "~/Documents/backup_mac20250207/Dev_manuscript/DEgeneList/E15.pseudobulk.significant.rds")   
DoHeatmap(scCeA.gaba.E15, features = top.de$gene, assay = "RNA") + NoLegend() + scale_fill_gradientn(colors=c("#6395c7", "white", "#e06ead"))   


#E15.cea.markers$diffexpressed[E15.cea.markers$avg_log2FC > 0.25 & FastFed.markers[[i]]$p_val_adj < 0.001] <- "yes"
DE$diffexpressed <- "NO"
DE$diffexpressed[DE$avg_logFC > 1] <- "Highlit_up"
DE$diffexpressed[DE$avg_logFC < -1] <- "Highlit_down"
#FastFed.markers[[i]]$diffexpressed[FastFed.markers[[i]]$gene %in% c("Tet2", "Mgll", "Grin2a", "Arfgap1", "Ush1c", "Add2")] <- "Highlit_special"

# Volcano plot colors
#mycolors <- c("red", "blue", "gray20")
#names(mycolors) <- c("Highlit_down", "Highlit_up", "NO")
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


########################################################################### Fig. 1J Manhatten distance (Chao)






















