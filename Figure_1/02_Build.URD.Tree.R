library("URD")
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library("Matrix")
library("pheatmap")


# set development colors
stage.colors <- brewer.pal(6, "PRGn")


# Load annotated data
se.data <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.1/CeA.dev.rds")

DefaultAssay(se.data) <- "RNA"

# redo normalization using sct
se.data.list <- SplitObject(se.data, split.by = "batch")
se.data.list <- lapply(X = se.data.list, FUN = function(x) {
  SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE)
}) # regress out percent.mito and nFeature_RNA, nCount_RNA is not so helpful for diffusion map calculation

se.data.features <- SelectIntegrationFeatures(object.list = se.data.list, nfeatures = 3000)
se.data.sct <- merge(se.data.list[[1]], se.data.list[2:length(se.data.list)])
se.data.sct$stage <- factor(se.data.sct$stage, levels = c("E15", "E18", "P0", "P4", "P10", "P21"))

# subset, only keep CeA cell populations
se.data.cea <- subset(se.data.sct, trajectories %in% c("Appetitve", "Aversive"))


################################### Build URD object #############################################
### using DATA slot 

count.data <- GetAssayData(object = se.data.cea[["SCT"]], layer = "counts")
#sct.data <- GetAssayData(object = se.data.sub.1[["SCT"]], layer = "data")
sct.scaledata <- GetAssayData(object = se.data.cea[["SCT"]], layer = "scale.data")
sct.scaledata <- as(sct.scaledata, "dgCMatrix")

metadata <- se.data.cea@meta.data

# Creat URD object
object <- createURD(count.data = count.data, meta = metadata, min.cells=0, min.genes=0, min.counts=0, gene.max.cut=50000, max.genes.in.ram=50000)
#object@count.data <- count.data
object@logupx.data <- sct.scaledata
object@var.genes <- se.data.features


# import UMAP from hamonized (harmony) data
object@tsne.y <- as.data.frame(se.data@reductions$umap@cell.embeddings)
colnames(object@tsne.y) <- c("tSNE1", "tSNE2")
plotDim(object, "stage", plot.title = "tSNE: Stage", reduction.use = "tsne", discrete.colors = stage.colors)


# Root cells (pseudotime 0): 
# cells from one E15 dataset collected at a slightly earlier developmental stage (E15_1_all) were assigned as the earliest developmental population 
root.cells <- whichCells(object, "batch", "E15_1_all")
# Tip clusters
P21.cells <- whichCells(object, "stage", "P21")
object@group.ids[P21.cells, "tip.clusters"] <- as.character(object@meta[P21.cells, "P21.annotated"])

# replace P21 cluster name (build tree function requir only cluster numbers as input): 
object@group.ids$tip.numbers[object@group.ids$tip.clusters == "CeL_Prkcd"] <- "1"
object@group.ids$tip.numbers[object@group.ids$tip.clusters == "CeC_Cdh9.Calcrl"] <- "2"
object@group.ids$tip.numbers[object@group.ids$tip.clusters == "CeM_Dlk1"] <- "3"

#object@group.ids$tip.numbers[object@group.ids$tip.clusters == "Ast_Dach2"] <- "4"
#object@group.ids$tip.numbers[object@group.ids$tip.clusters == "Ast_Synpo2"] <- "5"

object@group.ids$tip.numbers[object@group.ids$tip.clusters == "CeM_Drd2.Rai14"] <- "6"
object@group.ids$tip.numbers[object@group.ids$tip.clusters == "CeM_Il1rapl2.Tafa1"] <- "7"
object@group.ids$tip.numbers[object@group.ids$tip.clusters == "CeM_Tac1.Sst"] <- "8"
object@group.ids$tip.numbers[object@group.ids$tip.clusters == "CeM_Vdr"] <- "9"

#object@group.ids$tip.numbers[object@group.ids$tip.clusters == "IPAC_Pde1c"] <- "10"

object@group.ids$tip.numbers[object@group.ids$tip.clusters == "CeL_Sst"] <- "11"
object@group.ids$tip.numbers[object@group.ids$tip.clusters == "CeL_Nts.Tac2"] <- "12"

#object@group.ids$tip.numbers[object@group.ids$tip.clusters == "ITC_Foxp2"] <- "13"

plotDim(object, "tip.numbers", plot.title = "tSNE: Stage", reduction.use = "tsne")



################################### Data Processing #####################################
### Calculate diffusion map
object <- calcDM(object, knn = 40, sigma.use = 45)
plotDimArray(object = object, reduction.use = "dm", dims.to.plot = 1:18, label = "stage", plot.title = "", outer.title = "STAGE - Diffusion Map Sigma 9", legend = F, alpha = 0.3, discrete.colors = stage.colors)
plotDim(object, "stage", transitions.plot = 10000, transitions.alpha = 0.1, plot.title="Developmental stage (with transitions)")

### Calculate pseudotime
flood.result <- floodPseudotime(object, root.cells = root.cells, n = 50, minimum.cells.flooded = 2, verbose = T)
object <- floodPseudotimeProcess(object, flood.result, floods.name="pseudotime")
pseudotimePlotStabilityOverall(object)
plotDim(object, "pseudotime")
plotDists(object, "pseudotime", "stage", plot.title="Pseudotime by stage")
plotDists(object, "pseudotime", "trajectories", plot.title="Pseudotime by stage")

### define direction bias
object.ptlogistic <- pseudotimeDetermineLogistic(object, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)
object.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(object, "pseudotime", logistic.params=object.ptlogistic))

### Perform random walk
object.walks <- simulateRandomWalksFromTips(object, tip.group.id="tip.numbers", root.cells=root.cells, transition.matrix = object.biased.tm, n.per.tip = 25000, root.visits = 1, verbose = T)
object <- processRandomWalksFromTips(object, object.walks, verbose = T)
plotDim(object, "tip.numbers", plot.title="Cells in each tip")
plotDim(object, "visitfreq.log.1", plot.title="CeL_Prkcd", transitions.plot=10000)
plotDim(object, "visitfreq.log.2", plot.title="CeC_Cdh9.Calcrl", transitions.plot=10000)
plotDim(object, "visitfreq.log.3", plot.title="CeM_Dlk1", transitions.plot=10000)
plotDim(object, "visitfreq.log.11", plot.title="CeL_Sst", transitions.plot=10000)
plotDim(object, "visitfreq.log.12", plot.title="CeL_Nts.Tac2", transitions.plot=10000)
plotDim(object, "visitfreq.log.7", plot.title="CeM_Il1rapl2.Tafa1", transitions.plot=10000)
plotDim(object, "visitfreq.log.9", plot.title="CeM_Vdr", transitions.plot=10000)

### Build tree
axial.tree <- loadTipCells(object, "tip.numbers")
axial.tree <- buildTree(axial.tree, pseudotime = "pseudotime", divergence.method = "preference", cells.per.pseudotime.bin = 50, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)

axial.tree <- nameSegments(axial.tree, segments=c("1", "2", "3", "4", "5", "6", 
                                                  "7", "8", "9", "10", "11", "12", "13"), 
                           segment.names = c("CeL_Prkcd", "CeC_Cdh9.Calcrl", "CeM_Dlk1", "Ast_Dach2", "Ast_Synpo2", "CeM_Drd2.Rai14", 
                                             "CeM_Il1rapl2.Tafa1", "CeM_Tac1.Sst", "CeM_Vdr", "IPAC_Pde1c", "CeL_Sst", "CeL_Nts.Tac2", "ITC_Foxp2"))

plotTree(axial.tree, "trajectories", title="Developmental Stage")
plotTree(axial.tree, "stage", title="Developmental Stage")
plotTree(axial.tree, "pseudotime", title="Developmental Stage")

#saveRDS(axial.tree, file="~/Documents/Dev_manuscript/fig.2_precursors/urd.tree.k40.s45.final.final.rds")
#saveRDS(object, file="~/Documents/Dev_manuscript/fig.2_precursors/urd.object.k40.s45.final.final.rds")


plotDim(axial.tree, "segment", plot.title="URD tree segment", point.size = 1.5, alpha=0.5)
plotTree(axial.tree, "segment")


#################################################### Plot gene expression on UDR tree ###########################################################
### If one gene is not presented in all stages, the UDR object cannot plot the gene expression on the tree (it throws an error), for example:

for (i in 1:length(se.data.list)) {
  print(paste0(i, " - ", unique(se.data.list[[i]]$batch)))
  print("Rbp1" %in% rownames(se.data.list[[i]]@assays$SCT@scale.data))
}

#### Function that retrieve gene expression form Seurat to the URD:
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



### add gene expression

axial.tree <- AddScaledExpression(axial.tree, se.data.list, "Rbp1")
axial.tree <- AddScaledExpression(axial.tree, se.data.list, "Tdtomato")
axial.tree <- AddScaledExpression(axial.tree, se.data.list, "Sox5")


### Note: Tdtomato is not present in some of the dataset (data sets that dissected from the wild-type mice)
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

plotTree(axial.tree, "Sox5", title="Sox5", color.limits = c(-7, 7),
         tree.alpha = 1, cell.alpha = 0.6, cell.size = 0.5, continuous.colors = c("#FAEE85", "#FFEDA0", "firebrick4")) + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=14))

plotTree(axial.tree, "expression.Rbp1", title="Rbp1", color.limits = c(-7, 7),
         tree.alpha = 1, cell.alpha = 0.6, cell.size = 0.5, continuous.colors = c("#FAEE85", "#FFEDA0", "firebrick4")) + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=14))


### Plot predefined trajectory on tree
plotTree(axial.tree, "trajectories", title="Tdtomato", discrete.colors = c("#377EB8","#E41A1C"),
         tree.alpha = 0.5, cell.alpha = 0.6, cell.size = 2.5, hide.y.ticks = F, tree.size = 2) + 
  theme(axis.text=element_text(size=32),
        axis.title=element_text(size=32),
        legend.text = element_text(size=32),  # Increase legend text size
        legend.title = element_text(size=32)) # Increase legend title size


### Plot developmental stages on tree
plotTree(axial.tree, "stage", title="Tdtomato", discrete.colors = stage.colors,
         tree.alpha = 0.5, cell.alpha = 0.7, cell.size = 2.5, hide.y.ticks = F, tree.size = 2) + 
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        legend.text = element_text(size=24),  # Increase legend text size
        legend.title = element_text(size=24)) # Increase legend title size

##############################################################################################################

### Annotation of tree segments
se.data$segments <- NA
se.data$segments[axial.tree@tree$cells.in.segment[["1"]]] <- "CeL_Prkcd"
se.data$segments[axial.tree@tree$cells.in.segment[["2"]]] <- "CeC_Cdh9.Calcrl"
se.data$segments[axial.tree@tree$cells.in.segment[["3"]]] <- "CeM_Dlk1"
se.data$segments[axial.tree@tree$cells.in.segment[["13"]]] <- "P_ave_LC" ### precursor for CeL_Prkcd and CeC_Cdh9.Calcrl
se.data$segments[axial.tree@tree$cells.in.segment[["18"]]] <- "P_ave" ### precursor for all aversive

se.data$segments[axial.tree@tree$cells.in.segment[["6"]]] <- "CeM_Drd2.Rai14"
se.data$segments[axial.tree@tree$cells.in.segment[["7"]]] <- "CeM_Il1rapl2.Tafa1"
se.data$segments[axial.tree@tree$cells.in.segment[["8"]]] <- "CeM_Tac1.Sst"
se.data$segments[axial.tree@tree$cells.in.segment[["9"]]] <- "CeM_Vdr"
se.data$segments[axial.tree@tree$cells.in.segment[["16"]]] <- "P_app_M" ### precursor for CeM appetitive

se.data$segments[axial.tree@tree$cells.in.segment[["11"]]] <- "CeL_Sst"
se.data$segments[axial.tree@tree$cells.in.segment[["12"]]] <- "CeL_Nts.Tac2"
se.data$segments[axial.tree@tree$cells.in.segment[["17"]]] <- "P_app_L" ### precursor for CeL appetitive

se.data$segments[axial.tree@tree$cells.in.segment[["19"]]] <- "P_app" ### precursor for all appetitive


########################################## Pseudo time plotting - URD related supplementary figure ####################################################################
### add pseudo time to seurat object
#se.data$pseudotime <- NA
#se.data$pseudotime[rownames(axial.tree@pseudotime)] <- axial.tree@pseudotime$pseudotime
#FeaturePlot(se.data, features = "pseudotime", cols = c("red", "blue"))
#saveRDS(se.data, "/Users/songwei/Documents/Dev_manuscript/fig.2_precursors/Seurat.pseudotime.tree.rds")

### Subset for CeA-only populations
se.data.cea <- subset(se.data, trajectories %in% c("Appetitve", "Aversive"))

VlnPlot(se.data.cea, features = "pseudotime", group.by = "stage", cols = stage.colors, alpha = 0.3, pt.size = 0.5) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=14))


############################################### Plotting Root Cells - URD related supplementary figure #################################
app.root.cells <- c(axial.tree@tree$cells.in.nodes$`19-1`, axial.tree@tree$cells.in.nodes$`19-2`)
ave.root.cells <- c(axial.tree@tree$cells.in.nodes$`18-1`, axial.tree@tree$cells.in.nodes$`18-2`)
plotTree(axial.tree, "segment", title="segments",
         tree.alpha = 0.7, cell.alpha = 0.6, cell.size = 0.5,
         cells.highlight = c(ave.root.cells,app.root.cells),
         cells.highlight.alpha = 0.1,
         cells.highlight.size = 5) + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=14))

se.data$root.cells <- "None"
se.data$root.cells[ave.root.cells] <- "ave.root"
se.data$root.cells[app.root.cells] <- "app.root"
DimPlot(se.data, group.by = "root.cells", cols = c("gray90",  "#E41A1C", "#377EB8"), pt.size = 1.5, alpha = 0.7, order = c("app.root", "ave.root", "None"))







