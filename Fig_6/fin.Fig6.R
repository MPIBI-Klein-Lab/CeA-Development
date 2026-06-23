library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(clue)
library('glmGamPoi')
library("sctransform")
library(MetaNeighbor)
library(SingleCellExperiment)
library(tidyverse)

################################################################## Mapping of human-mouse adult cell types -- Supplementary Fig.12D
### loading data
cca.se.fin <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.6/cleaned.cca.Human.Mouse.Integrated.rds")

cca.se.fin$integrated.clusters <- cca.se.fin$cca.resolution0.8.number
Idents(cca.se.fin) <- cca.se.fin$integrated.clusters

DimPlot(cca.se.fin,  label = F, pt.size = 0.5, group.by = "species", alpha = 0.5, cols = c("#D8B365", "#35978F"), label.size = 20)
DimPlot(cca.se.fin,  label = F, pt.size = 0.5, group.by = "integrated.clusters", label.size = 20)

################################################################## Mapping of human-mouse adult cell types -- Fig.6A
### integrated.clusters: human-mouse shared clusters
### human.clean.cluster.final: human cluster numbers
### adult.annotation.new: adult mouse CeA annotations

human.shared <- as.matrix(table(cca.se.fin$human.clean.cluster.final, cca.se.fin$integrated.clusters))
human.shared.freq <- sweep(human.shared, 1, rowSums(human.shared), FUN = "/")

mouse.shared <- as.matrix(table(cca.se.fin$adult.annotation.new, cca.se.fin$integrated.clusters))
mouse.shared.freq <- sweep(mouse.shared, 1, rowSums(mouse.shared), FUN = "/")


human.mouse <- human.shared %*% t(mouse.shared)
human.mouse.freq <- human.shared.freq %*% t(mouse.shared.freq)


######## ordering
# Step 1: Get the best matching mouse cluster for each human cluster
best_mouse <- apply(human.mouse.freq, 1, function(x) colnames(human.mouse.freq)[which.max(x)])

# Step 2: Build an ordering by sorting human clusters by their top matching mouse cluster index
best_mouse_indices <- match(best_mouse, colnames(human.mouse.freq))
human_order <- order(best_mouse_indices)
mouse_order <- order(match(colnames(human.mouse.freq), best_mouse[human_order]))

# Step 3: Reorder matrix
reordered_matrix <- human.mouse.freq[human_order, mouse_order]

rereordered_matrix <- reordered_matrix[c("7", "2", "22", "19", "20", "9", "12", "8", "11", "24", "18", "3", "15", "17", "14", "21", "6", "10", "13", "23", "0", "1", "5", "4", "16"),
                                       c("Ast_Dach2", "Ast_Synpo2", "BLA_Gbx1.Lhx8",  "BLA_Nkx2-1", "CeL_Prkcd", "CeC_Cdh9.Calcrl", "CeM_Dlk1", 
                                         "CeM_Il1rapl2.Tafa1", "CeM_Drd2.Rai14", "CeM_Tac1.Sst", "CeM_Vdr", "CeL_Nts.Tac2", "CeL_Sst", "IPAC_Pde1c", "ITC_Foxp2")]

# Step 4: Plot
pheatmap(rereordered_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "black"))(100),
         cellwidth = 22,
         cellheight = 22,
         fontsize_col=15,
         fontsize_row=15,
         angle_col = "45"
)

#write.csv(rereordered_matrix, file = "human_mouse_freq.csv", quote = FALSE)

################################################################## Calculation of empirical p-value -- Fig.6A
# Store original similarity matrix
original_matrix <- rereordered_matrix

n_iter <- 1000
perm_max_values <- matrix(0, nrow = nrow(original_matrix), ncol = ncol(original_matrix))
rownames(perm_max_values) <- rownames(original_matrix)
colnames(perm_max_values) <- colnames(original_matrix)

pb <- txtProgressBar(0, n_iter, style = 3)

for (i in 1:n_iter) {
  # Shuffle only human cluster labels
  shuffled_human <- sample(cca.se.fin$human.clean.cluster.final)
  
  # Keep mouse fixed
  human.shared.perm <- as.matrix(table(shuffled_human, cca.se.fin$integrated.clusters))
  human.shared.freq.perm <- sweep(human.shared.perm, 1, rowSums(human.shared.perm), FUN = "/")
  
  mouse.shared <- as.matrix(table(cca.se.fin$adult.annotation.new, cca.se.fin$integrated.clusters))
  mouse.shared.freq <- sweep(mouse.shared, 1, rowSums(mouse.shared), FUN = "/")
  
  # Calculate permuted human-mouse similarity
  perm_matrix <- human.shared.freq.perm %*% t(mouse.shared.freq)
  
  # Track how often permuted similarity exceeds observed
  perm_max_values <- perm_max_values + (perm_matrix >= original_matrix)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Empirical p-values
empirical_pvals <- perm_max_values / n_iter

### Visualization

# Step 1: Floor p-values to avoid log(0)
empirical_pvals[empirical_pvals == 0] <- 1e-6

# Step 2: Transform to -log10 scale
log_pvals <- -log10(empirical_pvals)

# Step 3: Mask non-significant values (e.g., p > 0.01)
log_pvals_masked <- log_pvals
log_pvals_masked[empirical_pvals > 0.01] <- NA

pheatmap(log_pvals,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "black"))(100),
         cellwidth = 22,
         cellheight = 22,
         fontsize_col=15,
         fontsize_row=15,
         angle_col = "45"
)

pheatmap(log_pvals_masked,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "black"))(100),
         cellwidth = 22,
         cellheight = 22,
         fontsize_col = 15,
         fontsize_row = 15,
         angle_col = "45",
         na_col = "yellow")  # Choose a color for NA (e.g., grey or transparent)


### find best matches using thresholding
### Apply threshold: keep only values > 0.15, mask others

threshold.matrix <- human.mouse.freq
threshold.matrix[threshold.matrix <= 0.15] <- NA

threshold.matrix <- threshold.matrix[c("7", "2", "22", "19", "20", "9", "12", "8", "11", "24", "18", "3", "15", "17", "14", "21", "6", "10", "13", "23", "0", "1", "5", "4", "16"),
                                     c("Ast_Dach2", "Ast_Synpo2", "BLA_Gbx1.Lhx8",  "BLA_Nkx2-1", "CeL_Prkcd", "CeC_Cdh9.Calcrl", "CeM_Dlk1", 
                                       "CeM_Il1rapl2.Tafa1", "CeM_Drd2.Rai14", "CeM_Tac1.Sst", "CeM_Vdr", "CeL_Nts.Tac2", "CeL_Sst", "IPAC_Pde1c", "ITC_Foxp2")]

pheatmap(threshold.matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "black"))(100),
         cellwidth = 22,
         cellheight = 22,
         fontsize_col = 15,
         fontsize_row = 15,
         angle_col = "45",
         na_col = "yellow")  # NA values appear as grey (or any neutral color you prefer)

################################################################## Human CeA consensus types -- Fig.6B-C, Suplementary Fig.12 A-C
### loading human annotated dataset
human.se <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.6/clean.human.se.consensus.rds")

DimPlot(human.se, label = T, pt.size = 0.5, label.size = 10, repel = T, group.by = "human.clean.cluster.final")
DimPlot(human.se, group.by = "human.mouse.consensus", alpha = 0.7, pt.size = 0.5, label = F, label.size = 7, repel = T,
        cols = c("#F4CAE4", "#B3E2CD","#FFF2AE", "#E41A1C", "#377EB8","#E6F5C9","#FDCDAC","gray90","gray90","gray90","gray90","gray90","gray90","gray90","gray90"))

FeaturePlot(human.se, features = "SOX5", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            max.cutoff = NA,
            pt.size = 1)
FeaturePlot(human.se, features = "HTR2A", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            max.cutoff = 2,
            pt.size = 1)
FeaturePlot(human.se, features = "PRKCD", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            max.cutoff = 1,
            pt.size = 1)
FeaturePlot(human.se, features = "PNOC", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            max.cutoff = 1,
            pt.size = 1
)


### Appetitive
FeaturePlot(human.se, features = "CACNA1G", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            max.cutoff = 2,
            pt.size = 1
)
FeaturePlot(human.se, features = "TSHZ2", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            #max.cutoff = 1.5,
            pt.size = 1
)

### Aversive
FeaturePlot(human.se, features = "CARTPT", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            #max.cutoff = 1.5,
            pt.size = 1
)
FeaturePlot(human.se, features = "CALCRL", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            #max.cutoff = 1.5,
            pt.size = 1
)

### ITC
FeaturePlot(human.se, features = "FOXP2", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            max.cutoff = 5,
            pt.size = 1
)
FeaturePlot(human.se, features = "TSHZ1", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            #max.cutoff = 1.5,
            pt.size = 1
)

### AST
FeaturePlot(human.se, features = "DRD2", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            #max.cutoff = 5,
            pt.size = 1
)
FeaturePlot(human.se, features = "RARB", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            #max.cutoff = 1.5,
            pt.size = 1
)

### 
FeaturePlot(human.se, features = "NTS", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            max.cutoff = 2.5,
            pt.size = 1
)
FeaturePlot(human.se, features = "PDE1C", cols = c("#FAEE85", "firebrick4"), 
            #min.cutoff = 'q1', 
            #max.cutoff = 1.5,
            pt.size = 1
)

#### DE expression heatmap
human.se$human.clean.cluster.final <- factor(human.se$human.clean.cluster.final, levels = c("7", "2", "22", "19", "20", "9", "12", "8", "11", "24", "18", "3", "15", "17", "14", "21", "6", "10", "13", "23", "0", "1", "5", "4", "16"))
Idents(human.se) <- human.se$human.clean.cluster.final
DefaultAssay(human.se) <- "RNA"
h.markers.rna <- FindAllMarkers(human.se, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)

h.top15.rna <- h.markers.rna %>%
  mutate(cluster = factor(cluster, levels = levels(human.se$human.clean.cluster.final))) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)  # better than top_n(), which is deprecated

DoHeatmap(human.se, features = h.top15.rna$gene, size = 7, raster = T) + 
  scale_fill_gradientn(colors=c("white", "gray95", "gray30"))

################################################################## Metaneighbor analysis of replicable cell types between human and mouse developmental trajectories -- Fig.6D, Suplementary Fig. 12E
### Load SCT normalized human - mouse trajectory merged dataset  (gene name converted to human):
human.mouseDev.sct <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.6/clean.human.mouseDev.merged.rds")

### Build summarizedExperiment object for MetaNeighbour
sample.id <- human.mouseDev.sct$stage
trajectory.id <- human.mouseDev.sct$trajectories
assay.sct <- GetAssayData(human.mouseDev.sct, assay = "SCT", layer = "data")    # SCT makes so much more sense than RNA assay
hm.MN <- SummarizedExperiment(assays = assay.sct)
### Find variable genes for MetaNeighbor
var_genes.1 <- variableGenes(dat = hm.MN, exp_labels = human.mouseDev.sct$stage)


hm.AUROC <- MetaNeighborUS(
  var_genes = var_genes.1,
  dat = hm.MN,
  study_id = sample.id,      # NOT quoted!
  cell_type = trajectory.id,     # NOT quoted!
  symmetric_output = TRUE,
  one_vs_best = FALSE,
  fast_version = TRUE # Fast version is better
)

### heatmap for Suplementary Fig. 12E
pheatmap(hm.AUROC,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "gray90", "gray10"))(100),
         cellwidth = 16,
         cellheight = 16,
         cutree_rows = 4,
         cutree_cols = 4,
         fontsize_col=15,
         fontsize_row=15,
         angle_col = "90"
)


### Barplot for  Fig. 6D

human.scores <- hm.AUROC[c("human|Appetitve", "human|Aversive"),]
# Convert matrix to data frame and move rownames into a column
human.scores.df <- human.scores %>%
  as.data.frame() %>%
  rownames_to_column(var = "From") %>%
  pivot_longer(
    cols = -From,
    names_to = "To",
    values_to = "Score"
  )

human.scores.df <- human.scores.df %>%
  separate(col = To, into = c("Stage", "Region"), sep = "\\|")

# Fix typo in Region and From
human.scores.df <- human.scores.df %>%
  mutate(
    Region = str_replace(Region, "Appetitve", "Appetitive"),
    From   = str_replace(From, "Appetitve", "Appetitive")
  )

human.scores.df <- human.scores.df[!human.scores.df$Stage == "human", ]

# Redefine Region and Stage order
Region_order <- c("Appetitive", "Aversive", "IPAC", "Ast region", "ITC")
stage_order <- c("E15", "E18", "P0", "P4", "P10", "P21")

ordered_RegionStage <- as.vector(outer(stage_order, Region_order, paste, sep = "."))

human.scores.df <- human.scores.df %>%
  mutate(
    Region = factor(Region, levels = Region_order),
    Stage = factor(Stage, levels = stage_order),
    RegionStage = paste(Stage, Region, sep = "."),
    RegionStage = factor(RegionStage, levels = ordered_RegionStage)
  )


cols <- c("#377EB8", "#E41A1C", "#E6F5C9", "#F4CAE4", "#FDCDAC")  # matches Region_order

ggplot(human.scores.df, aes(x = RegionStage, y = Score, fill = Region)) +
  geom_col() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  facet_wrap(~ From, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = cols) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Similarity Scores by Region and Stage",
    x = "Region and Developmental Stage",
    y = "AUROC",
    fill = "Region"
  ) +
  theme_minimal(22) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


################################################################## Conserved expression of gene families between human and mouse -- Fig.6G

hmP4.z.results <- read.csv("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.6/human.mouse.conserved.families.csv")


hmP4.z.results <- hmP4.z.results %>%
  mutate(gene_set = recode(gene_set, !!!name_map))

hmP4.z.results <- hmP4.z.results %>%
  filter(gene_set != "Synaptotagmin")

### Visualization
### genes to label: p<0.001 for both aversive and appetitive
interested.labels <- union(hmP4.z.results$gene_set[hmP4.z.results$appetitive_p < 0.001], 
                           hmP4.z.results$gene_set[hmP4.z.results$aversive_p < 0.001])


ggplot(hmP4.z.results, aes(x = appetitive_delta, y = aversive_delta)) +
  # Background points (unlabeled)
  geom_point(color = "gray20", alpha = 0.35, size = 2) +
  
  # Highlighted points (labeled)
  geom_point(
    data = subset(hmP4.z.results, gene_set %in% interested.labels),
    color = "black",
    size = 2.5
  ) +
  
  # Labels with ggrepel
  geom_text_repel(
    data = subset(hmP4.z.results, gene_set %in% interested.labels),
    aes(label = gene_set),
    size = 5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 1,
    segment.color = "grey50"
  ) +
  
  # Threshold lines - thicker
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "#377EB8", size = 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "#E41A1C", size = 1) +
  
  # Axis labels and title
  labs(
    x = "Appetitive   ΔAUROC",
    y = "Aversive   ΔAUROC",
    title = "ΔAUROC Comparison (Selective Labeling with ggrepel)"
  ) +
  
  # Theme tweaks
  theme_minimal(base_size = 18) +  # increase overall base size
  theme(
    axis.title = element_text(size = 14),
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1.1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )





