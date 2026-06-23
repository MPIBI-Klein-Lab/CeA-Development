library(Seurat)
library(MetaNeighbor)
library(SummarizedExperiment)
library(RColorBrewer)
library(svMisc)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Nebulosa)
library(rlang)
source("~/Documents/backup_mac20250207/Dev_manuscript/fig.4/Final2025.MN.TimeSeries.source.R")

MN_data <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/Seurat.pseudotime.tree.0301.rds")

### SCT Normalization
DefaultAssay(MN_data) <- "RNA"
MN_data.list <- SplitObject(MN_data, split.by = "batch")
MN_data.list <- lapply(X = MN_data.list, FUN = function(x) {
  SCTransform(x, vst.flavor = "v2", vars.to.regress = c("percent.mito", "nFeature_RNA", "nCount_RNA"))
})
# regress our percent.mito and nFeature_RNA is important, regress out nCount_RNA also helps

MN_data.sct <- merge(MN_data.list[[1]], MN_data.list[2:length(MN_data.list)])


################################## Run metaneighbour for individual gene set -- Fig.3B and Supplementary Fig.6C ################################## 

############# Define custome gene sets
Semaphorins <- grep(pattern = "^Sema", x = rownames(MN_data.sct), value = TRUE)
Semaphorin.receptors <- c("Nrp2", "Nrp1", "Plxna2", "Plxna1", "Plxnd1", "Plxnb1", "Plxna4", "Plxna3", "Plxnb2", "Plxnc1", "Plxnc2", "Plxnb3", "Plxdc2", "Plxdc1")
Ephorins <- grep(pattern = "^Efn", x = rownames(MN_data.sct), value = TRUE)
RGS <- grep(pattern = "^Rgs", x = rownames(MN_data.sct), value = TRUE)
neuropeptides <- c('Adcyap1', 'Adm', 'Agrp', 'Avp', 'Bdnf', 'Bmp2', 'Bmp4', 'Calca', 'Cartpt', 'Cck', 'Cort', 'Crh', 'Edn1', 'Edn3', 'Fgf9', 'Gal', 'Ghrh', 'Grp', 'Hcrt', 'Igf1', 'Inhba', 'Inhbb', 'Kiss1', 'Nmb', 'Nms', 'Nmu', 'Nppc', 'Npvf', 'Npw', 'Npy', 'Nrtn', 'Nts', 'Oxt', 'Pdyn', 'Penk', 'Pnoc', 'Pomc', 'Prok2', 'Pthlh', 'Ptn', 'Qrfp', 'Rln1', 'Rspo1', 'Rspo2', 'Rxfp1', 'Sst', 'Tac1', 'Tac2', 'Tgfb3', 'Trh', 'Tshb', 'Vip', 'Wnt2', 'Wnt4', 'Wnt5a')
Copines <-	c("Cpne1",	"Cpne2",	"Cpne3",	"Cpne4",	"Cpne5",	"Cpne6",	"Cpne7",	"Cpne8",	"Cpne9")
Granins <-	c("Chga",	"Chgb",	"Gnas",	"Pcsk1n",	"Scg2",	"Scg3",	"Scg5",	"Vgf")

############# running test by calling 'GeneSetTest' function
GeneSet.list <- list("Semaphorins"=Semaphorins, 
                     "Semaphorin.receptors"=Semaphorin.receptors, 
                     "Ephorins"=Ephorins, 
                     "RGS"=RGS, 
                     "neuropeptides"=neuropeptides, 
                     "Copines"=Copines, 
                     "Granins"=Granins)

result.list <- list()
for (i in 1:length(GeneSet.list)) {
  meta.result <- GeneSetTest(genes = GeneSet.list[[i]], seu.obj = MN_data.sct,seu.assay = "SCT", seu.layer = "data",n.run = 15, split.data = c(0.8, 0.2),chance.use = "ALL")
  result.list[[i]] <- meta.result
  names(result.list)[i] <- names(GeneSet.list)[i]
}



### Calculate gene family expression

for (i in 1:length(GeneSet.list)) {
  MN_data.sct[[paste0(names(result.list)[[i]], ".perc")]] <- PercentageFeatureSet(MN_data.sct, 
                                                                                  features = intersect(rownames(MN_data.sct), GeneSet.list[[i]]), 
                                                                                  assay = "SCT")
}

MN_data.sct.subset <- subset(MN_data.sct, trajectories %in% c("Aversive", "Appetitve"))
df.expression <- MN_data.sct.subset[[c("trajectories", "stage",
                                       "Semaphorins.perc", 
                                       "Semaphorin.receptors.perc", 
                                       "Ephorins.perc", 
                                       #"Ephrin.receptors.perc", 
                                       #"Roundabout.perc",
                                       "RGS.perc", 
                                       #"AMPA.subunits.perc", 
                                       "neuropeptides.perc", 
                                       #"neuropeptides.receptors.perc", 
                                       #"cytoskeleton.binding.perc", 
                                       #"CamKIIs.perc", 
                                       #"Synaptotagmins.perc", 
                                       "Copines.perc", 
                                       "Granins.perc")]]



#### Plots for visualization
for (i in seq_along(GeneSet.list)) {
  # p1
  p1 <- ggplot(result.list[[i]], aes(x = stage, y = scores, group = chance, fill = cell_type)) +
    geom_point(aes(col = cell_type), position = position_jitter(width = 0.2), alpha = 0.5) +
    facet_wrap(~cell_type, scales = "fixed", ncol = 2) +
    geom_smooth(aes(linetype = chance, fill = cell_type), color = "black") +
    scale_fill_manual(values = c(Appetitve = "#377EB8", Aversive = "#E41A1C")) +
    scale_color_manual(values = c(Appetitve = "#377EB8", Aversive = "#E41A1C")) +
    theme_bw(27) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_blank(),
          legend.title = element_blank(),
          legend.position = c(-2, -2),
          legend.key = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 20)) +
    labs(title = "AUROC", x = NULL, y = NULL) +
    ylim(0.48, 0.93)
  
  # p2
  gene_name <- names(result.list)[[i]]
  perc_col_name <- paste0(gene_name, ".perc")
  upper_limit <- quantile(df.expression[[perc_col_name]], 0.99, na.rm = TRUE)
  
  p2 <- ggplot(df.expression, aes(x = stage, y = !!sym(perc_col_name), fill = trajectories)) +
    geom_boxplot(outlier.colour = "grey80", outlier.shape = 16, outlier.size = 2, notch = TRUE, width = 0.6) +
    geom_smooth(aes(col = trajectories, group = trajectories), method = "loess", linetype = "dashed", se = FALSE, size = 0.8) +
    scale_fill_manual(values = c(Appetitve = "#377EB8", Aversive = "#E41A1C")) +
    scale_color_manual(values = c(Appetitve = "#377EB8", Aversive = "#E41A1C")) +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0, upper_limit)) +
    theme_bw(27) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 20)) +
    labs(title = "Percent expressed", x = NULL, y = NULL)
  
  # Combine and save
  combined_plot <- plot_grid(p1, p2, rel_widths = c(7, 6))
  
  ggsave(
    filename = file.path("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.3/Single.GeneSets",
                         paste0(names(GeneSet.list)[i], ".pdf")),
    plot = combined_plot,
    width = 14.86,
    height = 4.86,
    units = "in"
  )
}


################################## Effect of gene set size -- Supplementary Fig.6A ################################## 
############# running test by calling 'GeneNumberTest' function

ngene.effect.zoomout <- GeneNumberTest(x=c(5, 10, 20, 30, 40, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300), MN_data, n.run=10)
ngene.effect.zoomin <- GeneNumberTest(x=c(5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35), MN_data, n.run=10)


ngene.effect.zoomin <- ngene.effect.zoomin %>% filter(stage %in% c("E15", "P0", "P21"))
ngene.effect.zoomout <- ngene.effect.zoomout %>% filter(stage %in% c("E15", "P0", "P21"))

# Define colors for cell types, be cautious about the typo of APPETITVE instead of APPETITIVE ?!!!
celltype.colors <- c(Aversive = "#E41A1C", Appetitve = "#377EB8")

p2 <- ggplot(ngene.effect.zoomin, aes(x=ngene, y=scores, group=cell_type, fill=cell_type)) + 
  geom_point(position = "jitter", alpha=0.3, size =1) +
  geom_smooth(aes(col=cell_type), linewidth=1) +
  scale_color_manual(values = celltype.colors) +
  scale_fill_manual(values = celltype.colors) +
  facet_wrap(~stage, scales = "fixed", ncol = 3) +
  theme_bw(24) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_blank(),
        legend.title=element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  labs(x = NULL, y = NULL) +
  ylim(0.45, 0.8)

### Residual plots
# Fit a model for each subgroup and calculate residuals
ngene.effect.zoomin <- ngene.effect.zoomin %>%
  group_by(cell_type, stage) %>%
  do({
    model <- lm(scores ~ ngene, data = .)
    data_frame(ngene = .$ngene, scores = .$scores, residuals = residuals(model))
  }) %>%
  ungroup()

ggplot(ngene.effect.zoomin, aes(x = ngene, y = residuals, color = cell_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # Adds a horizontal line at zero
  geom_point(alpha = 0.6) +  # Add points
  facet_wrap(~stage) +  # Facet by stage
  geom_smooth(method = "lm", se = FALSE) +  # Adds a smooth line, no confidence interval
  labs(title = "Residuals vs. ngene", x = "ngene", y = "Residuals") +
  scale_color_manual(values = c("Appetitve" = "#377EB8", "Aversive" = "#f74747")) +  # Color settings
  theme_minimal() +  # Minimal theme for clarity
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size = 12),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )


p1 <- ggplot(ngene.effect.zoomout, aes(x=ngene, y=scores, group=cell_type, fill=cell_type)) + 
  geom_point(position = position_jitter(width = 6), alpha=0.3) + 
  geom_vline(xintercept = 35, linetype = "dashed", color = "black", linewidth = 1) +
  facet_wrap(~stage, scales = "fixed",  ncol = 3) + 
  theme_bw(24) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(size = 5,
                                    face = "italic", angle = -90),
        legend.title=element_blank(),
        legend.position = "right",
        legend.key = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  labs(x = NULL, y = NULL) +
  #ylim(0.48, 0.9) +
  geom_smooth(aes(col=cell_type)) +
  scale_color_manual(values = celltype.colors) +
  scale_fill_manual(values = celltype.colors) 


plot_grid(p1, p2, ncol = 1, rel_heights = c(4, 3))


################################## Screening for the HGNC gene family list -- Fig.3C ################################## 
############# getting gene family list
library(readxl)
data <- read_excel("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.3/HGNC_genes.xlsx")
gene_list <- lapply(data, as.character)
names(gene_list) <- make.names(names(gene_list))
gene_list <- lapply(gene_list, function(x) x[!is.na(x)])

hist(sapply(gene_list, length), breaks = 35, xlab = '#Genes', cex.axis = 2)

gene_list.50 <- gene_list[sapply(gene_list, function(x) length(x) <= 50, simplify = TRUE)]

############# running test by calling 'GeneNumberTest' function
####### GeneNumberTest test a list of gene sets, return the delta-AUROC scores and output a matrix for heat map

DefaultAssay(MN_data.sct) <- "SCT"
HGNCsets.results <- GeneSetTestMatrix(gene_list.50, MN_data.sct, 
                                          seu.assay = "SCT",           
                                          seu.layer = "data",
                                          n.run = 7,
                                          split.data = c(0.8, 0.2),
                                          chance.use = "ALL" )


App.HGNC.results <- na.omit(HGNCsets.results$Appetitive)
Ave.HGNC.results <- na.omit(HGNCsets.results$Aversive)

combined.HGNC.results <- cbind(App.HGNC.results, Ave.HGNC.results)


################################## Screening from the manually curated gene set -- Fig.3C ################################## 
############# getting gene family list
library(readxl)
custom.data <- read_excel("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.3/Final2025.Families.axon.synapse.Rinput.xlsx")
custom.data <- lapply(custom.data, as.character)
names(custom.data) <- make.names(names(custom.data))
custom.data <- lapply(custom.data, function(x) x[!is.na(x)])

hist(sapply(custom.data, length), breaks = 10, xlab = '#Genes', cex.axis = 2)

############# running test by calling 'GeneNumberTest2' function
####### GeneNumberTest2 test a list of gene sets, return both the delta-AUROC scores and p-value comparing the paired gene set and randomly selected genes

DefaultAssay(MN_data.sct) <- "SCT"
custom.geneset.result <- GeneSetTestMatrix2(custom.data, MN_data.sct, 
                                                seu.assay = "SCT",           
                                                seu.layer = "data",
                                                n.run = 7,
                                                split.data = c(0.8, 0.2),
                                                chance.use = "ALL" )

combined.custom.results <- cbind(custom.geneset.result$Appetitive.mean, custom.geneset.result$Aversive.mean)
combined.custom.results.p <- cbind(custom.geneset.result$Appetitive.pValue, custom.geneset.result$Aversive.pValue)

###combining results
###find potential duplicated gene set names
duplicated.names <- intersect(rownames(combined.custom.results), rownames(combined.HGNC.results))
dAUROC.2025.final <- rbind(combined.custom.results, combined.HGNC.results[! rownames(combined.HGNC.results) %in% duplicated.names, ])
###further clean of duplicated gene sets
dAUROC.2025.final <- dAUROC.2025.final[! rownames(dAUROC.2025.final) %in% c("Synaptotagmin", "non.clustered.protocadherins", "Glycine_receptors", "EPH_receptors", "Calcium_channels._voltage.dependent"), ]
### save results
write.csv(dAUROC.2025.final, "~/Documents/backup_mac20250207/Dev_manuscript/data.submission/dAUROC.481GeneSets.mean.csv", row.names = T)

