library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(clue)
library(MetaNeighbor)
library(SingleCellExperiment)
library(pbapply)
library(future)
options(future.globals.maxSize = 8 * 1024^3)  # Set max size to 8GB
plan("multicore")  # Or "multisession" depending on your OS

######################################################### Build summarized experiment objects for Metaneighbour
### Prepare mouse developmental P4 datasets (gene name converted into human)
mouse.dev.se <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.6/clean.human.mouseDev.merged.rds")

### recalculate quality matrics for SCT normalization
mouse.dev.se$nCount_RNA <- Matrix::colSums(mouse.dev.se[["RNA"]]$counts)
mouse.dev.se$nFeature_RNA <- Matrix::colSums(mouse.dev.se[["RNA"]]$counts > 0)
mouse.dev.se[["percent.mito"]] <- PercentageFeatureSet(mouse.dev.se, pattern = "^MT-")
### Normalization
DefaultAssay(mouse.dev.se) <- "RNA"
mouse.P4.se <- subset(mouse.dev.se, subset = stage == "P4")
mouse.P4.list <- SplitObject(mouse.P4.se, split.by = "batch")
mouse.P4.list <- lapply(X = mouse.P4.list, FUN = function(x) {
  SCTransform(x, vst.flavor = "v2", vars.to.regress = c("percent.mito", "nFeature_RNA", "nCount_RNA"))
})
# regress out percent.mito and nFeature_RNA is important, regress out nCount_RNA also helps
mouse.P4.list <- lapply(X = mouse.P4.list, FUN = function(x) {
  subset(x, trajectories %in% c("Aversive", "Ast region", "Appetitve", "IPAC", "ITC"))
})


### Prepare human annotated datasets 
human.se <- readRDS("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.6/clean.human.se.consensus.rds")
### recalculate quality matrics for SCT normalization
human.se$nCount_RNA <- Matrix::colSums(human.se[["RNA"]]$counts)
human.se$nFeature_RNA <- Matrix::colSums(human.se[["RNA"]]$counts > 0)
human.se[["percent.mito"]] <- PercentageFeatureSet(human.se, pattern = "^MT-")

### Add metadata
human.se$trajectories <- human.se$human.mouse.consensus
human.se$stage <- "human"

### Normalization
DefaultAssay(human.se) <- "RNA"
human.list <- SplitObject(human.se, split.by = "donor_id")
human.list <- lapply(X = human.list, FUN = function(x) {
  SCTransform(x, vst.flavor = "v2", vars.to.regress = c("percent.mito", "nFeature_RNA", "nCount_RNA"))
})
human.list <- lapply(X = human.list, FUN = function(x) {
  subset(x, trajectories %in% c("Aversive", "Ast region", "Appetitve", "IPAC", "ITC"))
})



### merge human-mouse.P4 data sets
human.P4.list <- c(mouse.P4.list, human.list)

human.P4.list <- lapply(human.P4.list, function(obj) {
  if ("integrated" %in% names(obj@assays)) {
    obj[["integrated"]] <- NULL
  }
  return(obj)
})

human.P4.sct <- merge(
  human.P4.list[[1]],
  y = human.P4.list[2:length(human.P4.list)]
)


# Exclude genes with "-AS" from a list
filtered_genes <- grep("-AS", rownames(human.P4.sct), invert = TRUE, value = TRUE)
human.P4.sct <- human.P4.sct[filtered_genes,]



### build summarized experiment objects and define inputs for MetaNeighbor
hmP4.sample.id <- human.P4.sct$stage
hmP4.trajectory.id <- human.P4.sct$trajectories

hmP4.assay.sct <- GetAssayData(human.P4.sct, assay = "SCT", layer = "data")    # SCT makes so much more sense than RNA assay
hmP4.MN <- SummarizedExperiment(assays = hmP4.assay.sct)

### save object
#saveRDS(hmP4.MN, "h.mP4.SummarizedExp.rds")

######################################################### Function to test different gene set sizes on MetaNeighbor and extract AUROCs

test_gene_set_size <- function(gene_set_sizes, SE_object, study_id, cell_type,
                               source_label = "P4|Appetitve", target_label = "human|Appetitve",
                               source_label_2 = "P4|Aversive", target_label_2 = "human|Aversive",
                               n_repeats = 3, seed = 42) {
  set.seed(seed)
  all_genes <- rownames(SE_object)
  
  # Create all test combinations
  param_grid <- expand.grid(
    size = gene_set_sizes,
    replicate = seq_len(n_repeats)
  )
  
  # Wrap loop with progress bar
  results <- pbapply::pbapply(param_grid, 1, function(row) {
    size <- row[1]
    i <- row[2]
    
    sampled_genes <- sample(all_genes, size)
    
    auroc_mat <- tryCatch({
      MetaNeighborUS(
        var_genes = sampled_genes,
        dat = SE_object,
        study_id = study_id,
        cell_type = cell_type,
        symmetric_output = TRUE,
        one_vs_best = FALSE,
        fast_version = TRUE
      )
    }, error = function(e) return(NULL))
    
    if (is.null(auroc_mat)) return(NULL)
    
    appetitive_score <- tryCatch({
      auroc_mat[source_label, target_label]
    }, error = function(e) NA)
    
    aversive_score <- tryCatch({
      auroc_mat[source_label_2, target_label_2]
    }, error = function(e) NA)
    
    data.frame(
      gene_set_size = size,
      replicate = i,
      appetitive_AUROC = appetitive_score,
      aversive_AUROC = aversive_score
    )
  })
  
  # Combine all into one dataframe
  results_df <- do.call(rbind, Filter(Negate(is.null), results))
  
  n_total <- nrow(param_grid)
  n_success <- length(Filter(Negate(is.null), results))
  message(sprintf("Finished: %d of %d successful (%d dropped due to errors)",
                  n_success, n_total, n_total - n_success))
  
  results_df <- na.omit(results_df)
  
  return(results_df)
  
}


### Test of gene set size on prediction
gene_set_sizes <- c(5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 350, 400, 450, 500)

auroc_results <- test_gene_set_size(
  gene_set_sizes = gene_set_sizes,
  SE_object = hmP4.MN,
  study_id = hmP4.sample.id,
  cell_type = hmP4.trajectory.id,
  n_repeats = 200  # more than 50 iterations yeild robust estimates 
)


### Visualization of the effects of gene set size
ggplot(auroc_results, aes(x = gene_set_size)) +
  geom_point(aes(y = appetitive_AUROC, color = "Appetitive")) +
  geom_point(aes(y = aversive_AUROC, color = "Aversive")) +
  geom_smooth(aes(y = appetitive_AUROC, color = "Appetitive"), method = "loess", se = FALSE) +
  geom_smooth(aes(y = aversive_AUROC, color = "Aversive"), method = "loess", se = FALSE) +
  labs(title = "Effect of Gene Set Size on AUROC",
       y = "AUROC Score",
       x = "Number of Genes") +
  theme_minimal()


################################################################ look at individual gene families -- Fig. 6F
### Get Semaphorins and Neuropeptides gene family
Semaphorin.genes <- grep(pattern = "^SEMA", x = rownames(hmP4.assay.sct), value = TRUE)
Neuropeptide.genes <- c("ADCYAP1", "ADM",     "AGRP",    "AVP",     "BDNF",    "BMP2",    "BMP4",    "CALCB",   "CARTPT",  "CCK",     "CORT",    "CRH",    
 "EDN1",    "EDN3",    "FGF9",    "GAL",     "GHRH",    "GRP",     "HCRT",    "IGF1",    "INHBA",   "INHBB",   "KISS1",   "NMB",    
 "NMS",     "NMU",     "NPPC",    "NPVF",    "NPW",     "NPY",     "NRTN",    "NTS",     "OXT",     "PDYN",    "PENK",    "PNOC",   
 "POMC",    "PROK2",   "PTHLH",   "PTN",     "QRFP",    "RLN1",    "RSPO1",   "RSPO2",   "RXFP1",   "SST",     "TAC1",    "TAC3",   
 "TGFB3",   "TRH",     "TSHB",    "VIP",     "WNT2",    "WNT4",    "WNT5A",   "UCN3",    "CRHBP")  

### Run null distributions for Semaphorins (19 genes) and Neuropeptides (57 genes)
Null.sema.neuropeptides <- test_gene_set_size(
  gene_set_sizes = c(57, 19),
  SE_object = hmP4.MN,
  study_id = hmP4.sample.id,
  cell_type = hmP4.trajectory.id,
  n_repeats = 1000  # more than 50 iterations yeild robust estimates 
)

### Getting null distributions for Semaphorins (19 genes) and Neuropeptides (57 genes)
null_semaphorins <- Null.sema.neuropeptides$appetitive_AUROC[Null.sema.neuropeptides$gene_set_size == 19]
null_neuropeptided <- Null.sema.neuropeptides$appetitive_AUROC[Null.sema.neuropeptides$gene_set_size == 57]


### Get semaphorin AUROC
sema_auroc_app <- MetaNeighborUS(
  var_genes = Semaphorin.genes,
  dat = hmP4.MN,
  study_id = hmP4.sample.id,
  cell_type = hmP4.trajectory.id,
  symmetric_output = TRUE,
  one_vs_best = FALSE,
  fast_version = TRUE
)["P4|Appetitve", "human|Appetitve"]

sema_auroc_ave <- MetaNeighborUS(
  var_genes = Semaphorin.genes,
  dat = hmP4.MN,
  study_id = hmP4.sample.id,
  cell_type = hmP4.trajectory.id,
  symmetric_output = TRUE,
  one_vs_best = FALSE,
  fast_version = TRUE
)["P4|Aversive", "human|Aversive"]

### Compute Z-score
(sema_auroc_ave - mean(null_semaphorins)) / sd(null_semaphorins)

### Compute empirical p-value
mean(null_semaphorins >= sema_auroc_ave)

### Visualization
hist(null_semaphorins, breaks = 50, col = "grey", main = "Null AUROC Distribution (19 genes)",
     xlab = "AUROC", xlim = c(0.2, 0.9), freq = TRUE, density = 70)
abline(v = mean(null_semaphorins), col = "gray50", lwd = 1)
abline(v = sema_auroc_app, col = "#377EB8", lwd = 4)
abline(v = sema_auroc_ave, col = "#E41A1C", lwd = 4)



### Get Neuropeptide AUROC
Neuropeptide_auroc_app <- MetaNeighborUS(
  var_genes = Neuropeptide.genes,
  dat = hmP4.MN,
  study_id = hmP4.sample.id,
  cell_type = hmP4.trajectory.id,
  symmetric_output = TRUE,
  one_vs_best = FALSE,
  fast_version = TRUE
)["P4|Appetitve", "human|Appetitve"]

Neuropeptide_auroc_ave <- MetaNeighborUS(
  var_genes = Neuropeptide.genes,
  dat = hmP4.MN,
  study_id = hmP4.sample.id,
  cell_type = hmP4.trajectory.id,
  symmetric_output = TRUE,
  one_vs_best = FALSE,
  fast_version = TRUE
)["P4|Aversive", "human|Aversive"]

# Compute Z-score
(Neuropeptide_auroc_app - mean(null_neuropeptided)) / sd(null_neuropeptided)

# Compute empirical p-value
mean(null_neuropeptided >= Neuropeptide_auroc_app)

### Visualization
hist(null_neuropeptided, breaks = 50, col = "grey", main = "Null AUROC Distribution (57 genes)",
     xlab = "AUROC", xlim = c(0.2, 0.9), freq = TRUE, density = 70)
abline(v = mean(null_neuropeptided), col = "gray50", lwd = 1)
abline(v = Neuropeptide_auroc_app, col = "#377EB8", lwd = 4)
abline(v = Neuropeptide_auroc_ave, col = "#E41A1C", lwd = 4)



################################################################## Prepare gene families for Human-mouse conserved gene family screening
### read gene families (same list as in Fig.3)
library(readxl)
custom.data <- read_excel("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.3/Final2025.Families.axon.synapse.Rinput.xlsx")
custom.data <- lapply(custom.data, as.character)
names(custom.data) <- make.names(names(custom.data))
custom.data <- lapply(custom.data, function(x) x[!is.na(x)])


HGNC.list <- read_excel("~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.3/HGNC_genes.Rinput.xlsx")
HGNC.list <- lapply(HGNC.list, as.character)
names(HGNC.list) <- make.names(names(HGNC.list))
HGNC.list <- lapply(HGNC.list, function(x) x[!is.na(x)])
HGNC.50.GeneList <- HGNC.list[sapply(HGNC.list, function(x) length(x) <= 50, simplify = TRUE)]

#### transfer mouse to human gene names (Only needed once)
if (!requireNamespace("homologene", quietly = TRUE)) {
  install.packages("homologene")
}
library(homologene)

map_mouse_to_human_genes <- function(mouse_genes) {
  # Query homologene
  result <- homologene::mouse2human(mouse_genes)
  
  # Return human symbols (unique)
  return(unique(result$humanGene))
}

### apply to the two gene list
custom.data.human <- lapply(custom.data, function(genes) {
  map_mouse_to_human_genes(unique(genes))
})

HGNC.data.human <- lapply(HGNC.50.GeneList, function(genes) {
  map_mouse_to_human_genes(unique(genes))
})


### clean gene name list
custom.data.human <- custom.data.human[!names(custom.data.human) %in% c("cytoskeleton.binding...8", "Calcium.channel.regulatory.subunits", "Semaphorins")]

gene.family.all <- c(custom.data.human, HGNC.data.human)
### If further duplecated, keep only the first occurrence of each name
gene.family.all <- gene.family.all[!duplicated(names(gene.family.all))]


########################################################## Large-scale screen of conserved gene families -- Fig. 6G
### build function
compute_metaneighbor_zscores <- function(gene_sets, SE_object, study_id, cell_type,
                                         source_label = "P4|Appetitve", target_label = "human|Appetitve",
                                         source_label_2 = "P4|Aversive", target_label_2 = "human|Aversive",
                                         n_null = 200, seed = 42) {
  set.seed(seed)
  all_genes <- rownames(SE_object)
  
  gene_set_names <- names(gene_sets)
  
  # Apply with progress bar
  result_list <- pbapply::pblapply(gene_set_names, function(gene_set_name) {
    gene_set <- intersect(gene_sets[[gene_set_name]], all_genes)
    gene_set_size <- length(gene_set)
    
    if (gene_set_size < 4) {
      warning(paste("Skipping", gene_set_name, "- fewer than 4 matching genes"))
      return(NULL)
    }
    
    # Observed AUROC
    observed_scores <- tryCatch({
      MetaNeighborUS(
        var_genes = gene_set,
        dat = SE_object,
        study_id = study_id,
        cell_type = cell_type,
        symmetric_output = TRUE,
        one_vs_best = FALSE,
        fast_version = TRUE
      )
    }, error = function(e) NULL)
    
    if (is.null(observed_scores)) return(NULL)
    
    obs_app <- tryCatch(observed_scores[source_label, target_label], error = function(e) NA)
    obs_avs <- tryCatch(observed_scores[source_label_2, target_label_2], error = function(e) NA)
    
    # Null AUROC distributions
    null_app <- numeric(n_null)
    null_avs <- numeric(n_null)
    
    for (i in seq_len(n_null)) {
      sampled_genes <- sample(all_genes, gene_set_size)
      
      null_scores <- tryCatch({
        MetaNeighborUS(
          var_genes = sampled_genes,
          dat = SE_object,
          study_id = study_id,
          cell_type = cell_type,
          symmetric_output = TRUE,
          one_vs_best = FALSE,
          fast_version = TRUE
        )
      }, error = function(e) NULL)
      
      null_app[i] <- if (!is.null(null_scores)) tryCatch(null_scores[source_label, target_label], error = function(e) NA) else NA
      null_avs[i] <- if (!is.null(null_scores)) tryCatch(null_scores[source_label_2, target_label_2], error = function(e) NA) else NA
    }
    
    # Clean NAs
    null_app <- na.omit(null_app)
    null_avs <- na.omit(null_avs)
    
    if (length(null_app) < 5 || length(null_avs) < 5) {
      warning(paste("Too few valid null scores for", gene_set_name))
      return(NULL)
    }
    
    # Stats
    z_app <- (obs_app - mean(null_app)) / sd(null_app)
    z_avs <- (obs_avs - mean(null_avs)) / sd(null_avs)
    p_app <- mean(null_app >= obs_app)
    p_avs <- mean(null_avs >= obs_avs)
    delta_app <- obs_app - mean(null_app)
    delta_avs <- obs_avs - mean(null_avs)
    
    data.frame(
      gene_set = gene_set_name,
      n_genes = gene_set_size,
      appetitive_AUROC = obs_app,
      appetitive_null_mean = mean(null_app),
      appetitive_null_sd = sd(null_app),
      appetitive_z = z_app,
      appetitive_p = p_app,
      appetitive_delta = delta_app,
      aversive_AUROC = obs_avs,
      aversive_null_mean = mean(null_avs),
      aversive_null_sd = sd(null_avs),
      aversive_z = z_avs,
      aversive_p = p_avs,
      aversive_delta = delta_avs
    )
  })
  
  # Combine non-null results
  results_df <- do.call(rbind, Filter(Negate(is.null), result_list))
  rownames(results_df) <- NULL
  return(results_df)
}


### Apply function to merged dataset screening for ALL gene families
hmP4.z.results <- compute_metaneighbor_zscores(
  gene_sets = gene.family.all,
  SE_object = hmP4.MN,
  study_id = hmP4.sample.id,
  cell_type = hmP4.trajectory.id,
  source_label = "P4|Appetitve", target_label = "human|Appetitve",
  source_label_2 = "P4|Aversive", target_label_2 = "human|Aversive",
  n_null = 200  # or more for smoother estimates
)


### Clean for long names
name_map <- c(
  "Neuropeptides..CRH.BP" = "Neuropeptides",
  "Potassium_channels._voltage.gated" = "VGKCs",
  "PRD_class_homeoboxes_and_pseudogenes" = "PRD_homeoboxes",
  "SRY_.sex_determining_region_Y..boxes" = "SRY_family_TFs",
  "ZF_class_homeoboxes_and_pseudogenes" = "ZF_homeoboxes",
  "Calcium.calmodulin.dependent.protein.kinase.type.II" = "CAMKIIa",
  "A.kinase_anchoring_proteins" = "PKA_anchoring",
  "LIM_class_homeoboxes" = "LIM_homeoboxes",
  "SV2.family.and.related" = "SV2.family.and.related"
)


hmP4.z.results <- hmP4.z.results %>%
  mutate(gene_set = recode(gene_set, !!!name_map))

hmP4.z.results <- hmP4.z.results %>%
  filter(gene_set != "Synaptotagmin")

#write.csv(hmP4.z.results, "~/Documents/backup_mac20250207/Dev_manuscript/data.submission/fig.6/human.mouse.conserved.families.csv")
