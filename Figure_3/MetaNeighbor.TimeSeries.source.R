library(Seurat)
library(MetaNeighbor)
library(SummarizedExperiment)
library(RColorBrewer)
library(svMisc)
library(ggplot2)
library(dplyr)

# NOTE: all the following functions takes Seurat object as input, 
# which contains the 'Appetitve', 'Aversive' populations in their Idents, and 'stage' column in their metadata.
# For reproducibility, a slow version of 'MetaNeighborUS' is presented here (with fast_version = FALSE). Typical running time: if n.run=10, 15-20min are required for 1 gene set

################################# Function for extract meta neighbor training-testing AUROC scores ##################################

NM_run <- function(genes, 
                   seu.obj, 
                   seu.assay = "SCT", 
                   seu.layer = "data",
                   do.sampling = FALSE,
                   split.data = NULL
                   # do.sampling: if set TRUE, the function will randomly sample training and test cells. Then split.data arguments needs to be set.
                   # do.sampling: if set FALSE, seu.obj$sampling metadata must be provided
                   # split.data: if do.sampling = TRUE, the proportion of training v.s. testing cells need to be set. for example: c(0.5, 0.5) means split the data half-half
                   ){
  
  seu.obj <- seu.obj[genes, ]
  
  if (do.sampling) {
    seu.obj$sampling <- sample(c("train", "test"), replace=TRUE, size=ncol(seu.obj), prob = split.data)
    print("Do random sampling inside MN_run")
  #} else {
  #  print("using sampling metadata from seurat object")
  }
  
  MN_data.list <- SplitObject(seu.obj, split.by = "stage")
  
  #pb = txtProgressBar(min = 0, max = length(MN_data.list), initial = 0) 
  AUROC.df <- data.frame(scores=numeric(),
                         cell_type=character(),
                         stage=character())
  
  for (i in 1:length(MN_data.list)) {
    sample_id <- colnames(MN_data.list[[i]])
    study_id.sampling <- as.character(MN_data.list[[i]]$sampling)
    #study_id.batch <- as.character(MN_data.list[[i]]$batch)
    cell_type <- as.character(MN_data.list[[i]]$trajectories)
    data.assay <- GetAssayData(MN_data.list[[i]], assay = seu.assay, layer = seu.layer)
    
    MN.df <- DataFrame("sample_id" = sample_id, 
                       "study_id.sampling" = study_id.sampling, 
                       #"study_id.batch" = study_id.batch,
                       "cell_type" = cell_type)
    
    MN.se <- SummarizedExperiment(assays=SimpleList(counts=data.assay), colData=MN.df)
    
    scores = MetaNeighborUS(var_genes = genes, dat = MN.se, study_id = study_id.sampling, cell_type = cell_type, symmetric_output=FALSE)
    #the training cell types are displayed as columns and the test cell types are displayed as rows
    
    AUROC.df <- rbind(AUROC.df, data.frame(scores=scores["test|Appetitve", "train|Appetitve"],
                                           cell_type="Appetitve",
                                           stage=names(MN_data.list)[i]))
    AUROC.df <- rbind(AUROC.df, data.frame(scores=scores["test|Aversive", "train|Aversive"],
                                           cell_type="Aversive",
                                           stage=names(MN_data.list)[i]))
    #AUROC.df <- rbind(AUROC.df, data.frame(scores=scores["test|Ast region", "train|Ast region"],
    #                                       cell_type="Ast region",
    #                                       stage=names(MN_data.list)[i]))
    #AUROC.df <- rbind(AUROC.df, data.frame(scores=scores["test|ITC", "train|ITC"],
    #                                       cell_type="ITC",
    #                                       stage=names(MN_data.list)[i]))
    
    #setTxtProgressBar(pb,i)
    
  }
  
  #close(pb)
  return(AUROC.df)
  
}



################################# ################################# ################################# 
################### Function for making random chance line and compare with functional gene set ###########################################################

GeneSetTest <- function(genes, 
                        seu.obj, 
                        seu.assay = "SCT", 
                        seu.layer = "data",
                        n.run = 5,
                        split.data = c(0.5, 0.5),
                        chance.use = "ALL"
                        # split.data: 50-50 works the best
                        # chance.use: choose one of the two for drawing ramdon genes, whether use Highly variable genes from the seurat object or use random genes sampled from all genes
                        # if "HVG" is chosen, variable genes must exist in the RNA assay of the seurat object 
){
  
  timecorse.df <- data.frame(scores=numeric(),
                             cell_type=character(),
                             stage=character(),
                             chance=logical(),
                             n_run=numeric())
  
  sampling.list <- replicate(n.run, sample(c("train", "test"), replace=TRUE, size=ncol(seu.obj), prob = split.data))
  
  if (chance.use == "HVG"){
    
    if (seu.assay == "RNA") {
      HVG.genes <- seu.obj@assays$RNA@var.features
    } else {
      seu.obj.list <- SplitObject(seu.obj, split.by = "batch")
      HVG.genes <- SelectIntegrationFeatures(object.list = seu.obj.list, nfeatures = 3000)
      print("selecting HVGs in SCT assay")
    }
    sel.genes <- replicate(n.run, HVG.genes[sample(length(HVG.genes), replace=FALSE, size=sum(genes %in% rownames(seu.obj)) )])
    
  } else {
    sel.genes <- replicate(n.run, rownames(seu.obj)[sample(nrow(seu.obj), replace=FALSE, size=sum(genes %in% rownames(seu.obj)))])
  }
  
  for (j in 1:n.run) {
    
    seu.obj$sampling <- sampling.list[ ,j]
    
    AUROC.df <- NM_run(genes = genes, seu.obj = seu.obj, seu.assay = seu.assay, seu.layer = seu.layer, do.sampling = FALSE)
    AUROC.df$chance <- "FALSE"
    AUROC.df$n_run <- j
    timecorse.df <- rbind(timecorse.df, AUROC.df)
    
    chance.df <- NM_run(genes = sel.genes[,j], seu.obj = seu.obj, seu.assay = seu.assay, seu.layer = seu.layer, do.sampling = FALSE)
    chance.df$chance <- "TRUE"
    chance.df$n_run <- j
    timecorse.df <- rbind(timecorse.df, chance.df)
    
    progress(j, n.run)
  }
  
  timecorse.df$stage <- factor(timecorse.df$stage, levels=c("E15", "E18", "P0", "P4", "P10", "P21"))
  #timecorse.df <- timecorse.df %>% mutate(real_time = recode(stage, "E15" = -5,
  #                                              "E18" = -2,
  #                                              "P0" = 0,
  #                                              "P4" = 4,
  #                                              "P10" = 10,
  #                                              "P21" = 21)) 
  return(timecorse.df)
  
}

### computing time: if n.run=10, 20min are required for 1 gene set



################################# ################################# ################################# 
####### Function that test the effect of gene numbers on decoding scores.

GeneNumberTest <- function(x,  # A numerical vector that decide the randomly chosen genes to be tested
                           seu.obj, 
                           seu.assay = "SCT",            # either RNA or SCT
                           seu.layer = "data",
                           n.run = 10,
                           split.data = c(0.8, 0.2),
                           chance.use = "ALL"           # either "ALL" or "HVG"
                           # split.data: 50-50 works the best
                           # chance.use: choose one of the two for drawing ramdon genes, whether use Highly variable genes from the seurat object or use random genes sampled from all genes
                           # if "HVG" is chosen, variable genes must exist in the RNA assay of the seurat object 
){
  
  timecorse.df <- data.frame(scores=numeric(),
                             cell_type=character(),
                             stage=character(),
                             ngene=numeric(),
                             n_run=numeric())
  
  sampling.list <- replicate(n.run, sample(c("train", "test"), replace=TRUE, size=ncol(seu.obj), prob = split.data))
  
  if (chance.use == "HVG"){
    
    if (seu.assay == "RNA") {
      HVG.genes <- seu.obj@assays$RNA@var.features
      print("Using HVGs in RNA assay")
    } else {
      seu.obj.list <- SplitObject(seu.obj, split.by = "batch")
      HVG.genes <- SelectIntegrationFeatures(object.list = seu.obj.list, nfeatures = 3000)
      print("selecting HVGs in SCT assay")
    }
    genes <- HVG.genes
    
  } else {
    genes <- rownames(seu.obj)
    print("Using all genes")
  }
  
  for (k in x) {
    sel.genes <- replicate(n.run, genes[sample(nrow(seu.obj), replace=FALSE, size=k)])
    
    for (j in 1:n.run) {
      
      seu.obj$sampling <- sampling.list[ ,j]
      
      chance.df <- NM_run(genes = sel.genes[,j], seu.obj = seu.obj, seu.assay = seu.assay, seu.layer = seu.layer, do.sampling = FALSE)
      chance.df$ngene <- k
      chance.df$n_run <- j
      timecorse.df <- rbind(timecorse.df, chance.df)
      
      #progress(j, n.run)
      message(paste0("round", j, " of ", n.run))
    }
    
    message(paste0(which(x==k), " in ",length(x)))
    if (k == x[length(x)]) message("Done!")
  }
  
  timecorse.df$stage <- factor(timecorse.df$stage, levels=c("E15", "E18", "P0", "P4", "P10", "P21"))
  return(timecorse.df)
  
}


################################# ################################# ################################# 
####### Function that test a list of gene sets, normalize the AUROC score with the chance line and output a matrix for heat map

### CAUTION: it is important that if the seu.assay = "SCT", the default assay of the seurat object should also set to "SCT"
### This ensures the correct number of genes are calculated for random sampling of control genes

GeneSetTestMatrix <- function(Geneset_list,  # A list of vectors, each vector is a gene set that's going to be tested, element name should be the gene set name
                           seu.obj, 
                           seu.assay = "SCT",            # either RNA or SCT
                           seu.layer = "data",
                           n.run = 10,
                           split.data = c(0.8, 0.2),
                           chance.use = "ALL"           # either "ALL" or "HVG"
                           # chance.use: choose one of the two for drawing ramdon genes, whether use Highly variable genes from the seurat object or use random genes sampled from all genes
                           # if "HVG" is chosen, variable genes must exist in the RNA assay of the seurat object 
){
  
  # sample the test vs training datasets for each of the reduplicate run 
  sampling.list <- replicate(n.run, sample(c("train", "test"), replace=TRUE, size=ncol(seu.obj), prob = split.data))
  
  if (chance.use == "HVG"){
    
    if (seu.assay == "RNA") {
      HVG.genes <- seu.obj@assays$RNA@var.features
      print("Using HVGs in RNA assay")
    } else {
      seu.obj.list <- SplitObject(seu.obj, split.by = "batch")
      HVG.genes <- SelectIntegrationFeatures(object.list = seu.obj.list, nfeatures = 3000)
      print("selecting HVGs in SCT assay")
    }
    genes <- HVG.genes
    
  } else {
    genes <- rownames(seu.obj)
    print("Using all genes")
  }
  
  # define the two output matrix/dataframe
  Appetitve.matrix <- data.frame(E15=numeric(),
                                 E18=numeric(),
                                 P0=numeric(),
                                 P4=numeric(),
                                 P10=numeric(),
                                 P21=numeric())
  
  Aversive.matrix <- data.frame(E15=numeric(),
                                 E18=numeric(),
                                 P0=numeric(),
                                 P4=numeric(),
                                 P10=numeric(),
                                 P21=numeric())
  
  for (k in 1:length(Geneset_list)) {
    
    if (sum(Geneset_list[[k]] %in% rownames(seu.obj)) < 3) {
      # Output a warning message
      message(paste("Less than 3 genes found in", names(Geneset_list)[k]))
      # If condition is not met, skip to the next iteration
      next
    }

    # sample genes for chance line
    sel.genes <- replicate(n.run, 
                           genes[sample(nrow(seu.obj), replace=FALSE, size=sum(Geneset_list[[k]] %in% rownames(seu.obj)))]
                           )
    
    set.genes <- Geneset_list[[k]]
    
    # define a dataframe for the output of each gene set
    timecorse.df <- data.frame(scores=numeric(),
                               cell_type=character(),
                               stage=character(),
                               ngene=numeric(),
                               n_run=numeric())
    
    for (j in 1:n.run) {
      seu.obj$sampling <- sampling.list[ ,j]
      
      chance.df <- NM_run(genes = sel.genes[,j], seu.obj = seu.obj, seu.assay = seu.assay, seu.layer = seu.layer, do.sampling = FALSE)
      chance.df$chance <- "TURE"
      chance.df$n_run <- j
      timecorse.df <- rbind(timecorse.df, chance.df)
      
      AUROC.df <- NM_run(genes = set.genes, seu.obj = seu.obj, seu.assay = seu.assay, seu.layer = seu.layer, do.sampling = FALSE)
      AUROC.df$chance <- "FALSE"
      AUROC.df$n_run <- j
      timecorse.df <- rbind(timecorse.df, AUROC.df)
      
      #progress(j, n.run)
      message(paste0("round", j, " of ", n.run))
    }
    
    timecorse.df %>% group_by(cell_type, stage, chance) %>% summarise(mean = mean(scores)) %>% as.data.frame() -> timecorse.reshape
    
    timecorse.reshape <- reshape(timecorse.reshape, idvar = c("cell_type", "stage"), timevar = "chance", direction = "wide")
    timecorse.reshape %>% mutate(dAUROC=mean.FALSE-mean.TURE) -> timecorse.reshape
    timecorse.reshape$stage <- factor(timecorse.reshape$stage, levels=c("E15", "E18", "P0", "P4", "P10", "P21"))
    timecorse.reshape %>% arrange(stage, .by_group = TRUE) -> timecorse.reshape
    
    timecorse.reshape %>% filter(cell_type == "Appetitve") -> Appetitve.scores
    Appetitve.scores <- t(Appetitve.scores[, c("stage", "dAUROC")])
    Appetitve.matrix[k, ] <- as.numeric(Appetitve.scores["dAUROC", ])
    rownames(Appetitve.matrix)[k] <- names(Geneset_list)[k]
    
    timecorse.reshape %>% filter(cell_type == "Aversive") -> Aversive.scores
    Aversive.scores <- t(Aversive.scores[, c("stage", "dAUROC")])
    Aversive.matrix[k, ] <- as.numeric(Aversive.scores["dAUROC", ])
    rownames(Aversive.matrix)[k] <- names(Geneset_list)[k]
    
    message(paste0(k, " in ", length(Geneset_list)))
    if (k == length(Geneset_list)) message("Done!")
  }
  
  dAUROC.scores <- list("Appetitive" = Appetitve.matrix,
                        "Aversive" = Aversive.matrix)
  return(dAUROC.scores)
  
}



####### Function that test a list of gene sets, return both the delta-AUROC scores and p-value comparing the paired gene set and randomly selected genes
### CAUTION: it is important that if the seu.assay = "SCT", the default assay of the seurat object should also set to "SCT"
### This ensures the correct number of genes are calculated for random sampling of control genes

GeneSetTestMatrix2 <- function(Geneset_list,  # A list of vectors, each vector is a gene set that's going to be tested, element name should be the gene set name
                              seu.obj, 
                              seu.assay = "SCT",            # either RNA or SCT
                              seu.layer = "data",
                              n.run = 10,
                              split.data = c(0.8, 0.2),
                              chance.use = "ALL"           # either "ALL" or "HVG"
                              # chance.use: choose one of the two for drawing ramdon genes, whether use Highly variable genes from the seurat object or use random genes sampled from all genes
                              # if "HVG" is chosen, variable genes must exist in the RNA assay of the seurat object 
){
  
  # sample the test vs training datasets for each of the reduplicate run 
  sampling.list <- replicate(n.run, sample(c("train", "test"), replace=TRUE, size=ncol(seu.obj), prob = split.data))
  
  if (chance.use == "HVG"){
    
    if (seu.assay == "RNA") {
      HVG.genes <- seu.obj@assays$RNA@var.features
      print("Using HVGs in RNA assay")
    } else {
      seu.obj.list <- SplitObject(seu.obj, split.by = "batch")
      HVG.genes <- SelectIntegrationFeatures(object.list = seu.obj.list, nfeatures = 3000)
      print("selecting HVGs in SCT assay")
    }
    genes <- HVG.genes
    
  } else {
    genes <- rownames(seu.obj)
    print("Using all genes")
  }
  
  # define the two output matrix/dataframe
  Appetitve.matrix <- data.frame(E15=numeric(),
                                 E18=numeric(),
                                 P0=numeric(),
                                 P4=numeric(),
                                 P10=numeric(),
                                 P21=numeric())
  
  Aversive.matrix <- data.frame(E15=numeric(),
                                E18=numeric(),
                                P0=numeric(),
                                P4=numeric(),
                                P10=numeric(),
                                P21=numeric())
  
  Appetitve.p.matrix <- data.frame(E15=numeric(),
                                 E18=numeric(),
                                 P0=numeric(),
                                 P4=numeric(),
                                 P10=numeric(),
                                 P21=numeric())
  
  Apersive.p.matrix <- data.frame(E15=numeric(),
                                E18=numeric(),
                                P0=numeric(),
                                P4=numeric(),
                                P10=numeric(),
                                P21=numeric())
  
  for (k in 1:length(Geneset_list)) {
    
    if (sum(Geneset_list[[k]] %in% rownames(seu.obj)) < 3) {
      # Output a warning message
      message(paste("Less than 3 genes found in", names(Geneset_list)[k]))
      # If condition is not met, skip to the next iteration
      next
    }
    
    # sample genes for chance line
    sel.genes <- replicate(n.run, 
                           genes[sample(nrow(seu.obj), replace=FALSE, size=sum(Geneset_list[[k]] %in% rownames(seu.obj)))]
    )
    
    set.genes <- Geneset_list[[k]]
    
    # define a dataframe for the output of each gene set
    timecorse.df <- data.frame(scores=numeric(),
                               cell_type=character(),
                               stage=character(),
                               ngene=numeric(),
                               n_run=numeric())
    
    for (j in 1:n.run) {
      seu.obj$sampling <- sampling.list[ ,j]
      
      chance.df <- NM_run(genes = sel.genes[,j], seu.obj = seu.obj, seu.assay = seu.assay, seu.layer = seu.layer, do.sampling = FALSE)
      chance.df$chance <- "TURE"
      chance.df$n_run <- j
      timecorse.df <- rbind(timecorse.df, chance.df)
      
      AUROC.df <- NM_run(genes = set.genes, seu.obj = seu.obj, seu.assay = seu.assay, seu.layer = seu.layer, do.sampling = FALSE)
      AUROC.df$chance <- "FALSE"
      AUROC.df$n_run <- j
      timecorse.df <- rbind(timecorse.df, AUROC.df)
      
      #progress(j, n.run)
      message(paste0("round", j, " of ", n.run))
    }
    
    # Calculate delta-AUROC
    
    timecorse.df %>% group_by(cell_type, stage, chance) %>% summarise(mean = mean(scores)) %>% as.data.frame() -> timecorse.reshape
    
    timecorse.reshape <- reshape(timecorse.reshape, idvar = c("cell_type", "stage"), timevar = "chance", direction = "wide")
    timecorse.reshape %>% mutate(dAUROC=mean.FALSE-mean.TURE) -> timecorse.reshape
    timecorse.reshape$stage <- factor(timecorse.reshape$stage, levels=c("E15", "E18", "P0", "P4", "P10", "P21"))
    timecorse.reshape %>% arrange(stage, .by_group = TRUE) -> timecorse.reshape
    
    timecorse.reshape %>% filter(cell_type == "Appetitve") -> Appetitve.scores
    Appetitve.scores <- t(Appetitve.scores[, c("stage", "dAUROC")])
    Appetitve.matrix[k, ] <- as.numeric(Appetitve.scores["dAUROC", ])
    rownames(Appetitve.matrix)[k] <- names(Geneset_list)[k]
    
    timecorse.reshape %>% filter(cell_type == "Aversive") -> Aversive.scores
    Aversive.scores <- t(Aversive.scores[, c("stage", "dAUROC")])
    Aversive.matrix[k, ] <- as.numeric(Aversive.scores["dAUROC", ])
    rownames(Aversive.matrix)[k] <- names(Geneset_list)[k]
    
    # Calculate p-value. using t-test by default
    
    timecorse.pvalue <- timecorse.df %>%
      group_by(cell_type, stage) %>%
      summarize(p_value = t.test(scores ~ chance)$p.value)
    timecorse.pvalue$stage <- factor(timecorse.pvalue$stage, levels=c("E15", "E18", "P0", "P4", "P10", "P21"))
    timecorse.pvalue %>% arrange(stage, .by_group = TRUE) -> timecorse.pvalue
    
    timecorse.pvalue %>% filter(cell_type == "Appetitve") -> Appetitve.pvalue
    Appetitve.pvalue <- t(Appetitve.pvalue[, c("stage", "p_value")])
    Appetitve.p.matrix[k, ] <- as.numeric(Appetitve.pvalue["p_value", ])
    rownames(Appetitve.p.matrix)[k] <- names(Geneset_list)[k]
    
    timecorse.pvalue %>% filter(cell_type == "Aversive") -> Aversive.pvalue
    Aversive.pvalue <- t(Aversive.pvalue[, c("stage", "p_value")])
    Apersive.p.matrix[k, ] <- as.numeric(Aversive.pvalue["p_value", ])
    rownames(Apersive.p.matrix)[k] <- names(Geneset_list)[k]
    
    message(paste0(k, " in ", length(Geneset_list)))
    if (k == length(Geneset_list)) message("Done!")
  }
  
  dAUROC.scores <- list("Appetitive.mean" = Appetitve.matrix,
                        "Aversive.mean" = Aversive.matrix,
                        "Appetitive.pValue" = Appetitve.p.matrix,
                        "Aversive.pValue" = Apersive.p.matrix)
  return(dAUROC.scores)
  
}








