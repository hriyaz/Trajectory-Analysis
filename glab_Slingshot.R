glab_Slingshot <- function(expression, meta, meta_trait, start.clus, nfeatures, npcs){

  suppressPackageStartupMessages({
    library(Seurat); library(slingshot); library(SingleCellExperiment)
    library(RColorBrewer); library(ggrepel); library(tidyverse); library(ggplot2)
  })

  #Create seurat object
  seur_obj <- CreateSeuratObject(counts = expression)

  #Variable Features
  seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = nfeatures)
  var_genes <- VariableFeatures(seur_obj)
  seur_obj <- seur_obj[rownames(seur_obj$RNA) %in% var_genes,]

  #PCA
  seur_obj <- ScaleData(seur_obj)
  seur_obj <- RunPCA(seur_obj, npcs = npcs)

  #Convert to SCE
  sce_obj <- as.SingleCellExperiment(seur_obj)
  rd <- as.data.frame(reducedDim(sce_obj))

  #Run Slingshot
  clusters <- meta[match(rownames(rd), rownames(meta)),meta_trait]
  sling_obj <- slingshot(sce_obj,  reducedDim = 'PCA', clusterLabels = clusters, start.clus = start.clus)

  #Curves
  curves <- slingCurves(sling_obj, as.df = TRUE)
  MST <- slingMST(sling_obj, as.df = TRUE)

  rd <<- rd
  clusters <<- clusters

  curves <<- curves
  MST <<- MST

  sds <<- SlingshotDataSet(sling_obj)
  meta <<- meta
  meta_trait <<- meta_trait
}
