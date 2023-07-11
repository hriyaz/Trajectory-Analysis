SlingshotNOCLUSTERINFO <- function(expression, meta, meta_trait, start.clus, nfeatures, npcs){

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

plot_SlingshotNOCLUSTERINFO <- function(rd, clusters, curves, MST, PCX1, PCX2, text_size){

#Run ggplot
normal_PCA <- ggplot(rd, aes_string(x = colnames(rd)[PCX1], y = colnames(rd)[PCX2])) + geom_point() 


sling_curve <- normal_PCA + geom_path(data = curves %>% arrange(Order), aes(group = Lineage), linewidth = 1, color = "gray") + geom_point(aes(color = as.factor(clusters)))
sling_curve_L <- normal_PCA + geom_path(data = curves %>% arrange(Order), aes(group = Lineage), linewidth = 1, color = "gray") + geom_point(aes(color = as.factor(clusters))) + geom_text_repel(label = row.names(rd), size = text_size)

print(sling_curve)
print(sling_curve_L)

sling_MST <- normal_PCA + geom_path(data = MST %>% arrange(Order), aes(group = Lineage), linewidth = 1, color = "gray") + geom_point(data = MST, size = 2.5) + geom_point(aes(color = as.factor(clusters)))
sling_MST_L <- normal_PCA + geom_path(data = MST %>% arrange(Order), aes(group = Lineage), linewidth = 1, color = "gray") + geom_point(data = MST, size = 2.5) + geom_point(aes(color = as.factor(clusters))) + geom_text_repel(label = row.names(rd), size = text_size)

print(sling_MST)
print(sling_MST_L)

#lin1 <- getLineages(as.data.frame(reducedDim(sds)), clusterLabels = slingshot_clusternames)  #sds
#crv1 <- getCurves(lin1)

g <- ggplot_build(sling_curve)
ggplot_colors <<- unique(g$data[[1]]["fill"])

}

plot3d_SlingshotNOCLUSTERINFO <- function(sds, meta, meta_trait){
  library(rgl)
  sds_redctns <- reducedDim(sds)
  groupcolor_ind <- match(rownames(sds_redctns), rownames(meta))
  clusters <- meta[match(rownames(sds_redctns), rownames(meta)),meta_trait]
  
  group_colors <- c("red", "green3", "blue", "purple", "black", "darkkhaki")
  color_vector <- group_colors[match(clusters, unique(clusters))]
  
  plot3d(reducedDim(sds), size = 4, colvar = clusters, col = color_vector)#, col = meta$Group) 
  plot3d.SlingshotDataSet(sds, lwd = 3, add = TRUE)
  text3d(sds_redctns[,1], sds_redctns[,2], sds_redctns[,3], texts = rownames(sds_redctns), add = TRUE, cex = 0.65, adj = c(1, 1))
  # color_legend <- c("red", "darkkhaki", "green3", "blue", "pink", "black")
  # legend3d("topright", legend = c("1","2","3","4c","4s","NA"), pch = 16, col = color_legend, cex=1, inset=c(0.02))
}


# plot(reducedDims(sling_obj)$PCA, col = as.factor(slingshot_clusternames), pch=16, asp = 1)
# lines(SlingshotDataSet(sling_obj), lwd=2, col='black')

#sds <- SlingshotDataSet(sling_obj)
#rd <- as.data.frame(reducedDim(sds))

#getLineages(reducedDims(sling_obj)$PCA, slingshot_clusternames, start.clus= '1')

#slingshot_clusternames2 <-  meta2$Group[match(rownames(reducedDim(sce_obj)), rownames(meta2))]
#slingshot_clusternames2 <-  new_meta$Group[match(rownames(reducedDim(sce_obj)), rownames(new_meta))]

