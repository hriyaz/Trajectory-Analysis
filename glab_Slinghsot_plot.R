glab_Slingshot_plot <- function(rd, clusters, curves, MST, PCX1, PCX2, text_size){

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
