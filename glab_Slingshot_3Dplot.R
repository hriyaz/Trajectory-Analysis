glab_Slingshot_plot_3D <- function(sds, meta, meta_trait){
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
