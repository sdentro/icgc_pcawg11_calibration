####################################################################
# 1_purity_ploidy
####################################################################
library(ggplot2)
library(reshape2)
plot.purity.heatmap = function(purity, pngfilename, title) {
  purity.m = melt(data.matrix(purity))
  p = ggplot(purity.m) + 
    aes(Var2, Var1) + 
    geom_tile(aes(fill=value), colour="white") + 
    scale_fill_gradient(low="white", high="steelblue") +
    xlab("Method") + 
    ylab("Sample") +
    theme_bw() + 
    ggtitle(title)
  
  png(pngfilename, width=500, height=500)
  print(p)
  dev.off()
}

####################################################################
# 2_mutation_assignments
####################################################################
library(RColorBrewer)
#' Plot the heatmap associated with the identity array
plotHeatmapFull = function(identity.array, outfilename) {
  no.muts = nrow(identity.array)
  # Plot every 10th point to speed up plotting
  png(outfilename)
  heatmap(identity.array, Rowv=NA, Colv=NA, scale="none")
  dev.off()
}

#' Plot the heatmap associated with the identity array
plotHeatmapMedium = function(identity.array, outfilename) {
  no.muts = nrow(identity.array)
  # Plot every 10th point to speed up plotting
  png(outfilename)
  heatmap(identity.array[seq(10,no.muts,10),seq(10,no.muts,10)], Rowv=NA, Colv=NA, scale="none")
  dev.off()
}