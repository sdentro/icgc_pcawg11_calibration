library(RColorBrewer)

####################################################################
# 2_mutation_assignments
####################################################################
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