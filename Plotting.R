####################################################################
# 1_purity_ploidy
####################################################################

plot.purity.heatmap = function(purity, pngfilename, title) {
  library(ggplot2)
  library(reshape2)
  purity.m = melt(purity, id.vars=c("sample"))
  purity.m$sample = as.factor(purity.m$sample)
  purity.m$variable = as.factor(purity.m$variable)
  purity.m$value = as.numeric(purity.m$value)
  
  p = ggplot(purity.m) + 
    aes(variable, sample) +
    geom_tile(aes(fill=value), colour="white") + 
    scale_fill_gradient(low="white", high="steelblue") +
    xlab("Method") + 
    ylab("Sample") +
    theme_bw() + 
    ggtitle(title)
  
  png(pngfilename, width=750, height=1000)
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

####################################################################
# 4_copy_number
####################################################################
getIdeogram = function(chromosome) {
  return(plotIdeogram(hg19IdeogramCyto, subchr=paste("chr",chromosome, sep=""), color="black", alpha=0.0))
}

create.segs.track = function(data, start_colname, end_colname, max_colname, min_colname) {
  if (nrow(data) > 0) {
    p = ggplot(data) + 
      geom_segment(mapping=aes_string(x=start_colname, xend=end_colname, y=max_colname, yend=max_colname), colour="purple", size=4) + 
      geom_segment(mapping=aes_string(x=start_colname, xend=end_colname, y=min_colname, yend=min_colname), colour="darkblue", size=4) + 
      ylim(0, 8)
  } else {
    p = ggplot(data) + ylim(0, 8)
  }
  return(p)
}
