args = commandArgs(TRUE)
samplename = toString(args[1])
vanloo_wedge_file = toString(args[2])
peifer_file = toString(args[3])

source("Parser.R")
source("Plotting.R")

suppressMessages(library(biovizBase))
suppressMessages(library(ggbio))
suppressMessages(library(GenomicRanges))

# samplename = "0c7af04b-e171-47c4-8be5-5db33f20148e"
# vanloo_wedge_infile = "data/vanloo_wedge/4_copy_number/0c7af04b-e171-47c4-8be5-5db33f20148e_segments.txt"
# #peifer_infile_subclonal = "data/peifer/Copy_Number/KICH_0c7af04b_subclonal_cn.txt"
# peifer_infile = "data/peifer/Copy_Number/UCEC_af96db5a_allelic_states.txt"

vanloo_wedge = parse.cn(vanloo_wedge_infile)
peifer = parse.cn.peifer(peifer_infile)

vector_of_names = c("vanloo_wedge", "peifer")
list_of_tables = list(vanloo_wedge, peifer)

# print("Loading BB")
# b = load.battenberg("../vanloo_wedge/train1/battenberg/ac757c6d-acc6-45f0-94e0-e6aa1f3709b9//ac757c6d-acc6-45f0-94e0-e6aa1f3709b9_subclones.txt")
# print("Loading Titan")
# t = load.titan("../broad/titan/ac757c6d-acc6-45f0-94e0-e6aa1f3709b9_cluster01//ac757c6d-acc6-45f0-94e0-e6aa1f3709b9_cluster01_segs.txt")
# print("Loading BAF")
# baf = load.baf("../vanloo_wedge/train1/battenberg/ac757c6d-acc6-45f0-94e0-e6aa1f3709b9/ac757c6d-acc6-45f0-94e0-e6aa1f3709b9.BAFsegmented.txt")
# #logr = load.baf("../vanloo_wedge/train1/battenberg/ac757c6d-acc6-45f0-94e0-e6aa1f3709b9/ac757c6d-acc6-45f0-94e0-e6aa1f3709b9.logRsegmented.txt")
# print("Loading LogR")
# start = system.time()
# logr = load.logr("../vanloo_wedge/train1/battenberg/ac757c6d-acc6-45f0-94e0-e6aa1f3709b9/ac757c6d-acc6-45f0-94e0-e6aa1f3709b9_mutantLogR.tab")
# print(system.time()-start)
# save(logr, filename="../vanloo_wedge/train1/battenberg/ac757c6d-acc6-45f0-94e0-e6aa1f3709b9/logr.RData")

# Load ideogram data
data(hg19IdeogramCyto)

for (chrom in unique(vanloo_wedge$chromosome)) {
  print(chrom)
  
  # Create the tracks
  p = getIdeogram(chrom)
  # p.baf = ggplot(baf[baf$Chromosome==chrom,]) + aes(x=Position, y=BAF) + geom_point(size=0.8, alpha=0.8) + ylim(0,1)
  # p.logr = ggplot(logr[logr$Chromosome==chrom,]) + aes(x=Position, y=LogR) + geom_point(size=0.8, alpha=0.8) + ylim(-2,4)
  vanloo_wedge_track = create.segs.track(vanloo_wedge[vanloo_wedge$chromosome==chrom,], start_colname="start", end_colname="end", max_colname="major_cn", min_colname="minor_cn")
  peifer_track = create.segs.track(peifer[peifer$Chromosome==chrom,], start_colname="Start", end_colname="End", max_colname="A", min_colname="B")
  
  # Plot
  png(filename=paste("4_copy_number/",samplename, "_chr", chrom, ".png", sep=""), width=1250, height=350) 
  print(tracks(p, vanloo_wedge=vanloo_wedge_track, peifer=peifer_track, heights=c(1.2, 2, 2)))
  dev.off()
}