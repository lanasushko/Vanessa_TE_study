setwd('/tmp/global2/ssushko/Vanessa_proj/4_mitochondrial')

library(Gviz)

# LR999940.1:1,717,400-1,736,114
chr <- "LR999940.1"
gen <- "ilVanCard2.2"

TEs=read.table('/tmp/global2/ssushko/Vanessa_proj/4_mitochondrial/TEs_in_region.bed')
colnames(TEs)=c('chrom','chromStart','chromEnd','feature','score','strand')
range=makeGRangesFromDataFrame(TEs)

options(ucscChromosomeNames=FALSE)
atrack <- AnnotationTrack(range, name = "TEs") # add annotation track

plotTracks(atrack)

gtrack <- GenomeAxisTrack() # add scale/coordinates track
plotTracks(list(gtrack, atrack))


# import BAM

# /ebio/abt6_projects7/small_projects/ssushko/Vanessa_proj/finished_analyses/selfalign2.2_mito_discovery/aln2.2_sorted.bam


library(Rsamtools)

afrom <- 1717400
ato <- 1736114

bamfile='/tmp/global2/ssushko/Vanessa_proj/4_mitochondrial/aln_pacbiovcardreads_sorted.bam'
alignmentTrack <- AlignmentsTrack(bamfile, 
                                  chromosome = chr, 
                                  from = afrom, 
                                  to = ato,
                                  isPaired = FALSE)  # Set FALSE for long reads

indexBam(bamfile)

plotTracks(alignmentTrack, 
           from = afrom, 
           to = ato, 
           chromosome = chr, 
           type = "coverage")




