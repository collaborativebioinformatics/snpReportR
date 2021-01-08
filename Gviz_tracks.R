#Jenny Smith 
#1/8/20
#Purpose: Create Gviz tracks for SNVs of interest. 

gene_tracks <- function(vcf_s4, transriptsGR) {
  
  vcf.ranges <- rowRanges(vcf_s4)
  vcf.ranges <- vcf.ranges[order(vcf.ranges)]
  chrs <- as.vector(seqnames(vcf.ranges)@values)
  # chrs <- chrs[grep("chr1$", chrs)]
  print(chrs)
  
  for (chr in chrs){
    ranges.sub <- vcf.ranges[seqnames(vcf.ranges)==chr]
    transcripts.sub <- suppressWarnings(subsetByOverlaps(transriptsGR, ranges.sub))
    
    itrack.1 <- IdeogramTrack(genome = "hg38", chromosome = chr)
    gtrack <- GenomeAxisTrack()
    grtrack <- GeneRegionTrack(transcripts.sub, 
                               genome = "hg38",
                               chromosome = "chr1", 
                               name = "Gene Model",
                               showId=TRUE)
    track <- AnnotationTrack(ranges.sub,
                             genome = "hg38",
                             name="SNVs",
                             ucscChromosomeNames=FALSE)
    
    plotTracks(c(itrack.1,gtrack,track,grtrack), collapse=FALSE)
  }
  
}