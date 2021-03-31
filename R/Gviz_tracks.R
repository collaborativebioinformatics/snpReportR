#Jenny Smith
#1/8/20
#Purpose: Create Gviz tracks for SNVs of interest.

gene_tracks <- function(vcf_s4, txdb) {

  #ensure seqlevels style is UCSC (eg chr1)
  # seqlevelsStyle(gff) <- "UCSC"
  seqlevelsStyle(txdb) <- "UCSC"

  #define ranges to visualize
  vcf.ranges <- rowRanges(vcf_s4)
  vcf.ranges <- vcf.ranges[order(vcf.ranges)]
  chrs <- as.vector(seqnames(vcf.ranges)@values)
  print(chrs)

  #define gene models
  genesGR <- genes(txdb)
  # IDmap <-  as.data.frame(gff)
  genes.sub <- suppressWarnings(subsetByOverlaps(genesGR, vcf.ranges))
  geneModel <- select(txdb,
                      keys=names(genes.sub),
                      columns=c("EXONCHROM","EXONSTART","EXONEND","EXONSTRAND" ,
                                "EXONNAME","GENEID","TXNAME","TXTYPE"),
                      keytype = "GENEID") %>%
    # left_join(., IDmap, by=c("GENEID"="gene_id")) %>%
    dplyr::select(chromosome=EXONCHROM,
                  start=EXONSTART,
                  end=EXONEND,
                  strand=EXONSTRAND,
                  gene=GENEID,
                  exon=EXONNAME,
                  transcript=TXNAME)
                  # symbol=geneSymbol) #will need to use a GFF merged in

  for (chr in chrs){
    #Subset each gene of interest by chromosome
    ranges.sub <- vcf.ranges[seqnames(vcf.ranges)==chr]
    ranges.sub$ID <- names(ranges.sub)

    #Define tracks
    itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
    gtrack <- GenomeAxisTrack()
    indeltrack <- AnnotationTrack(ranges.sub,
                                  genome = "hg38",
                                  name="SNVs",
                                  group=ranges.sub$ID,
                                  id=ranges.sub$ID,
                                  showId=TRUE,
                                  cex.id=1.0,
                                  fontsize=10)
    grtrack <- GeneRegionTrack(geneModel,
                               genome = "hg38",
                               chromosome = chr,
                               name = "Gene Model",
                               transcriptAnnotation = "transcript",
                               # collapseTranscripts="longest",
                               showId=TRUE)

    plotTracks(c(itrack,gtrack,indeltrack,grtrack), collapse=FALSE)
  }
}
