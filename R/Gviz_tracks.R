#Jenny Smith
#1/8/20
#Purpose: Create Gviz tracks for SNVs of interest.

#' Title
#'
#' @param vcf_s4 is a s4 vectors object from bioconductor VariantAnnotation package
#' @param txdb is a transcript database from
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' txdb <- GenomicFeatures::makeTxDbFromGFF("gencode.v22.chr_patch_hapl_scaff.annotation.gtf.gz")
#' AnnotationDbi::saveDb(txdb, "gencode.v22.chr_patch_hapl_scaff.annotation.sqlite")
#' vcf.s4 <- VariantAnnotation::readVcf("/path/to/vcf")
#' gene_tracks(vcf_s4 = vcf.s4[1:3], txdb=txdb)
#' }
#'
#' @import dplyr
gene_tracks <- function(vcf_s4, txdb, genome="hg38") {

  #ensure seqlevels style is UCSC (eg chr1)
  # GenomeInfoDb::seqlevelsStyle(gff) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(txdb) <- "UCSC"

  #define ranges to visualize
  vcf.ranges <- SummarizedExperiment::rowRanges(vcf_s4)
  vcf.ranges <- vcf.ranges[order(vcf.ranges)]
  chrs <- as.vector(GenomeInfoDb::seqnames(vcf.ranges)@values)
  # print(chrs)

  #define gene models
  genesGR <- GenomicFeatures::genes(txdb)
  # IDmap <-  as.data.frame(gff)
  genes.sub <- suppressWarnings(IRanges::subsetByOverlaps(genesGR, vcf.ranges))
  geneModel <- AnnotationDbi::select(txdb,
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
    ranges.sub <- vcf.ranges[GenomeInfoDb::seqnames(vcf.ranges)==chr]
    ranges.sub$ID <- names(ranges.sub)

    #Define tracks
    itrack <- Gviz::IdeogramTrack(genome = genome, chromosome = chr)
    gtrack <- Gviz::GenomeAxisTrack()
    indeltrack <- Gviz::AnnotationTrack(ranges.sub,
                                  genome = genome,
                                  name="SNVs",
                                  group=ranges.sub$ID,
                                  id=ranges.sub$ID,
                                  showId=TRUE,
                                  cex.id=1.0,
                                  fontsize=10)
    grtrack <- Gviz::GeneRegionTrack(geneModel,
                               genome = genome,
                               chromosome = chr,
                               name = "Gene Model",
                               transcriptAnnotation = "transcript",
                               # collapseTranscripts="longest",
                               showId=TRUE)

    Gviz::plotTracks(c(itrack,gtrack,indeltrack,grtrack), collapse=FALSE)
  }
}
