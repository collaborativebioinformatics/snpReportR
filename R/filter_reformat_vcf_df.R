#Jenny Smith
#Jan. 8, 2021
#Purpose: reformat the VCF output from v2.5 CTAT Mutations Pipeline

#' Filter and Reformat the VCF from v2.5 CTAT Mutations Pipeline
#'
#' @param vcf.df is a dataframe derived from vcfR package
#' @param vcf.s4 is a s4 vectors object from bioconductor VariantAnnotation package
#'
#' @return
#' @export
#'
#' @examples
#'\dontrun{
#' vcf <- vcfR::read.vcfR("/path/to/vcf")
#' vcf.df <- cbind(as.data.frame(vcfR::getFIX(vcf)), vcfR::INFO2df(vcf))
#' vcf.s4 <- VariantAnnotation::readVcf("/path/to/vcf")
#' vcf.filtered <- filter_ctat_vcf(vcf.df,vcf.s4)
#'}
#'
#' @import dplyr
filter_ctat_vcf <- function(vcf.df, vcf.s4){

  #Define annotations in VCF header lines
  functional_annots_names <- VariantAnnotation::header(vcf.s4) %>%
    VariantAnnotation::info(.) %>%
    as.data.frame(.)

  functional_annots_names <- functional_annots_names["ANN","Description"] %>%
    gsub("^.+\\'(.+)\\'","\\1",.) %>%
    str_split(., pattern = "\\|") %>%
    unlist() %>%
    gsub("\\s", "", .)

  functional_annots_names <- functional_annots_names[-length(functional_annots_names)]


  # expland the annotations from ANN attribute of the VCF file for the first 3 transcripts.
  #  https://www.biostars.org/p/226965/
  functional_annots.df <- data.frame(do.call(rbind, strsplit(as.vector(vcf.df$ANN), split = "\\|")))
  functional_annots.df <- functional_annots.df[,1:45]  #keep only the first 3 transcripts
  colnames(functional_annots.df) <- paste(functional_annots_names,rep(1:3, each=15), sep="_")

  #Run the filtering function
  vcf.df.filter <- vcf.df %>%
    mutate(S4_Vector_IDs=names(SummarizedExperiment::rowRanges(vcf.s4))) %>%
    bind_cols(., functional_annots.df) %>%

    mutate(rsID=ifelse(!is.na(RS), paste0("rs", RS), RS)) %>%
    mutate_at(vars(chasmplus_pval,vest_pval), ~as.numeric(.)) %>%
    group_by(GENE) %>%
    mutate(Number_SNVs_per_Gene=n()) %>%
    ungroup() %>%

    dplyr::select(GENE,Number_SNVs_per_Gene, COSMIC_ID,
                  rsID,CHROM:ALT,
                  FATHMM,SPLICEADJ,
                  matches("chasmplus_(pval|score)"),
                  matches("vest_(pval|score)"),
                  TISSUE,TUMOR,
                  Annotation_1,Annotation_Impact_1,Feature_Type_1,
                  Transcript_BioType_1,
                  coding_DNA_change_1=HGVS.c_1,
                  protein_change_1=HGVS.p_1,
                  -ANN, everything(), ANN) %>%
    dplyr::filter(grepl("PATHOGENIC", FATHMM) | !is.na(SPLICEADJ)) %>%
    dplyr::filter(grepl("HIGH|MODERATE",Annotation_Impact_1) | !is.na(SPLICEADJ)) %>%
    arrange(desc(chasmplus_score), desc(vest_score),
            desc(Number_SNVs_per_Gene), Annotation_Impact_1)


  return(vcf.df.filter)

}

