#Jenny Smith
#3/30/21

#' Create TxDB and Genes GRanges object from a GFF file
#'
#' @param URL The URL of the reference genes in gff3 format
#' @param save boolean - whether to save the txDB as sqlite and the
#'
#' @return
#' @export
#'
#' @examples
#' objects <- download_gen_refs()
#'
#' head(objects$gtf)
#' objects$txdb
#'
download_gen_refs <- function(URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gff3.gz",
                              save=FALSE){

  #Define input and output filenames
  filename <- basename(URL)
  outname <- gsub(".gz","", filename)
  rds <- here::here(paste0(outname,".RDS"))
  sqlite <- here::here(paste0(outname,".sqlite"))

  #If the files already exist in the current working director
  if(file.exists(rds) & file.exists(sqlite)){
    message(paste("Reference files exist. Using ",sqlite, "and", rds, "."))
    refs <- list("gtf"=readRDS(rds),
                 "txdb"=AnnotationDbi::loadDb(sqlite))

    return(refs)
  }

  #Download the gff file to the current working directory
  gtf_file <- utils::download.file(url=URL,
                            destfile=filename)

  #The "phase" metadata column contains non-NA values for features of type stop_codon. This information was ignored.
  txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGFF(filename))
  gtf <- suppressWarnings(rtracklayer::import.gff3(filename))

  #optionally save the large data files
  if(save){
    base::saveRDS(gtf, rds)
    AnnotationDbi::saveDb(txdb, sqlite)
  }

  refs <- list("gtf"=gtf,
               "txdb"=txdb)
  return(refs)
}

