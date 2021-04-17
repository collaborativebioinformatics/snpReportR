## code to prepare `DATASET` dataset goes here
# https://r-pkgs.org/data.html
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(VariantAnnotation))
library(tibble)
library(snpReportR)


#How can I parameterize this so that its not a `for loop`?
sample_names <- c()
for (input_file in dir("data-raw", pattern = "applied.cancer.vcf", full.names = TRUE)){
  sample_name <- paste0(gsub("_applied.cancer.vcf(.gz)?|init.wAnnot.vcf(.gz)?","", basename(input_file)))
  sample_names <- c(sample_names, sample_name)

  #As a dataframe
  vcf <-  suppressMessages(vcfR::read.vcfR(input_file))
  vcf_df <- cbind(as.data.frame(vcfR::getFIX(vcf)),
                  vcfR::INFO2df(vcf)) %>%
    tibble::as_tibble(.,  .name_repair = "unique")

  #as a GRanges/S4 bioconductor compatible object
  #Easiest would be to transform this into a dataframe rather than use the two different types of objects (df vs S4)
  vcf_s4 <- suppressMessages(VariantAnnotation::readVcf(input_file))

  #Filter the dataset for relevant mutations
  filtered_df <- snpReportR::filter_ctat_vcf(vcf.df=vcf_df, vcf.s4=vcf_s4)

  #Create a list of the input fileformats to be used in the report generation
  input_data <- list(data=filtered_df,
                     vcf.df=vcf_df,
                     vcf.s4=vcf_s4)

  # https://stackoverflow.com/questions/49673667/how-to-use-devtoolsuse-data-on-a-list-of-data-frames
  assign(sample_name, value=input_data)
  do.call("use_data", list(as.name(sample_name), overwrite = TRUE))
}

#Save vector of sample IDs/ sample names
use_data(sample_names, overwrite = TRUE)

#Save RNAseq Data
counts <- read.delim("data-raw/edgeR_normcounts.tsv")
use_data(counts,overwrite = TRUE)

degs <- read.delim("data-raw/edgeR_test_results.tsv")
use_data(degs, overwrite = TRUE)
