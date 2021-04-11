#!/usr/bin/env Rscript

library(jsonlite)
library(dplyr)
library(tidyr)
library(rmarkdown)
library(here)
library(snpReportR)


#Configure gmailr package.
#https://cran.r-project.org/web/packages/gmailr/vignettes/gmailr.html
gmailr::gm_auth_configure() #this needs to be run at the top of the script to send emails


#---------------------------------Define Parameters----------------------------------------------
vcf_files <- dir(here("data-raw"),
             pattern="applied.cancer.vcf",
             full.names = TRUE, recursive = T)

rnaseq_counts <- here("data-raw/edgeR_normcounts.tsv")
rnaseq_de <- here("data-raw/edgeR_test_results.tsv")
genome_ref <- NULL #set to different URL for a different genomic reference (default: gencode.v22.annotation.gff3.gz)

#Email Config
sender         <- 'jennyl.smith.workonly@gmail.com'
recipients     <- c('jennyl.smith.workonly@gmail.com')
smtp.username  <- 'jennyl.smith.workonly@gmail.com'


#Must provide path to pandoc since Rscript does not natively know where to find it,
#unlike running rmakrdown::render() interactively
Sys.setenv(RSTUDIO_PANDOC=rmarkdown::find_pandoc(cache = FALSE)$dir)

#---------------------------------Define Parameters----------------------------------------------

for (input_file in vcf_files){
  #Define the output HTML report file
  filename <- paste0(gsub("_applied.cancer.vcf|.wAnnot.vcf.gz","", basename(input_file)), "_snpReport.html")

  #Create the genomic references Bioconductor objects
  source("R/make_txDB.R")
  gen_refs <- download_gen_refs()

  #Define the input files
  Params <- list(data=input_file,
                 counts=rnaseq_counts,
                 degs=rnaseq_de)
  Params <- c(Params, gen_refs)

  #Render the Rmarkdown HTML report
  rmarkdown::render(input = here("R/Report_HTML_v2_JSmith.Rmd"),
                    output_format="all",
                    output_file=filename,
                    output_dir = here("Reports"),
                    params=Params)


  email_message <- blastula::render_email('R/Email_Body_v2_JSmith.Rmd')
  outfile <- here("Reports", filename)


  email <- gmailr::gm_mime() %>%
    gmailr::gm_to(recipients) %>%
    gmailr::gm_from(sender) %>%
    gmailr::gm_subject("snpReportR") %>%
    gmailr::gm_html_body(email_message$html_html)
    # gmailr::gm_attach_file(outfile)
    # gmailr::gm_attach_file(here::here("snpReporter_logo.png"),id = "foobar")


  gmailr::gm_create_draft(email)
  gmailr::gm_send_message(mail=email)
}


print("Finished Sending Email")


