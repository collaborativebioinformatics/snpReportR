#!/usr/bin/env Rscript
#Jenny Smith
#4/11/21

library(dplyr)
library(tidyr)
library(rmarkdown)
library(here)

#Step 0.
#Prepare the vcf files in data-raw directory
#Data must be saved with {{{ sample_name }}}.vcf or {{{ sample_name }}}.vcf.gz


#Step 1.
#Configure gmailr package and pandoc
#https://cran.r-project.org/web/packages/gmailr/vignettes/gmailr.html
gmailr::gm_auth_configure() #this needs to be run at the top of the script to send emails
sender         <- 'jennyl.smith.workonly@gmail.com'
recipients     <- c('jennyl.smith.workonly@gmail.com')

#pandoc path
#Must provide path to pandoc since Rscript does not natively know where to find it,
#unlike running rmakrdown::render() interactively
Sys.setenv(RSTUDIO_PANDOC=rmarkdown::find_pandoc(cache = FALSE)$dir)


#Step 2.
#Define Input Parameters and files
  #4A. Read in the VCF data
  vcf_obj <- {{{ sample_name }}}

  #4B.Create the genomic references Bioconductor objects
  #Default: "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gff3.gz"
  gen_refs <- suppressMessages(snpReportR::download_gen_refs())

  #4C. define the parameterized inputs as a list
  rnaseq <- list(counts=read.delim(here("data-raw/edgeR_normcounts.tsv"),sep="\t"),
                 degs=read.delim(here("data-raw/edgeR_test_results.tsv"), sep="\t"))

  Params <- vcf_obj
  Params <- c(Params, gen_refs, rnaseq)


  #4D. Define the output HTML report file
  filename <- paste0({{{ sample_name }}}, "_snpReport.html")



#Step 4.
#Render the Rmarkdown HTML report
# here("inst/rmarkdown/templates/expressed_variant_reports/skeleton/skeleton.Rmd")
rmarkdown::render(input = here("Report_HTML_v3_JSmith.Rmd"),
                    output_format="html", #change to `all`, but having an error in latex pdf conversion
                    output_file=filename,
                    output_dir = {{{ directory }}},
                    params=Params)


#Step 5.
#Send the automated email message
#here('inst/rmarkdown/templates/email_generation/skeleton/skeleton.Rmd')
email_message <- blastula::render_email(here("Email_Body_v2_JSmith.Rmd"))
outfile <- here({{{ directory }}}, filename)


email <- gmailr::gm_mime() %>%
    gmailr::gm_to(recipients) %>%
    gmailr::gm_from(sender) %>%
    gmailr::gm_subject("snpReportR") %>%
    gmailr::gm_html_body(email_message$html_html)
    # gmailr::gm_attach_file(outfile)
    # gmailr::gm_attach_file(here::here("snpReporter_logo.png"),id = "foobar")


  # gmailr::gm_create_draft(email)
  gmailr::gm_send_message(mail=email)


print("Finished Sending Email")


