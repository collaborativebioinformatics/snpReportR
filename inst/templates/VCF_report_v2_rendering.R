#!/usr/bin/env Rscript
#Jenny Smith
#4/11/21

library(dplyr)
library(tidyr)
library(rmarkdown)
library(here)

#Step 0.
#Prepare the vcf files  and the RNAseq counts in data/ directory.
#VCF Data must be saved with {{{ sample_name }}}.vcf or {{{ sample_name }}}.vcf.gz
#RNAseq data must have generic names degs and counts.

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
#Define Common Input Parameters and files

  #4A.Create the genomic references Bioconductor objects
  #Default: "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gff3.gz"
  gen_refs <- suppressMessages(snpReportR::download_gen_refs())

  #4B. define the parameterized inputs as a list
  rnaseq <- list(counts=if(class({{{ counts }}}) == "character"){
                          read.delim(counts, sep="\t")
                        }else{
                          {{{ counts }}}
                        },
                 degs=if(class({{{ degs }}})=="character"){
                          read.delim(degs, sep="\t")
                      }else{
                         {{{ degs }}}
                      })

  Params <- c(gen_refs, rnaseq)

#Step3.
# Create the Analysis reports "in batch" for the sample IDs provided
# Data should already be prepared according to  step 0 as an .rda object (need to workout the logic on how to sub in an rda object)
outfiles <- list()
render_reports <- lapply({{{ sample_names }}}, function(sample_name){

  #3A. Read in the VCF data
  vcf_obj <- if (exists(sample_name)){
    suppressMessages(get(sample_name))
  }else{
    suppressMessages(get(load(paste0(sample_name, ".rda"))))
  }

  #3B. update the params for this individual sample
  Params <- c(vcf_obj, Params)


  #3C. Define the output HTML report file
  filename <- paste0(sample_name, "_snpReport.html")
  outfile <- here({{{ directory }}}, filename)
  outfiles[[sample_name]] <- outfile

  rmarkdown::render(input = here("inst/templates/Report_HTML_v3_JSmith.Rmd"),
                    # output_format="html", #change to `all`, but having an error in latex pdf conversion
                    output_file=filename,
                    output_dir = {{{ directory }}},
                    params=Params)


})



#Step 5.
#Send the automated email message
#here('inst/rmarkdown/templates/email_generation/skeleton/skeleton.Rmd')
email_message <- blastula::render_email(here("Email_Body_v2_JSmith.Rmd"))

send_emails <- lapply(names(outfiles), function(sample_name){
      email <- gmailr::gm_mime() %>%
        gmailr::gm_to(recipients) %>%
        gmailr::gm_from(sender) %>%
        gmailr::gm_subject("snpReportR") %>%
        gmailr::gm_html_body(email_message$html_html)
        gmailr::gm_attach_file(outfiles[[sample_name]])
        # gmailr::gm_attach_file(here::here("snpReporter_logo.png"),id = "foobar")


      # gmailr::gm_create_draft(email)
      gmailr::gm_send_message(mail=email)
      print(paste0("Finished Sending Email for ", sample_name))
})






