#!/usr/bin/env Rscript
#Jenny Smith
#4/11/21

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(rmarkdown))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(snpReportR))

#Step 0.
#Prepare the vcf files  and the RNAseq counts in data/ directory.
#VCF Data must be saved with `{ sample_name }.vcf` or `{ sample_name }.vcf.gz`
#RNAseq data must have generic names degs_results and counts.


#Step 1.
#Configure gmailr package and pandoc
#https://cran.r-project.org/web/packages/gmailr/vignettes/gmailr.html
gmailr::gm_auth_configure() #this needs to be run at the top of the script to send emails
sender         <- 'my_email@gmail.com'
recipients     <- c('person1@gmail.com')

#pandoc path: Must provide path to pandoc since Rscript does not natively know where to find it, unlike running rmakrdown::render() interactively
Sys.setenv(RSTUDIO_PANDOC=rmarkdown::find_pandoc(cache = FALSE)$dir)
print(Sys.getenv("RSTUDIO_PANDOC"))


#Step 2.
#Define Common Input Parameters and files

  #4A.Create the genomic references Bioconductor objects
  #Default: "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gff3.gz"
  gen_refs <- suppressMessages(snpReportR::download_gen_refs())

  #4B. define the parameterized inputs as a list
  counts_var <- "{{{ counts_results }}}"
  if(exists(counts_var)){
    counts_results <- suppressMessages(get(counts_var))
  }else{
    counts_results <- read.delim(counts_var, sep="\t")
  }

  degs_var <- "{{{ degs_results }}}"
  if(exists(degs_var)){
    degs_results <- suppressMessages(get(degs_var))
  }else{
    degs_results <- read.delim(degs_var, sep="\t")
  }

  rnaseq <- list(counts=counts_results,
                 degs=degs_results)

  Params <- c(gen_refs, rnaseq)


#Step3.
# Create the Analysis reports "in batch" for the sample IDs provided
# Data should already be prepared according to  step 0 as an .rda object (need to workout the logic on how to sub in an rda object)
outfiles <- list()
sample_list <- unlist(strsplit("{{{ sample_names }}}",split = ","))

rendered_reports <- lapply(sample_list, function(sample_name){

  #3A. Read in the VCF data
  vcf_obj <- if (exists(sample_name)){
    suppressMessages(get(sample_name))
  }else{
    suppressMessages(get(load(paste0("data/",sample_name, ".rda"))))
  }

  #3B. update the params for this individual sample
  Params <- c(vcf_obj, Params)


  #3C. Define the output HTML report file
  filename <- paste0(sample_name, "_snpReport.html")
  outfile <- here("{{{ directory }}}", filename)
  outfiles[[sample_name]] <- outfile

  rmarkdown::render(input = here("inst/templates/Report_HTML_v3_JSmith.Rmd"),
                    output_format="html_document", #change to `all`, but having an error in latex pdf conversion
                    output_file=filename,
                    output_dir = "{{{ directory }}}",
                    params=Params)
  return(outfiles)

})



#Step 5.
#Send the automated email message

send_emails <- lapply(1:length(rendered_reports), function(i){

  sample_name <- suppressMessages(names(rendered_reports[[i]]))
  report <- suppressMessages(unlist(rendered_reports[[i]]))
  email_message <- blastula::render_email(here("inst/templates/Email_Body_v2_JSmith.Rmd"))
  email <- gmailr::gm_mime() %>%
        gmailr::gm_to(recipients) %>%
        gmailr::gm_from(sender) %>%
        gmailr::gm_subject("snpReportR") %>%
        gmailr::gm_html_body(email_message$html_html) %>%
        gmailr::gm_attach_file(report)
        # gmailr::gm_attach_file(here::here("snpReporter_logo.png"),id = "foobar")


      # gmailr::gm_create_draft(email)
      gmailr::gm_send_message(mail=email)
      print(paste0("Finished Sending Email for ", sample_name))
})


# Notes and Issues

#ISSUE: if using the pre-built dataset in snpReportR::setup(),
#inputting the `snpReportR::counts_results` and snpReporR::degs_results object results with the entire dataframe printed out to this R script.
#This makes the script too large (>5Mb) to open in Rstudio...
#Not sure how to address this, bc that actually makes the script more reproducible...
#but makes it much less editable and readable...
# so this works  setup(sample_names=sample_names, counts_results="counts_results",degs_results="degs_results", directory="my_reports")
# and this does not b/c filesize setup(sample_names=sample_names, counts_results=counts_results,degs_results=degs_results, directory="my_reports")


