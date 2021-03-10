#!/usr/bin/env Rscript 

library(jsonlite)
library(dplyr)
library(tidyr)
library(rmarkdown)
library(here)

here::here()

#Configure gmailr package. 
#https://cran.r-project.org/web/packages/gmailr/vignettes/gmailr.html
gmailr::gm_auth_configure() #this needs to be run at the top of the script to send emails
# gm_auth()

#---------------------------------Define Parameters----------------------------------------------
currentDate    <- Sys.Date()
files <- dir("Data", pattern="wAnnot.vcf.gz", full.names = TRUE, recursive = T)

#Email Config
sender         <- 'jennyl.smith.workonly@gmail.com'
recipients     <- c('jennyl.smith.workonly@gmail.com')
smtp.username  <- 'jennyl.smith.workonly@gmail.com'

Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")

#---------------------------------Define Parameters----------------------------------------------

for (input_file in files){
  filename <- paste0(gsub("_applied.cancer.vcf|.wAnnot.vcf.gz","", basename(input_file)), "_snpReport.html")
  Params <- list(data=here(input_file))

  rmarkdown::render(input = here:::here("R_Code/Report_v2_JSmith.Rmd"),
                    output_file=filename, 
                    output_dir = here("Reports"),
                    params=Params)
  
  email_message <- blastula::render_email('R_Code/Email_Body_JSmith.Rmd')
  outfile <- here("Reports", filename) 

  
  email <- gmailr::gm_mime() %>%
    gmailr::gm_to(recipients) %>%
    gmailr::gm_from(sender) %>%
    gmailr::gm_subject("snpReportR") %>%
    gmailr::gm_html_body(email_message$html_str) %>% 
    gmailr::gm_attach_file(outfile) %>% 
    gmailr::gm_attach_file(here::here("snpReporter_logo.png"), id = "foobar")
  
  
  # gmailr::gm_create_draft(email)
  gmailr::gm_send_message(mail=email)
}

  
print("Finished Sending Email")


