#!/usr/bin/env Rscript 

library(jsonlite)
library(dplyr)
library(tidyr)
library(rmarkdown)
library(here)



setwd(here())

#Configure gmailr package. 
#https://cran.r-project.org/web/packages/gmailr/vignettes/gmailr.html
# Sys.getenv("GMAILR_APP")
gmailr::gm_auth_configure() #this needs to be run at the top of the script to send emails
# gm_auth()

#---------------------------------Define Parameters----------------------------------------------
currentDate    <- Sys.Date()
files <- dir("Data", pattern="wAnnot.vcf.gz", full.names = TRUE, recursive = T)

#Email Config
sender         <- 'jennyl.smith.workonly@gmail.com'
recipients     <- c('jennyl.smith.workonly@gmail.com')
smtp.username  <- 'jennyl.smith.workonly@gmail.com'

#---------------------------------Define Parameters----------------------------------------------

for (input_file in files){
  filename <- paste0(gsub("_applied.cancer.vcf|.wAnnot.vcf.gz","", basename(input_file)), "_snpReport.html")
  Params <- list(data=here(input_file))
  Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
  rmarkdown::render(input = here("R_Code/Report_v2_JSmith.Rmd"),
                    output_file=filename, 
                    output_dir = here("Reports"),
                    params=Params)
  
  outfile <- here(file.path("Reports", filename))
  message <- paste(readLines(outfile, warn = FALSE), collapse = "\n")
  
  email <- gmailr::gm_mime() %>%
    gmailr::gm_to(recipients) %>%
    gmailr::gm_from(sender) %>%
    gmailr::gm_subject("snpReportR") %>%
    gmailr::gm_html_body(message) 
  
  gmailr::gm_create_draft(email)
  gmailr::gm_send_message(mail=email)
}

  
print("Finished Sending Email")


### Testing 

#https://cran.r-project.org/web/packages/gmailr/vignettes/sending_messages.html
#This works!
# email <- gm_mime() %>%
#   gm_to(recipients) %>%
#   gm_from(sender) %>%
#   gm_subject("A Test Email") %>%
#   gm_html_body(
#     '<h1>A plot of <b>MotorTrend</b> data <i>(1974)</i></h1>
#     <br><img src="cid:foobar">') 
# gm_create_draft(email)
# gm_send_message(mail=email)
