#' Set-up the project analysis directories
#'
#' @param sample_names a character vector of the sample names/ sample IDs of the VCF files.
#' @param counts a character vector of the variable name of the data.frame object OR a filepath with the tab separated RNAseq counts.
#' @param degs a character vector of the variable name of the data.frame object OR a filepath with the tab separated results of EdgeR DE analysis.
#' @param directory a character string for the name of the directory where reports will be saved.
#'
#' @return
#' @export
#'
#' @examples
#' library(snpReportR)
#' setup(sample_names=sample_names, counts_results="counts_results",degs_results="degs_results", directory="my_reports")
setup <- function(sample_names, counts_results, degs_results, directory="reports"){
  #https://sharla.party/post/usethis-for-reporting/
  #creates a directory if it doesn’t already exist - does not recognize capitaliztion differences
  usethis::use_directory(directory)
  outdir <-  paste(Sys.Date(), "VCF_report_rendering.R",sep="_") %>%
    paste(directory, ., sep="/")

  #its opening the VCF_report_v2_rendering.R not in the reports directory? why?
  usethis::use_template(template = "VCF_report_rendering.R",
                        data = list(sample_names=lapply(sample_names, as.character),
                                    counts_results=counts_results,
                                    degs_results=degs_results,
                                    directory=directory),
                        package = "snpReportR",
                        save_as = outdir,
                        open=TRUE)
}
