#' Set-up the project directories
#'
#' @param sample_names a character vector of the sample names/ sample IDs of the VCF files.
#' @param counts a character vector of the variable name or a filepath with the tab separated RNAseq counts.
#' @param degs a character vector of the variable name or a filepath with the tab separated results of EdgeR DE analysis.
#' @param directory a character string for the name of the directory where reports will be saved.
#'
#' @return
#' @export
#'
#' @examples
#' setup(sample_names=c("DRR131561_dx.variants.HC_hard_cutoffs",counts=counts,degs=degs, directory="my_reports")
setup <- function(sample_names, counts, degs, directory="reports"){
  #https://sharla.party/post/usethis-for-reporting/
  #creates a directory if it doesn’t already exist - does not recognize capitaliztion differences
  usethis::use_directory(directory)

  # #its opening the VCF_report_v2_rendering.R not in the reports directory? why?
  usethis::use_template(template = "VCF_report_v2_rendering.R",
                        data = list(sample_names=sample_names,
                                    counts=counts,
                                    degs=degs,
                                    directory=directory),
                        package = "snpReportR",
                        open=TRUE)
}
