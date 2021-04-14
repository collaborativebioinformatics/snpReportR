#' Set-up the project directories
#'
#' @param directory a character string for the name of the directory where reports will be saved.
#' @param sample_name sample_name is a character string for the sample identifier, which must be included in the VCF filename
#'
#' @return
#' @export
#'
#' @examples
#' setup("my_reports")
setup <- function(sample_name, directory="reports"){
  #https://sharla.party/post/usethis-for-reporting/
  #creates a directory if it doesn’t already exist - does not recognize capitaliztion differences
  usethis::use_directory(directory)
  usethis::use_template(template = "VCF_report_v2_rendering.R",
                        data = list(sample_name=sample_name,
                                    directory=directory),
                        package = "snpReportR",
                        open=TRUE)
}
