#' Generate report
#'
#' Generates RSeqR Quality Report
#'
#' @param data A named list containing `annotated_peaks`, `feature_overlaps`,
#' and `rlfs_results`
#' @param outputFile A path indicating the report output HTML file. 
#' @param ... Arguments passed to `rmarkdown::render()`
#' @return NULL
#' @examples
#' 
#' RSeqR::makeReport(RSeqR::SRX1025890, output_file = "report.html")
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
makeReport <- function(data, output_file = "RSeqReport.html", ...) {
  template <- system.file("Rmd", "reportTemplate.Rmd", package = "RSeqR")
  rmarkdown::render(template, 
                    params = data, 
                    output_format = "html_document",
                    output_dir = normalizePath(dirname(output_file)),
                    output_file = output_file, ...)
  return(NULL)
}
