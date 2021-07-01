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
#' RSeqR::makeReport(RSeqR::SRX1025890, outputFile = "report.html")
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
makeReport <- function(data, outputFile = "RSeqReport.html", ...) {
  template <- system.file("Rmd", "reportTemplate.Rmd", package = "RSeqR")
  data <- data[names(data) %in% c(
    "corr_data", "anno_data", "rlfs_data", "bam_stats", 
    "read_qc_data", "configlist", "total_peaks"
  )]
  rmarkdown::render(template, 
                    params = data, 
                    output_format = "html_document",
                    output_dir = normalizePath(dirname(outputFile)),
                    output_file = outputFile, ...)
  return(NULL)
}
