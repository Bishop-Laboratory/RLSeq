#' Generate an RLSeq Report
#'
#' @param object An RLRanges object.
#' @param output A path indicating the report output HTML file. Default: "rlreport.html"
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @param ... Arguments passed to `rmarkdown::render()`
#' @return TRUE
#' @examples
#' 
#' rlbase <- "https://rlbase-data.s3.amazonaws.com"
#' pks <- file.path(rlbase, "peaks", "SRX1025890_hg38.broadPeak")
#' cvg <- file.path(rlbase, "coverage", "SRX1025890_hg38.bw")
#' rlr <- RLRanges(pks, coverage=cvg,  genome="hg38", mode="DRIP")
#' rlr <- RLSeq(rlr)
#' 
#' report(rlr)
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
report <- function(object, 
                   output = "rlreport.html", 
                   quiet=FALSE,
                   ...) {
  
  # Get the template
  template <- system.file("Rmd", "report.Rmd", package = "RLSeq")

  # Render template  
  rmarkdown::render(
    template,
    params = list(
      "object"=object
    ),
    output_format = "html_document",
    output_dir = normalizePath(dirname(output)),
    output_file = output
  )
  
  # Return
  return(TRUE)
}
