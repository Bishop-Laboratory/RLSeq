#' Generate an RLSeq Report
#'
#' @param object An RLRanges object.
#' @param reportPath A path indicating the report output HTML file. Default: "rlreport.html"
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @param ... Arguments passed to `rmarkdown::render()`
#' @return TRUE
#' @examples 
#' \dontrun{
#' 
#' # Example dataset
#' rlbase <- "https://rlbase-data.s3.amazonaws.com"
#' pks <- file.path(rlbase, "peaks", "SRX1025890_hg38.broadPeak")
#' cvg <- file.path(rlbase, "coverage", "SRX1025890_hg38.bw")
#' 
#' # Run RLSeq
#' rlr <- RLRanges(pks, coverage = cvg, genome = "hg38", mode = "DRIP")
#' rlr <- RLSeq(rlr)
#' 
#' # Generate HTML report
#' report(rlr, reportPath="report.html")
#' 
#' }
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
report <- function(object,
    reportPath = "rlreport.html",
    quiet = FALSE,
    ...) {

    # Get the template
    template <- system.file("Rmd", "report.Rmd", package = "RLSeq")
    
    # Render template
    # note: We use callr::r here because if you don't, then you get
    # this issue: 
    # https://community.rstudio.com/t/rmarkdown-displays-a-plot-when-its-not-supposed-to/93757/2
    callr::r(
        function(template, object, reportPath, quiet, ...) { 
            rmarkdown::render(
                template,
                params = list(
                    "object" = object
                ),
                output_format = "html_document",
                output_dir = normalizePath(dirname(reportPath)),
                output_file = reportPath,
                quiet = quiet, 
                ...
            )
        }, 
        args  = list(template=template, 
                     object=object, 
                     reportPath=reportPath,
                     quiet=quiet, 
                     ...), 
        show = ! quiet
    )
    
    # Return
    return(TRUE)
}
