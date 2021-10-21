#' RLSeq Report
#'
#' Builds an HTML report to showcase the results available in the supplied
#' [RLRanges] object (see also [RLResults]). 
#' 
#' @param object An RLRanges object.
#' @param reportPath A path indicating the report output HTML file.
#'  Default: "rlreport.html"
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @param ... Arguments passed to `rmarkdown::render()`
#' @return TRUE
#' @examples
#'
#' # Example data with RLSeq() already run.
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Get a TMP file (only for example usae)
#' tmp <- tempfile()
#'
#' # Generate the report
#' report(rlr, reportPath = tmp)
#' @export
report <- function(object,
    reportPath = "rlreport.html",
    quiet = FALSE,
    ...) {

    # Check for missing packages and stop if found
    pkgs <- vapply(
        X = c("kableExtra", "DT", "rmarkdown"),
        FUN = requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
    if (any(!pkgs)) {
        pks <- paste0(names(pkgs)[which(!pkgs)], collapse = '", "')
        stop(
            'Packages needed for report() but not installed: "',
            pks, '"'
        )
    }

    # Get the template
    template <- system.file("Rmd", "report.Rmd", package = "RLSeq")

    # Render template
    # note: We use callr::r here because if you don't, then you get
    # this issue:
    # https://community.rstudio.com/t/
    # rmarkdown-displays-a-plot-when-its-not-supposed-to/93757/2
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
        args = list(
            template = template,
            object = object,
            reportPath = reportPath,
            quiet = quiet,
            ...
        ),
        show = !quiet
    )

    # Return
    return(TRUE)
}
