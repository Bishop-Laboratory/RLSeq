#' Run RLSeq
#'
#' Runs the full RLSeq pipeline.
#'
#' @param object An RLRanges object.
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @param ... Arguments passed to analyzeRLFS.
#' @return An RLRanges object with all results available.
#' @details 
#' 
#' The \code{RLSeq()} function does all of the following:
#' 
#' \enumerate{
#'   \item RLFS Perm Test. Runs the \code{analyzeRLFS()} function to test the 
#'   enrichment of user-supplied ranges within R-loop-forming sequences. 
#'   \item Predict Condition. Runs the \code{predictCondition()} function to 
#'   predict whether the user-supplied sample robustly maps R-loops or not.
#'   \item Feature enrichment test. Runs the \code{featureEnrich()} function to
#'   test the enrichment of user-supplied ranges within R-loop-relevant genomic features.
#'   \item Correlation Analysis. Runs the \code{corrAnalyze()} function to test
#'   the correlation of user-supplied R-loop signal with other samples in RLBase around
#'   "gold-standard" R-loop regions. 
#'   \item Gene annotation. Runs the \code{geneAnnotation()} function to find the overlap
#'   of genes with the user-supplied ranges. 
#'   \item R-loop Region Analysis. Runs the \code{rlRegionTest()} function to find
#'   the overlap of user-supplied ranges with consensus R-loop sites (RL-Regions).
#' }
#' 
#' @examples
#'
#' # Example RLRanges
#' rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Run RLSeq
#' rlr <- RLSeq(rlr)
#' 
#' @export
RLSeq <- function(object, quiet = FALSE, ...) {
    if (!quiet) message("[1/6] RLFS Perm Test")
    object <- analyzeRLFS(object, quiet = TRUE, ...)

    if (!quiet) message("[2/6] Predict Condition")
    object <- predictCondition(object)

    if (!quiet) message("[3/6] Feature Enrichment Test")
    if (GenomeInfoDb::genome(object)[1] %in% c("hg38", "mm10")) {
        object <- featureEnrich(object, quiet = TRUE)
    } else {
        warning(
            "RLSeq only contains built-in annotations for 'hg38' and 'mm10'",
            " genomes. Please lift-over or run featureEnrich() separately ",
            "with custom annotations. Skipping."
        )
        Sys.sleep(2)
    }


    if (!quiet) message("[4/6] Correlation Analysis")
    if (object@metadata$coverage != "") {
        if (GenomeInfoDb::genome(object)[1] == "hg38") {
            if (.Platform$OS.type != "windows") {
                object <- corrAnalyze(object)
            } else {
                if (! quiet) {
                    warning("corrAnalyze does not work on windows OS. Please",
                            " run corrAnalyze() directly with `force=TRUE` to override.")
                }
            }
        } else {
            warning(
                "Only 'hg38' genome ranges are available",
                " for correlation analysis. Skipping."
            )
            Sys.sleep(2)
        }
    } else {
        message("No coverage provided... skipping.")
        corrRes <- NA
    }

    if (!quiet) message("[5/6] Gene Annotation")
    objectanno <- try(geneAnnotation(object, quiet = TRUE), silent = TRUE)
    if ("try-error" %in% class(objectanno)) {
        warning(objectanno)
    } else {
        object <- objectanno
    }

    if (!quiet) message("[6/6] R-loop Region Analysis")
    if (GenomeInfoDb::genome(object)[1] == "hg38") {
        object <- rlRegionTest(object)
    } else {
        warning(
            "Only 'hg38' genome ranges",
            " are available for RL-region analysis. Skipping."
        )
        Sys.sleep(1)
    }

    return(object)
}
