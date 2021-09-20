#' Run RLSeq
#'
#' A convenience function which runs the full RLSeq pipeline.
#'
#' @param object An RLRanges object.
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @return An RLRanges object with all results available.
#' @examples
#' 
#' rlbase <- "https://rlbase-data.s3.amazonaws.com"
#' pks <- file.path(rlbase, "peaks", "SRX1025890_hg38.broadPeak")
#' cvg <- file.path(rlbase, "coverage", "SRX1025890_hg38.bw")
#' rlr <- RLRanges(pks, coverage=cvg,  genome="hg38", mode="DRIP")
#' rlr <- RLSeq(rlr)
#' 
#' @export
RLSeq <- function(object, quiet=FALSE) {

  if (! quiet) message("[1/6] RLFS Perm Test")
  object <- analyzeRLFS(object, quiet = TRUE)

  if (! quiet) message("[2/6] Predict Condition")
  object <- predictCondition(object)

  if (! quiet) message("[3/6] Feature Enrichment Test")
  object <- featureEnrich(object, quiet = TRUE)

  if (! quiet) message("[4/6] Correlation Analysis")
  if (object@metadata$coverage != "") {
    if (GenomeInfoDb::genome(object)[1] == "hg38") {
      object <- corrAnalyze(object)
    } else {
      warning("Only 'hg38' is available for correlation analysis. Skipping.")
      Sys.sleep(1)
    }
  } else {
    message("No coverage provided... skipping.")
    corrRes <- NA
  }

  if (! quiet) message("[5/6] Gene Annotation")
  object <- geneAnnotation(object, quiet = TRUE)

  if (! quiet) message("[6/6] R-loop Region Analysis")
  object <- rlRegionTest(object)
  
  return(object)
}
