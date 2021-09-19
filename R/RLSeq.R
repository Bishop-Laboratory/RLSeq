#' Run RLSeq and generate report
#'
#' Runs the full RLSeq pipeline and generates an analysis report in HTML format.
#'
#' @param peaks A GRanges object containing the R-loop ranges to check.
#' @param genome UCSC genome which peaks were generated from.
#' Only "hg38" and "mm10" currently available.
#' Use RLSeq::liftUtil() to convert to the correct format, if needed.
#' @param coverage (optional) The path to the bigWig file corresponding to `peaks`.
#' If supplied, correlation analysis will be performed.
#' @param outputFile A path indicating the report output HTML file.
#' @param dataOnly If TRUE, the HTML report will not be created.
#' @param sampleName The name to give to this sample in the report. Default: "User-supplied Sample".
#' @param mode The mode type this sample belongs to. See options at RLSeq::modes.
#' @param ... Arguments passed to `rmarkdown::render()`
#' @return A Named List with the datasets passed to RLSeq::makeReport().
#' @examples
#'
#' URL <- paste0(RLSeq:::RLBASE_BW_URL, "SRX1025890_hg38.bw")
#' RLSeq::RLSeq(RLSeq::SRX1025890_peaks,
#'   coverage = BW_FILE,
#'   genome = "hg38", outputFile = "report.html"
#' )
#' file.remove(BW_FILE)
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
RLSeq <- function(peaks, genome,
                  coverage = NULL,
                  sampleName = "User-supplied Sample",
                  mode = "User-supplied Sample",
                  outputFile = "report.html",
                  dataOnly = FALSE,
                  ...) {

  # Verify mode
  if (mode != "User-supplied Sample") {
    stopifnot(!mode %in% RLSeq::modes$mode)
  }

  # Get data sets
  message("[1] RLFS Perm Test...")
  rlfsRes <- analyzeRLFS(peaks = peaks, genome = genome)

  message("[2] Predict Condition...")
  pred <- predictCondition(rlfsRes = rlfsRes)

  message("[3] Feature Enrichment Test...")
  featTest <- featureEnrich(peaks = peaks, genome = genome)

  message("[4] Correlation Analysis...")
  if (!is.null(coverage)) {
    if (genome == "hg38") {
      corrRes <- corrAnalyze(coverage = coverage, genome = genome)
    } else {
      warning("Only hg38 is currently available for correlation analysis. Skipping...")
      corrRes <- NA
    }
  } else {
    message("No coverage provided... skipping.")
    corrRes <- NA
  }

  message("[5] Gene Annotation...")
  annoGenes <- geneAnnotation(peaks, genome = genome)

  message("[6] R-loop Region Analysis...")
  rlRegions <- rlRegionTest(peaks, genome = genome)

  message("[7] Collating results...")
  resLst <- list(
    "rlfsRes" = rlfsRes,
    "predictedCondition" = pred,
    "featureTest" = featTest,
    "corrRes" = corrRes,
    "annoGenes" = annoGenes,
    "RLoopRegions" = rlRegions,
    "sampleName" = sampleName,
    "mode" = mode,
    "genome" = genome
  )

  if (!dataOnly) {
    message("[8] Generating HTML Report...")
    report(resLst, outputFile = outputFile)
  }

  return(resLst)
}
