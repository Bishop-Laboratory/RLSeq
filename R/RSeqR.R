#' Run RSeqR and generate report
#'
#' Runs the full RSeqR pipeline and generates an analysis report in HTML format.
#'
#' @param peaks A GRanges object containing the R-loop ranges to check.
#' @param genome UCSC genome which peaks were generated from. 
#' Only "hg38" and "mm10" currently available. 
#' Use RSeqR::liftUtil() to convert to the correct format, if needed.
#' @param coverage (optional) The path to the bigWig file corresponding to `peaks`. 
#' If supplied, correlation analysis will be performed.
#' @param outputFile A path indicating the report output HTML file. 
#' @param ... Arguments passed to `rmarkdown::render()`
#' @return A Named List with the datasets passed to RSeqR::makeReport().
#' @examples
#' 
#' URL <- paste0("https://rmapdb-data.s3.us-east-2.amazonaws.com/bigwigs/",
#'               "rseq-coverage-unstranded/SRX1025890_TC32_NT_DRIP.hg38.bw")
#' BW_FILE <- "SRX1025890.bw"
#' download.file(URL, destfile=BW_FILE)
#' RSeqR::RSeqR(RSeqR::SRX1025890_peaks, coverage=BW_FILE,
#'              genome="hg38", outputFile = "report.html")
#' file.remove(BW_FILE)
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
RSeqR <- function(peaks, genome, coverage=NULL, outputFile = "RSeqR_Report.html", ...) {
  
  message("[1] RLFS Perm Test...")
  rlfsRes <- analyzeRLFS(peaks=peaks, genome=genome)
  
  message("[2] Predict Condition...")
  pred <- predictCondition(rlfsRes = rlfsRes)
  
  message("[3] Feature Enrichment Test...")
  featTest <- featureEnrich(peaks=peaks, genome=genome)
  
  message("[4] Correlation Analysis...")
  if (! is.null(coverage)) {
    if (genome == "hg38") {
      corrRes <- corrAnalyze(coverage=coverage, genome = genome)
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
  rlRegions <- rlRegionTest(peaks, genome=genome)
  
  message("[7] Make Report...")
  resLst <- list(
    "rlfsRes" = rlfsRes,
    "predictedCondition" = pred,
    "featureTest" = featTest,
    "corrRes" = corrRes,
    "annoGenes" = annoGenes,
    "RLoopRegions" = rlRegions
  )
  makeReport(resLst, outputFile = outputFile)
  
  return(resLst)
}

