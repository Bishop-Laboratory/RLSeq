#' Analyze Correlations
#'
#' Finds the pairwise correlation in signal around gold-standard R-Loop sites 
#' between the query sample and the RMapDB database.
#'
#' @param coverage The path to a coverage (.bigWig/.bw) file (can be a URL).
#' @param genome The UCSC genome ID to use. (Currently only "hg38" is supported)
#' @return A named list containing the results of correlation analysis.
#' @examples
#'
#' BW_URL <- paste0("https://rmapdb-data.s3.us-east-2.amazonaws.com/bigwigs/",
#'                  "rseq-coverage-unstranded/SRX1025890_hg38.bw")
#' result <- RLSeq::corrAnalyze(BW_URL, genome="hg38")
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
corrAnalyze <- function(coverage, genome="hg38") {
  
  # Get the signal around GS R-loop sites
  bw <- getGSSignal(coverage, genome=genome)
  
  # Wrangle BW into tibble
  bw <- bw %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::rename(chrom = .data$seqnames) %>%
    dplyr::mutate(chrom = as.character(.data$chrom))
  
  # Get positions
  positions <- RLSeq::gsSignalRMapDB %>%
    dplyr::select(.data$location) %>%
    dplyr::mutate(chrom = gsub(.data$location, pattern = "(.+)_(.+)_(.+)",
                               replacement = "\\1"),
                  start = gsub(.data$location, pattern = "(.+)_(.+)_(.+)",
                               replacement = "\\2"),
                  end = gsub(.data$location, pattern = "(.+)_(.+)_(.+)",
                             replacement = "\\3")) %>%
    dplyr::select(-.data$location) %>%
    dplyr::mutate(dplyr::across(c("start", "end"), as.integer))
  
  # Summarize across gs intervals with valr
  bwMap <- valr::bed_map(x = positions, y = bw, value = sum(score))
  
  # Combine with the original matrix
  combinedMat <- gsSignalRMapDB %>%
    dplyr::inner_join(
      bwMap %>%
        dplyr::mutate(location = paste0(.data$chrom, "_",
                                        .data$start, "_", 
                                        .data$end)) %>%
        dplyr::select(.data$location, user_supplied=.data$value),
      by="location"
    ) %>%
    dplyr::distinct(.data$location, .keep_all = TRUE) %>%
    tibble::column_to_rownames("location") %>%
    as.matrix() 
  
  # Find the correlation
  corMat <- cor(combinedMat)
  
  return(corMat)
  
}
