#' Analyze Correlations
#'
#' Finds the pairwise correlation in signal around gold-standard R-Loop sites 
#' between the query sample and the RMapDB database.
#'
#' @param coverage The results list from running analyzeRLFS().
#' @param genome The UCSC genome ID to use. (Currently only "hg38" is supported)
#' @return A named list containing the results of correlation analysis.
#' @examples
#'
#' URL <- paste0("https://rmapdb-data.s3.us-east-2.amazonaws.com/bigwigs/",
#'               "rseq-coverage-unstranded/SRX1025890_TC32_NT_DRIP.hg38.bw")
#' BW_FILE <- "SRX1025890.bw"
#' download.file(URL, destfile=BW_FILE)
#' result <- RSeqR::corrAnalyze(BW_FILE, genome="hg38")
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
corrAnalyze <- function(coverage, genome="hg38") {
  
  # Get the locations of the gs sites
  positions <- rownames(RSeqR::gsSignalRMapDB) 
  positions <- tibble::tibble(location = positions) %>%
    dplyr::mutate(seqnames = gsub(.data$location, pattern = "(.+)_(.+)_(.+)",
                                  replacement = "\\1"),
                  start = gsub(.data$location, pattern = "(.+)_(.+)_(.+)",
                               replacement = "\\2"),
                  end = gsub(.data$location, pattern = "(.+)_(.+)_(.+)",
                             replacement = "\\3")) %>%
    dplyr::select(-.data$location) %>%
    GenomicRanges::makeGRangesFromDataFrame()
  
  # Read in the bigWig file using these locations
  bw <- rtracklayer::import(con = rtracklayer::BigWigFile(coverage), 
                            selection = positions)
  
  # Wrangle BW into tibble
  bw <- bw %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::rename(chrom = .data$seqnames) %>%
    dplyr::mutate(chrom = as.character(.data$chrom))
  
  # Wrangle positions into tibble
  posTbl <- positions %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::rename(chrom = .data$seqnames) %>%
    dplyr::mutate(chrom = as.character(.data$chrom))
  
  # Summarize across gs intervals with valr
  bwMap <- valr::bed_map(x = posTbl, y = bw, value = sum(score))
  
  # Combine with the original matrix
  combinedMat <- gsSignalRMapDB %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "location") %>%
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
