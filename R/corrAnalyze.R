#' Analyze Correlations
#'
#' Finds the pairwise correlation in signal around gold-standard R-Loop sites
#' between the query sample and the RMapDB database.
#'
#' @param object An RLRanges object.
#' @return A named list containing the results of correlation analysis.
#' @examples
#'
#' pks <- file.path(rlbase, "peaks", "SRX1025890_hg38.broadPeak")
#' cvg <- file.path(rlbase, "coverage", "SRX1025890_hg38.bw")
#' rlr <- RLRanges(pks, genome="hg38", mode="DRIP")
#' 
#' rlr <- RLSeq::corrAnalyze(rlr)
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
corrAnalyze <- function(object) {
  
  # Get genome
  genome <- GenomeInfoDb::genome(object)[1]
  stopifnot(genome == "hg38")
  
  # Get coverage
  coverage <- object@metadata$coverage
  if (coverage == "" || (! file.exists(coverage) && ! urlExists(coverage))) {
    stop("No coverage found. Content of 'coverage' slot: ", 
         coverage, ". Set coverage with object@metadata$coverage <-")
  }
  
  # TODO: This MUST go through RLHub
  gssig <- file.path(rlbase, "RLHub", "gsSignalRLBase.rda")
  tmp <- tempfile()
  download.file(gssig, destfile = tmp, quiet = TRUE)
  load(tmp)
  
  # Get the signal around GS R-loop sites
  bw <- getGSSignal(coverage, gssignal = gsSignalRLBase)

  # Wrangle BW into tibble
  bw <- bw %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::rename(chrom = .data$seqnames) %>%
    dplyr::mutate(chrom = as.character(.data$chrom))

  # Get positions
  positions <- gsSignalRLBase %>%
    dplyr::select(.data$location) %>%
    dplyr::mutate(
      chrom = gsub(.data$location,
        pattern = "(.+)_(.+)_(.+)",
        replacement = "\\1"
      ),
      start = gsub(.data$location,
        pattern = "(.+)_(.+)_(.+)",
        replacement = "\\2"
      ),
      end = gsub(.data$location,
        pattern = "(.+)_(.+)_(.+)",
        replacement = "\\3"
      )
    ) %>%
    dplyr::select(-.data$location) %>%
    dplyr::mutate(dplyr::across(c("start", "end"), as.integer))

  # Summarize across gs intervals with valr
  bwMap <- valr::bed_map(x = positions, y = bw, value = sum(score))

  # Combine with the original matrix
  combinedMat <- gsSignalRLBase %>%
    dplyr::inner_join(
      bwMap %>%
        dplyr::mutate(location = paste0(
          .data$chrom, "_",
          .data$start, "_",
          .data$end
        )) %>%
        dplyr::select(.data$location, a_ = .data$value),
      by = "location"
    ) %>%
    dplyr::distinct(.data$location, .keep_all = TRUE) %>%
    tibble::column_to_rownames("location") %>%
    as.matrix()
  
  # Rename column
  colnames(combinedMat)[
    colnames(combinedMat) == "a_"
  ] <- object@metadata$sampleName

  # Find the correlation
  corMat <- cor(combinedMat)
  
  # Add back to object
  slot(object@metadata$results, "correlationMat") <- corMat

  return(object)
}
