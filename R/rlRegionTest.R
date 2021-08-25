#' R-Loop region test
#'
#' Tests the enrichment of peaks within R-loop regions
#'
#' @param peaks A GRanges object containing R-loop peaks
#' @param genome UCSC genome identifier to use. Can be only "hg38" currently. 
#' If you have "hg19" peaks, please use RSeqR::liftUtil() to convert them.
#' @return A named list containing the results of testing.
#' @examples
#' 
#' RSeqR::rlRegionTest(RSeqR::SRX1025890_peaks, genome="hg38")
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
rlRegionTest <- function(peaks, genom="hg38") {
  
  # Wrangle the peaks
  toTest <- peaks %>%
    tibble::as_tibble() %>%
    dplyr::select(seqnames, start, end)
  
  # Get the RL Regions
  
  
  # Get the genome
  chromSizes <- getChromSizes(genome) %>%
    dplyr::rename(chrom = X1, size = X2) 
  
  # Test on all annotations
  
  
  return(annoRes)
  
}



