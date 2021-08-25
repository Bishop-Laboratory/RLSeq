#' R-Loop region test
#'
#' Tests the overlap of peaks within R-loop regions
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
    as.data.frame() %>%
    tibble::rownames_to_column(var = "name") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(seqnames = as.character(seqnames)) %>%
    dplyr::select(chrom=seqnames, start, end, name)
  
  # Get the RL Regions
  rlReg <- RSeqR::rlRegions %>%
    dplyr::mutate(
      chrom = as.character(gsub(Location, pattern = "(.+):(.+)\\-(.+)", replacement = "\\1")),
      start = as.numeric(gsub(Location, pattern = "(.+):(.+)\\-(.+)", replacement = "\\2")),
      end = as.numeric(gsub(Location, pattern = "(.+):(.+)\\-(.+)", replacement = "\\3"))
    ) %>%
    dplyr::select(
      chrom, start, end, name=`RL Region`
    ) 
  
  # Get shared seqnames
  sharedSeqs <- intersect(toTest$chrom, rlReg$chrom)
  toTest <- dplyr::filter(toTest, chrom %in% sharedSeqs)
  rlReg <- dplyr::filter(rlReg, chrom %in% sharedSeqs)
  
  # Get the genome
  chromSizes <- RSeqR:::getChromSizes(genome) %>%
    dplyr::rename(chrom = X1, size = X2) 
  
  # Test on all annotations
  olap <- valr::bed_intersect(toTest, rlReg, suffix = c("__peaks", "__RLFS"))
  sig <- valr::bed_projection(toTest, rlReg, genome=chromSizes)
  
  # Return results
  return(
    list(
      "Overlapp" = olap,
      "Test_resutls" = sig
    )
  )
}
