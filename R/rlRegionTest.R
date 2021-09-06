#' R-Loop region test
#'
#' Tests the overlap of peaks within R-loop regions
#'
#' @param peaks A GRanges object containing R-loop peaks
#' @param genome UCSC genome identifier to use. Can be only "hg38" currently. 
#' If you have "hg19" peaks, please use RLSeq::liftUtil() to convert them.
#' @return A named list containing the results of testing.
#' @examples
#' 
#' RLSeq::rlRegionTest(RLSeq::SRX1025890_peaks, genome="hg38")
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
rlRegionTest <- function(peaks, genome="hg38") {
  
  # Wrangle the peaks
  toTest <- peaks %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "name") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(seqnames = as.character(.data$seqnames)) %>%
    dplyr::select(chrom = .data$seqnames, .data$start, .data$end, .data$name)
  
  # Get the RL Regions
  rlReg <- RLSeq::rlRegions %>%
    dplyr::mutate(
      chrom = as.character(gsub(.data$Location, 
                                pattern = "(.+):(.+)\\-(.+)",
                                replacement = "\\1")),
      start = as.numeric(gsub(.data$Location,
                              pattern = "(.+):(.+)\\-(.+)",
                              replacement = "\\2")),
      end = as.numeric(gsub(.data$Location,
                            pattern = "(.+):(.+)\\-(.+)",
                            replacement = "\\3"))
    ) %>%
    dplyr::select(
      .data$chrom, .data$start, .data$end, name=.data$`RL Region`
    ) 
  
  # Get shared seqnames
  sharedSeqs <- intersect(toTest$chrom, rlReg$chrom)
  toTest <- dplyr::filter(toTest, .data$chrom %in% sharedSeqs)
  rlReg <- dplyr::filter(rlReg, .data$chrom %in% sharedSeqs)
  
  # Get the genome
  chromSizes <- getChromSizes(genome) %>%
    dplyr::rename(chrom = .data$X1, size = .data$X2) 
  
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
