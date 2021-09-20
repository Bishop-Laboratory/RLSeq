#' R-Loop region test
#'
#' Tests the overlap of peaks within R-loop regions
#'
#' @param object An RLRanges object with genome "hg38".
#' @return An RLRanges object with test results included.
#' @examples
#' pks <- file.path(rlbase, "peaks", "SRX1025890_hg38.broadPeak")
#' rlr <- RLRanges(pks, genome="hg38", mode="DRIP", quiet=TRUE)
#' 
#' rlr <- rlRegionTest(rlr)
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
rlRegionTest <- function(object) {

  # Wrangle the peaks
  pkName <- names(object)
  toTest <- object %>%
    tibble::as_tibble() %>%
    dplyr::mutate(seqnames = as.character(.data$seqnames),
                  name = {{ pkName }}) %>%
    dplyr::select(chrom = .data$seqnames, .data$start, .data$end, .data$name)
  
  # Get RLRegions 
  # TODO: NEEDS to be in RLHub 
  rlregions_table <- file.path(rlbase, "RLHub", "rlregions_table.rda")
  tmp <- tempfile()
  download.file(rlregions_table, destfile = tmp, quiet = TRUE)
  load(tmp)

  # Get the RL Regions
  rlReg <- tableToRegions(rlregions_table)

  # Get shared seqnames
  sharedSeqs <- intersect(toTest$chrom, rlReg$chrom)
  toTest <- dplyr::filter(toTest, .data$chrom %in% sharedSeqs)
  rlReg <- dplyr::filter(rlReg, .data$chrom %in% sharedSeqs)

  # Get the genome
  chromSizes <- getChromSizes(object)

  # Test on all annotations
  olap <- valr::bed_intersect(toTest, rlReg, suffix = c("__peaks", "__rlregion"))
  sig <- valr::bed_fisher(toTest, rlReg, genome = chromSizes)
  
  # Return to object
  slot(object@metadata$results, "rlRegionRes") <- list(
    "Overlap" = olap,
    "Test_results" = sig
  )

  # Return results
  return(object)
}
