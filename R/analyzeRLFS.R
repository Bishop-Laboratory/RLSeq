#' Analyze RLFS
#'
#' Analyzes the enrichment of peaks within R-loop forming sequences.
#'
#' @param object An RLRanges object.
#' @param mask GRanges object containing masked genomic ranges.
#'  Not needed unless masked genome unavailable (see RLSeq::genomeMasks).
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @param ... Arguments passed to `regioneR::permTest()`
#' @return A named list containing the results of RegioneR,
#' peaks annotated with the RLFS they overlap with,
#' and the Z score results within 3000 BP of an RLFS.
#' @examples
#' 
#' pks <- file.path(rlbase, "peaks", "SRX1025890_hg38.broadPeak")
#' rlr <- RLRanges(pks, genome="hg38", mode="DRIP")
#' rlr <- analyzeRLFS(rlr)
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @importFrom utils capture.output
#' @export
analyzeRLFS <- function(object,
                        mask = NULL,
                        quiet = FALSE,
                        ...) {

  # Check RLFS, chrom_sizes, and mask
  if (!quiet) {
    message("+ [i] Evaluating Inputs.")
  }
  genome <- GenomeInfoDb::genome(object)[1]
  RLFS <- getRLFSAnno(object)
  chrom_sizes <- getChromSizes(object)
  if (is.null(mask)) {
    available_masks <- gsub(names(RLSeq::genomeMasks),
      pattern = "\\.masked",
      replacement = ""
    )
    if (! genome %in% available_masks) {
      stop(genome, " is not available in mask list. You may generate it with ")
    } else {
      mask <- RLSeq::genomeMasks[[paste0(genome, ".masked")]]
    }
  }

  # Prevent stranded assignment
  GenomicRanges::strand(RLFS) <- "*"
  RLFS <- GenomicRanges::reduce(RLFS)

  # Final check
  stopifnot("data.frame" %in% class(chrom_sizes) &
    "GRanges" %in% class(mask) &
    "GRanges" %in% class(RLFS))

  # Run RLFS perm test
  genomeNow <- GenomicRanges::GRanges(
    dplyr::mutate(chrom_sizes, start = 1, end = .data$size)
  )
  GenomeInfoDb::seqlevels(genomeNow) <- GenomeInfoDb::seqlevels(mask)
  GenomeInfoDb::seqinfo(genomeNow) <- GenomeInfoDb::seqinfo(mask)

  # Run RegioneR
  if (!quiet) {
    message("+ [ii] Running permTest.")
  }
  pt <- suppressWarnings(
    regioneR::permTest(
      A = object, B = RLFS,
      genome = genomeNow,
      mask = mask,
      randomize.function = regioneR::circularRandomizeRegions,
      evaluate.function = regioneR::numOverlaps,
      alternative = "greater",
      ...
    )
  )

  # Return Z scores
  if (!quiet) {
    message("+ [iii] Extracting pileup.")
  }
  z <- suppressWarnings(
    regioneR::localZScore(A = object, B = RLFS, pt, window = 5000, step = 50)
  )
  
  # Add results to object
  slot(object@metadata$results, "rlfsRes") <- list(
    "perTestResults" = pt,
    "Z-scores" = z
  )

  return(object)
}
