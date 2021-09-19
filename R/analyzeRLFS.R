#' Analyze RLFS
#'
#' Analyzes the enrichment of peaks within R-loop forming sequences.
#'
#' @param peaks A GRanges object containing the R-loop ranges to check.
#' @param genome UCSC genome which peaks were generated from. Ignored if "chrom_sizes" or "RLFS" is specified.
#' @param chrom_sizes A data.frame containing two columns with chromosome names
#' and chromosome sizes.
#' @param mask A GRanges object containing the regions of the genome to mask. Can be generated from
#' regioneR::getMask().
#' @param RLFS A GRanges object containing the R-loop forming sequences to
#' compare against.
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @param ... Arguments passed to `regioneR::permTest()`
#' @return A named list containing the results of RegioneR,
#' peaks annotated with the RLFS they overlap with,
#' and the Z score results within 3000 BP of an RLFS.
#' @examples
#'
#' result <- RLSeq::analyzeRLFS(RLSeq::SRX1025890_peaks, genome = "hg38")
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @importFrom utils capture.output
#' @export
analyzeRLFS <- function(peaks,
                        genome = "hg38",
                        chrom_sizes = NULL,
                        mask = NULL,
                        RLFS = NULL,
                        quiet = FALSE,
                        ...) {

  # TODO: Validate genome to avoid TRUE/FALSE Needed error from permTest

  # Check RLFS, chrom_sizes, and mask
  if (!quiet) {
    message("+ [i] Evaluating Inputs.")
  }
  if (is.null(genome) &
    (is.null(RLFS) | is.null(chrom_sizes) | is.null(mask))) {
    stop("Must provide genome UCSC org ID or chrom_sizes, mask, and RLFS")
  }
  if (is.null(RLFS)) {
    n_ <- capture.output(RLFS <- getRLFSAnno(genome))
  }
  if (is.null(chrom_sizes)) {
    n_ <- capture.output(chrom_sizes <- getChromSizes(genome))
  }
  if (is.null(mask)) {
    available_masks <- gsub(names(RLSeq::genomeMasks),
      pattern = "\\.masked",
      replacement = ""
    )
    if (!genome %in% available_masks) {
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
  genomeNow <- GenomicRanges::GRanges(chrom_sizes %>% dplyr::mutate(start = 1) %>%
    dplyr::select(
      seqnames = .data$X1,
      .data$start, end = .data$X2
    ))
  GenomeInfoDb::seqlevels(genomeNow) <- GenomeInfoDb::seqlevels(mask)
  GenomeInfoDb::seqinfo(genomeNow) <- GenomeInfoDb::seqinfo(mask)

  # Run RegioneR
  if (!quiet) {
    message("+ [ii] Running permTest.")
  }
  pt <- suppressWarnings(
    regioneR::permTest(
      A = peaks, B = RLFS,
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
    regioneR::localZScore(A = peaks, B = RLFS, pt, window = 5000, step = 50)
  )

  return(list(
    "perTestResults" = pt,
    "Z-scores" = z
  ))
}
