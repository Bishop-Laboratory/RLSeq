#' Analyze RLFS
#'
#' Analyzes the enrichment of ranges within R-loop forming sequences (RLFS).
#'
#' @param object An RLRanges object.
#' @param mask GRanges object containing masked genomic ranges.
#'  Not needed unless masked genome unavailable (see RLSeq::genomeMasks).
#'  Custom masks can be generated using regioneR::getMask().
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @param ... Arguments passed to `regioneR::permTest()`
#' @return An RLRanges object with RLFS analysis results included.
#' @examples 
#' \dontrun{
#' 
#' # Example dataset
#' rlbase <- "https://rlbase-data.s3.amazonaws.com"
#' pks <- file.path(rlbase, "peaks", "SRX1025890_hg38.broadPeak")
#' cvg <- file.path(rlbase, "coverage", "SRX1025890_hg38.bw")
#' 
#' # Get RLRanges object
#' rlr <- RLRanges(pks, coverage = cvg, genome = "hg38", mode = "DRIP")
#' 
#' # Perform RLFS analysis
#' rlr <- analyzeRLFS(rlr)
#' 
#' }
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
        message(" - Evaluating Inputs...")
    }
    genome <- GenomeInfoDb::genome(object)[1]
    RLFS <- getRLFSAnno(object)
    chrom_sizes <- getChromSizes(object)
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
    genomeNow <- GenomicRanges::GRanges(
        dplyr::mutate(chrom_sizes, start = 1, end = .data$size)
    )
    GenomeInfoDb::seqlevels(genomeNow) <- GenomeInfoDb::seqlevels(mask)
    GenomeInfoDb::seqinfo(genomeNow) <- GenomeInfoDb::seqinfo(mask)

    # Run RegioneR
    if (!quiet) {
        message(" - Running permTest...")
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
        message(" - Extracting pileup...")
    }
    z <- suppressWarnings(
        regioneR::localZScore(A = object, B = RLFS, pt, window = 5000, step = 50)
    )

    # Add results to object
    slot(object@metadata$results, "rlfsRes") <- list(
        "perTestResults" = pt,
        "Z-scores" = z
    )

    if (!quiet) {
        message(" - Done.")
    }

    return(object)
}
