#' Analyze RLFS
#'
#' Analyzes the enrichment of ranges within R-loop forming sequences (RLFS).
#'
#' @param object An RLRanges object.
#' @param mask GRanges object containing masked genomic ranges.
#'  Not needed unless masked genome unavailable (see \code{RLSeq:::genomeMasks}).
#'  Custom masks can be generated using regioneR::getMask().
#' @param useMask If FALSE, masked genome is not used. This is
#' not recommended unless a mask is unavailable as it can lead to spurious
#' results. Default: TRUE.
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @param noZ If TRUE, Z-score distribution is not calculated. Default: FALSE.
#' @param ntimes Number of permutations to perform (default: 100).
#' @param stepsize The step size for calculating the Z score distribution.
#' Default: 50. See also [regioneR::localZScore()].
#' @param ... Arguments passed to [regioneR::permTest()].
#' @return An RLRanges object with RLFS analysis results included.
#' @examples
#'
#' # Example dataset
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Perform RLFS analysis (remove ntimes=2 and noZ=TRUE for a typical analysis)
#' rlr <- analyzeRLFS(rlr, ntimes = 2, noZ = TRUE)
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
#' @export
analyzeRLFS <- function(object,
    mask = NULL,
    quiet = FALSE,
    useMask = TRUE,
    noZ = FALSE,
    ntimes = 100,
    stepsize = 50,
    ...) {
    if (!quiet) message(" - Evaluating Inputs...")

    # Check RLFS, chrom_sizes, and mask
    genome <- GenomeInfoDb::genome(object)[1]
    RLFS <- getRLFSAnno(object)
    chrom_sizes <- getChromSizes(object)
    if (is.null(mask) & useMask) {
        available_masks <- gsub(
            names(genomeMasks),
            pattern = "\\.masked",
            replacement = ""
        )
        if (!genome %in% available_masks) {
            stop(
                genome, " is not available in mask list. See 'mask' param."
            )
        } else {
            mask <- genomeMasks[[paste0(genome, ".masked")]]
        }
    } else if (useMask) {
        stopifnot(methods::is(mask, "GRanges"))
    }

    # Prevent stranded assignment
    GenomicRanges::strand(RLFS) <- "*"
    RLFS <- GenomicRanges::reduce(RLFS)

    # Final check
    stopifnot(
        methods::is(chrom_sizes, "data.frame") & methods::is(RLFS, "GRanges")
    )

    # Run RLFS perm test
    genomeNow <- GenomicRanges::GRanges(
        dplyr::mutate(chrom_sizes, start = 1, end = .data$size)
    )

    if (!is.null(mask) & useMask) {
        GenomeInfoDb::seqlevels(
            genomeNow,
            pruning.mode = "coarse"
        ) <- GenomeInfoDb::seqlevels(mask)
        GenomeInfoDb::seqinfo(genomeNow) <- GenomeInfoDb::seqinfo(mask)
    }

    if (!quiet) message(" - Running permTest...")

    # Run RegioneR
    gr <- GenomicRanges::GRanges(object)
    gr <- GenomicRanges::trim(gr)
    if (useMask) {
        pt <- regioneR::permTest(
            A = gr, B = RLFS,
            genome = genomeNow,
            mask = mask,
            randomize.function = regioneR::circularRandomizeRegions,
            evaluate.function = regioneR::numOverlaps,
            alternative = "greater",
            ntimes = ntimes,
            ...
        )
    } else {
        pt <- regioneR::permTest(
            A = gr, B = RLFS,
            genome = genomeNow,
            randomize.function = regioneR::circularRandomizeRegions,
            evaluate.function = regioneR::numOverlaps,
            alternative = "greater",
            ntimes = ntimes,
            ...
        )
    }


    if (!quiet) message(" - Extracting pileup...")

    if (!noZ) {
        # Return Z scores
        z <- regioneR::localZScore(
            A = gr, B = RLFS, pt, window = 5000, step = stepsize
        )
    } else {
        z <- NA
    }
    # Add results to object
    methods::slot(object@metadata$results, "rlfsRes") <- list(
        "perTestResults" = pt,
        "Z-scores" = z
    )

    if (!quiet) message(" - Done.")

    return(object)
}
