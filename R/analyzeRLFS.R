#' Analyze RLFS
#'
#' Analyzes the enrichment of ranges within R-loop forming sequences (RLFS).
#' See *details*.
#'
#' @param object An [RLRanges] object.
#' @param mask GRanges object containing masked genomic ranges.
#'  Not needed unless masked genome unavailable (see [genomeMasks]).
#'  Custom masks can be generated using [regioneR::getMask].
#' @param useMask If FALSE, masked genome is not used. This is
#' not recommended unless a mask is unavailable as it can lead to spurious
#' results. Default: TRUE.
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @param noZ If TRUE, Z-score distribution is not calculated. Default: FALSE.
#' @param ntimes Number of permutations to perform (default: 100).
#' @param stepsize The step size for calculating the Z score distribution.
#' Default: 50. See also [regioneR::localZScore].
#' @param ... Arguments passed to [regioneR::permTest].
#' @details
#'
#' R-loop forming sequences are regions of the genome with sequences that are
#' favorable for R-loop formation. They are computationally
#' predicted with the
#' [QmRLFS-finder](https://github.com/piroonj/QmRLFS-finder)
#' software program and serve as a data-independent test of whether a sample
#' has mapped R-loops robustly or not.
#'
#' ## Method
#'
#' Permutation testing is implemented via [regioneR::permTest] such that,
#' for each permutation, R-loop peaks were randomized using
#' [regioneR::circularRandomizeRegions] and then the number of overlaps with
#' RLFS are counted. 100 permutations are used by default to build an empirical
#' distribution for peak/RLFS overlap. Then the true number of overlaps from
#' non-randomized peaks and RLFS are compared to the null distribution to
#' calculate Z-score and significance of enrichment. Finally, a Z-score
#' distribution was calculated (using [regioneR::localZScore])
#' 5kb upstream and downstream of the average RLFS midpoint.
#'
#' These results are subsequently used in the binary classification of
#' the sample as "POS" (maps R-loops) or "NEG" (does not map R-loops). See also
#' [predictCondition].
#'
#' @return An RLRanges object with RLFS analysis results
#' accessible via `RLSeq::rlresult(object, "rlfsRes")`. Contains the
#' following structure:
#'
#' - `perTestResults`
#'   * An object of the class `permTestResultsList` from `regioneR` with
#'   the results of permutation testing. See also
#'   [regioneR::permTest] for full description.
#' - `Z-scores`
#'   * An object of the class `localZScoreResultsList` from `regioneR`.
#'   Contains the results of local Z-score analysis +/-5kb around each RLFS.
#'   See also [regioneR::localZScore].
#'
#' @examples
#'
#' # Example dataset
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Perform RLFS analysis (remove ntimes=2 and noZ=TRUE for a typical analysis)
#' rlr <- analyzeRLFS(rlr, ntimes = 2, noZ = TRUE)
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
            min.parallel = Inf,
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
            min.parallel = Inf,
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
