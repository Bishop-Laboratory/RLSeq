#' Analyze RLFS
#'
#' Analyzes the enrichment of ranges within R-loop forming sequences (RLFS).
#'
#' @param object An RLRanges object.
#' @param mask GRanges object containing masked genomic ranges.
#'  Not needed unless masked genome unavailable (see genomeMasks).
#'  Custom masks can be generated using regioneR::getMask().
#' @param quiet If TRUE, messages are suppressed. Default: FALSE.
#' @param ... Arguments passed to `regioneR::permTest()`
#' @return An RLRanges object with RLFS analysis results included.
#' @examples
#'
#' # Example dataset
#' rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Perform RLFS analysis
#' rlr <- analyzeRLFS(rlr)
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
#' @export
analyzeRLFS <- function(object,
    mask = NULL,
    quiet = FALSE,
    ...) {
    
    if (!quiet) message(" - Evaluating Inputs...")

    # Check RLFS, chrom_sizes, and mask
    genome <- GenomeInfoDb::genome(object)[1]
    RLFS <- getRLFSAnno(object)
    chrom_sizes <- getChromSizes(object)
    if (is.null(mask)) {
        available_masks <- gsub(names(genomeMasks),
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
    }

    # Prevent stranded assignment
    GenomicRanges::strand(RLFS) <- "*"
    RLFS <- GenomicRanges::reduce(RLFS)

    # Final check
    stopifnot(
        methods::is(chrom_sizes, "data.frame") &
            methods::is(mask, "GRanges") &
            methods::is(RLFS, "GRanges")
    )

    # Run RLFS perm test
    genomeNow <- GenomicRanges::GRanges(
        dplyr::mutate(chrom_sizes, start = 1, end = .data$size)
    )
    GenomeInfoDb::seqlevels(genomeNow) <- GenomeInfoDb::seqlevels(mask)
    GenomeInfoDb::seqinfo(genomeNow) <- GenomeInfoDb::seqinfo(mask)

    if (!quiet) message(" - Running permTest...")
    
    # Run RegioneR
    gr <- GenomicRanges::GRanges(object)
    pt <- suppressWarnings(
        regioneR::permTest(
            A = gr, B = RLFS,
            genome = genomeNow,
            mask = mask,
            randomize.function = regioneR::circularRandomizeRegions,
            evaluate.function = regioneR::numOverlaps,
            alternative = "greater",
            ...
        )
    )
    
    if (!quiet) message(" - Extracting pileup...")

    # Return Z scores
    z <- suppressWarnings(
        regioneR::localZScore(A = gr, B = RLFS, pt, window = 5000, step = 50)
    )

    # Add results to object
    methods::slot(object@metadata$results, "rlfsRes") <- list(
        "perTestResults" = pt,
        "Z-scores" = z
    )

    if (!quiet) message(" - Done.")

    return(object)
}
