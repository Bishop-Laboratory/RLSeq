#' R-Loop region test
#'
#' Tests the overlap of user-supplied ranges with R-loop regions (RL regions).
#'
#' @param object An RLRanges object with genome "hg38".
#' @details
#'
#' R-loop regions (RL regions) are consensus sites of R-loop formation. For
#' more information, see [RLHub::rlregions]. The `rlRegionTest` is a simple
#' function which finds the overlap of user-supplied samples with RL regions and
#' calculates Fisher's exact test via [valr::bed_fisher].
#'
#' @return An RLRanges object with test results accessible via
#' `rlresult(object, "rlRegionRes")`.
#'
#' ### Structure
#'
#' The structure of the results is a named `list` containing the following:
#'
#' * `Overlap`
#'   - A `tbl` showing the overlap between RL regions and user-supplied ranges.
#'   - Column description:
#'     * `chrom` - The chromosome name
#'     * `start__peaks` - The starting position of the user-supplied peak in
#'     the overlap.
#'     * `end__peaks` - Same as above for end position.
#'     * `name__peaks` - The name of the user-supplied peak in the overlap
#'     (from `names(object)`).
#'     * `start/end/name__rlregion` - Same as above for RL regions.
#'     * `strand__rlregion` - The genomic strand of the RL region in the overlap.
#'     * `.overlap` - The size of the overlap.
#' * `Test_results`
#'   - A `tbl` showing the results of the Fisher's exact test.
#'   See [valr::bed_fisher].
#'
#' @examples
#'
#' # Example RLRanges data
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # RL Region Test
#' rlRegionTest(rlr)
#' @export
rlRegionTest <- function(object) {

    # Stop if not genome == hg38
    if (!GenomeInfoDb::genome(object)[1] == "hg38") {
        stop(
            "Only 'hg38' ranges are allowed. Please convert to 'hg38' ",
            "using a lift-over or skip this step."
        )
    }

    # Wrangle the peaks
    pkName <- names(object)
    toTest <- object %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(
            seqnames = as.character(.data$seqnames),
            name = {{ pkName }}
        ) %>%
        dplyr::select(
            chrom = .data$seqnames,
            .data$start,
            .data$end,
            .data$name
        )

    # Get RLRegions
    rlregions_table <- RLHub::rlregions_meta(quiet = TRUE)

    # Get the RL Regions
    rlReg <- tableToRegions(rlregions_table)

    # Get shared seqnames
    sharedSeqs <- intersect(toTest$chrom, rlReg$chrom)
    toTest <- dplyr::filter(toTest, .data$chrom %in% sharedSeqs)
    rlReg <- dplyr::filter(rlReg, .data$chrom %in% sharedSeqs)

    # Get the genome
    chromSizes <- getChromSizes(object)

    # Test on all annotations
    olap <- valr::bed_intersect(
        toTest,
        rlReg,
        suffix = c("__peaks", "__rlregion")
    )
    sig <- valr::bed_fisher(toTest, rlReg, genome = chromSizes)

    # Return to object
    methods::slot(object@metadata$results, "rlRegionRes") <- list(
        "Overlap" = olap,
        "Test_results" = sig
    )

    # Return results
    return(object)
}
