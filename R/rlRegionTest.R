#' R-Loop region test
#'
#' Tests the overlap of user-supplied ranges with R-loop regions (R-loop consensus
#' sites derived from meta-analysis of RLBase)
#'
#' @param object An RLRanges object with genome "hg38".
#' @return An RLRanges object with test results included.
#' @examples
#'
#' # Example RLRanges data
#' rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))
#'
#' # RL Region Test
#' rlRegionTest(rlr)
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
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
