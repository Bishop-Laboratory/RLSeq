#' R-Loop region test
#'
#' Tests the overlap of user-supplied ranges with R-loop regions.
#'
#' @param object An RLRanges object with genome "hg38".
#' @return An RLRanges object with test results included.
#' @examples 
#' \dontrun{
#' 
#' # Get example data
#' rlbase <- "https://rlbase-data.s3.amazonaws.com"
#' pks <- file.path(rlbase, "peaks", "SRX1025890_hg38.broadPeak")
#' cvg <- file.path(rlbase, "coverage", "SRX1025890_hg38.bw")
#' rlr <- RLRanges(pks, coverage = cvg, genome = "hg38", mode = "DRIP")
#' 
#' # RL Region Test
#' rlRegionTest(rlr)
#' 
#' }
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
rlRegionTest <- function(object) {
    
    # Stop if not genome == hg38
    if (! GenomeInfoDb::genome(object)[1] == "hg38") {
        stop("Only 'hg38' ranges are allowed. Please convert to 'hg38' ",
             "using a lift-over or skip this step.")
    }

    # Wrangle the peaks
    pkName <- names(object)
    toTest <- object %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
            seqnames = as.character(.data$seqnames),
            name = {{ pkName }}
        ) %>%
        dplyr::select(chrom = .data$seqnames, .data$start, .data$end, .data$name)

    # Get RLRegions
    # TODO: NEEDS to be in RLHub
    rlbase <- "https://rlbase-data.s3.amazonaws.com"
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
