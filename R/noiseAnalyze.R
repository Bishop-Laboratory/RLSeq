#' Analyze sample noise
#'
#' Analyzes the noiseness of the supplied sample using the method described by
#' *Diaz et al.*. See *details*.
#'
#' Currently, this does not work on windows.
#'
#' @param object An RLRanges object.
#' @param windows Genomics windows to use for quantifying signal.
#' Will be automatically supplied if not provided. It is recommended NOT
#' to specify this option for most analysis types, as doing so will impair
#' the ability to compare to RLBase samples. Default: NULL.
#' @param force Force `noiseAnalyze` to run, even if on Windows. Default: FALSE.
#' @return An RLRanges object with noise analysis results included as a `tbl`.
#' The result is accessed via
#' `rlresults(object, "noiseAnalysis")`.
#' @details
#'
#' ## Method
#'
#' The method used for noise analysis is a minor modification of the method
#' developed by [Diaz et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22499706/)
#' and also implemented by the
#' [deepTools](https://deeptools.readthedocs.io/en/develop/) function,
#' [plotFingerprint](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html).
#'
#' Briefly, if user-supplied [RLRanges] contain a bigWig coverage file,
#' then the coverage is quantified within random genomic regions
#' (`randomWindows`). The regions are then ranked. A good signal-to-noise
#' ratio will yield a distribution where most bins have little coverage
#' but a few have very high coverage. Use downstream tools like
#' `plotNoise` and `plotCompareNoise` to visualize these results.
#'
#'
#' @examples
#'
#' # Example RLRanges object
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # noiseAnalyze does not work on Windows OS
#' if (.Platform$OS.type != "windows") {
#'     # run noiseAnalyze
#'     rlr <- noiseAnalyze(rlr)
#' }
#'
#' @export
noiseAnalyze <- function(object,
    windows = NULL,
    force = FALSE) {

    # Check os
    if (.Platform$OS.type == "windows" & !force) {
        stop(
            "noiseAnalyze does not work on Windows.",
            " Override with force=TRUE."
        )
    }

    # Get genome
    genome <- GenomeInfoDb::genome(object)[1]
    stopifnot(genome %in% c("hg38", "hg19", "mm10"))

    # Get coverage
    coverage <- object@metadata$coverage
    if (coverage == "" || (!file.exists(coverage) && !urlExists(coverage))) {
        stop(
            "No coverage found. Content of 'coverage' slot: ",
            coverage, ". Set coverage with object@metadata$coverage <-"
        )
    }

    # Get random windows
    if (is.null(windows)) {
        windows <- randomWindows[[genome]]
    } else {
        warning(
            "Supplying custom genomic windows impairs",
            " ability to compare to RLBase samples. It is",
            " typically recommended to not set this option."
        )
    }

    # Get the coverage sum within each
    bw <- rtracklayer::import(
        con = rtracklayer::BigWigFile(object@metadata$coverage),
        selection = GenomicRanges::makeGRangesFromDataFrame(windows)
    )

    # Sum within ranges
    bdg <- valr::gr_to_bed(bw)
    bdg$score <- bw$score
    mappd <- valr::bed_map(windows, bdg, value = sum(.data$score))
    noise <- mappd %>%
        dplyr::arrange(dplyr::desc(.data$value)) %>%
        dplyr::mutate(
            value = .data$value / max(.data$value),
            rank = dplyr::row_number(.data$value)
        )

    # Add back to object
    methods::slot(object@metadata$results, "noiseAnalysis") <- noise

    return(object)
}
