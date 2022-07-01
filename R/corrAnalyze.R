#' Analyze Correlations
#'
#' Finds the pairwise correlation in signal around gold-standard R-Loop sites
#' between the query sample and the coverage tracks in the
#' [RLBase](https://gccri.bishop-lab.uthscsa.edu/rlbase/) database.
#' See *details*.
#'
#' Currently, this does not work on windows.
#'
#' @param object An RLRanges object.
#' @param force Force `corrAnalyze` to run, even if on Windows. Default: FALSE.
#' @return An RLRanges object with correlation results included as a `matrix`.
#' The correlation matrix is accessed via
#' `rlresults(object, "correlationMat")`.
#' @details
#'
#' ## Method
#'
#' The `corrAnalyze` function performs a correlation test that can
#' be used to assess sample-sample similarity by calculating coverage signal
#' (from genomic alignments) around “gold standard”
#' R-loop sites (PMID: [33411340](https://pubmed.ncbi.nlm.nih.gov/33411340/)).
#' The resulting correlation matrix is useful for determining how well a
#' supplied sample correlates with previously-published datasets.
#'
#' During the
#' [RLBase-data workflow](https://github.com/Bishop-Laboratory/RLBase-data),
#' the signal for each R-loop mapping sample
#' within “gold standard” R-loop sites was calculated see [RLHub::gs_signal].
#'
#' The `corrAnalyze` function loads [RLHub::gs_signal] and accepts an [RLRanges]
#' object with a valid coverage slot. It then does the following:
#'
#' 1. The coverage is quantified within the “gold standard”
#' sites and added as a column to the signal matrix from [RLHub::gs_signal].
#' 2. Then, the [stats::cor] function is used to calculate the Pearson
#' correlation pairwise between all samples, yielding a correlation matrix
#' 3. Finally, the correlation matrix is stashed in the in the
#' `correlationMat` slot of the [RLResults] and returned.
#'
#' @examples
#'
#' # Example RLRanges object
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # corrAnalyze does not work on Windows OS
#' if (.Platform$OS.type != "windows") {
#'
#'     # run corrAnalyze
#'     rlr <- corrAnalyze(rlr)
#' }
#' @export
corrAnalyze <- function(object, force = FALSE) {

    # Check os
    if (.Platform$OS.type == "windows" & !force) {
        stop(
            "corrAnalyze does not work on Windows.",
            " Override with force=TRUE."
        )
    }

    # Get genome
    genome <- GenomeInfoDb::genome(object)[1]
    stopifnot(genome == "hg38")

    # Get coverage
    coverage <- object@metadata$coverage
    if (coverage == "" || (!file.exists(coverage) && !urlExists(coverage))) {
        stop(
            "No coverage found. Content of 'coverage' slot: ",
            coverage, ". Set coverage with object@metadata$coverage <-"
        )
    }

    # Load GS signal
    gsSignalRLBase <- RLHub::gs_signal(quiet = TRUE)

    # Get the signal around GS R-loop sites
    bw <- getGSSignal(coverage, gssignal = gsSignalRLBase)

    # Wrangle BW into tibble
    bw <- bw %>%
        as.data.frame() %>%
        dplyr::as_tibble() %>%
        dplyr::rename(chrom = .data$seqnames) %>%
        dplyr::mutate(chrom = as.character(.data$chrom))

    # Get positions
    positions <- gsSignalRLBase %>%
        dplyr::select(.data$location) %>%
        dplyr::mutate(
            chrom = gsub(.data$location,
                pattern = "(.+)_(.+)_(.+)",
                replacement = "\\1"
            ),
            start = gsub(.data$location,
                pattern = "(.+)_(.+)_(.+)",
                replacement = "\\2"
            ),
            end = gsub(.data$location,
                pattern = "(.+)_(.+)_(.+)",
                replacement = "\\3"
            )
        ) %>%
        dplyr::select(-.data$location) %>%
        dplyr::mutate(dplyr::across(c("start", "end"), as.integer))

    # Summarize across gs intervals with valr
    bwMap <- valr::bed_map(x = positions, y = bw, value = sum(.data$score))

    # Remove any colnames that are shared with object sampleName
    gsSignalRLBase <- gsSignalRLBase[
        , which(colnames(gsSignalRLBase) != object@metadata$sampleName)
    ]

    # Combine with the original matrix
    combinedMat <- gsSignalRLBase %>%
        dplyr::inner_join(
            bwMap %>%
                dplyr::mutate(location = paste0(
                    .data$chrom, "_",
                    .data$start, "_",
                    .data$end
                )) %>%
                dplyr::select(.data$location, a_ = .data$value),
            by = "location"
        ) %>%
        dplyr::distinct(.data$location, .keep_all = TRUE)
    combinedMat <- as.data.frame(combinedMat)
    rownames(combinedMat) <- combinedMat$location
    combinedMat <- combinedMat[, -which(colnames(combinedMat) == "location")]
    combinedMat <- as.matrix(combinedMat)

    # Rename column
    colnames(combinedMat)[
        colnames(combinedMat) == "a_"
    ] <- object@metadata$sampleName

    # Find the correlation
    corMat <- stats::cor(combinedMat)

    # Add back to object
    methods::slot(object@metadata$results, "correlationMat") <- corMat

    return(object)
}
