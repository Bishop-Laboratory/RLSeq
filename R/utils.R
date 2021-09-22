#' Check if URL exists
#' @param urlcon URL to check
#' @return logical. TRUE if status code 200, FALSE if not
urlExists <- function(urlcon) {
    res <- tryCatch(
        {
            con <- url(urlcon, open = "rb")
            close(con)
            return(TRUE)
        },
        error = function(e) {
            close(con)
            return(NULL)
        }
    )
    return(!is.null(res))
}


#' Get Chrom Sizes
#' Helper function which extracts chrom sizes from an RLRanges object.
#' @param object An RLRanges object.
#' @return A tibble containing chrom sizes
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
getChromSizes <- function(object) {
    cs <- GenomeInfoDb::seqinfo(object) %>%
        as.data.frame()
    cs$chrom <- rownames(cs)
    dplyr::as_tibble(cs) %>%
        dplyr::select(.data$chrom, size = .data$seqlengths)
}


#' Check RLFS
#' Helper function that checks whether a genome has RLFS available
#' @param genome the UCSC genome name to check
#' @return A logical, TRUE if available, FALSE if not
checkRLFSAnno <- function(genome) {
    return(
        urlExists(
            paste0(
                file.path(
                    RLBASE_URL,
                    "rlfs-beds/"
                ),
                genome,
                ".rlfs.bed"
            )
        )
    )
}


#' Get RLFS
#' Helper function that retrieves RLFS ranges
#' @param object An RLRanges object.
#' @return A GRanges object with RLFS for that species.
getRLFSAnno <- function(object) {

    # Get genome
    genome <- GenomeInfoDb::genome(object)[1]

    # Check if annotations available first
    if (!checkRLFSAnno(genome)) {
        stop("No RLFS annotations available for ", genome)
    }

    # Read in RLFS
    rlfs <- rtracklayer::import.bed(
        paste0(
            file.path(
                RLBASE_URL, "rlfs-beds/"
            ),
            genome,
            ".rlfs.bed"
        )
    )

    # Return as a GRanges object
    return(rlfs)
}


#' Get GS Signal
#'
#' Extract signal around GS R-loop sites
#' @param coverage The path to a .bigWig file (can be a URL)
#' @param gssignal The GS signal obtained from RLHub.
#' @return A named list containing the results of correlation analysis.
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
getGSSignal <- function(coverage, gssignal) {
    # Get the locations of the gs sites
    positions <- gssignal$location
    positions <- dplyr::tibble(location = positions) %>%
        dplyr::mutate(
            seqnames = gsub(.data$location,
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
        GenomicRanges::makeGRangesFromDataFrame()

    # Read in the bigWig file using these locations
    bw <- rtracklayer::import(
        con = rtracklayer::BigWigFile(coverage),
        selection = positions
    )
}


#' Table to Regions
#'
#' Helper function to Convert "table" format to "regions" format.
#' @param table A tibble in "Table" format from RLHub.
#' @return A tibble in "regions" format.
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
#' @export
tableToRegions <- function(table) {
    locpat <- "(.+):(.+)\\-(.+):(.+)"
    table %>%
        dplyr::mutate(
            chrom = as.character(
                gsub(.data$location, pattern = locpat, replacement = "\\1")
            ),
            start = as.numeric(
                gsub(.data$location, pattern = locpat, replacement = "\\2")
            ),
            end = as.numeric(
                gsub(.data$location, pattern = locpat, replacement = "\\3")
            ),
            strand = as.character(
                gsub(.data$location, pattern = locpat, replacement = "\\4")
            ),
            strand = dplyr::case_when(
                .data$strand == "." ~ "*",
                TRUE ~ .data$strand
            )
        ) %>%
        dplyr::select(
            .data$chrom,
            .data$start,
            .data$end,
            .data$strand,
            name = .data$rlregion
        ) %>%
        dplyr::distinct()
}
