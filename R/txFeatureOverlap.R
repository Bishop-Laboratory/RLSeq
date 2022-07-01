#' Calculate overlap with transcript features
#'
#' Tests the overlap of transcript features with supplied peaks. See *details*.
#'
#' @param object An RLRanges object.
#' @param quiet If TRUE, messages will be suppressed. Default: FALSE
#' @details
#'
#' ## Method
#'
#' Transcript annotations were curated as part of the
#' [RLBase-data](https://github.com/Bishop-Laboratory/RLBase-data) workflow
#' and are provided via [RLHub::annotations].
#'
#' In `txFeatureOverlap`, each annotation "type" (e.g., "Exons", "Introns", etc)
#' is compared to the supplied RLRanges, yielding overlap statistics with
#' the following procedure:
#'
#' 1. For each annotation type, the peaks are overlapped with the annotations.
#' 2. Then the number of overlapping peaks is counted and summarised using a
#' priority order. This order determines which feature is assigned to a peak
#' when that peak overlaps multiple features. The order is "TSS", "TTS", "5'UTR",
#' "3'UTR", "Exon", "Intron", "Intergenic".
#'
#' @return An RLRanges object containing the results of the enrichment test
#' accessed via `rlresult(object, "txFeatureOverlap")`. The results
#' are in `tbl` format.
#'
#' @examples
#'
#' # Example RLRanges dataset
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # RL Region Test
#' txFeatureOverlap(rlr)
#'
#' @export
txFeatureOverlap <- function(object,
    quiet = FALSE) {

    # Check genome
    genome <- GenomeInfoDb::genome(object)[1]
    if (!genome %in% c("hg38", "mm10")) {
        stop(
            "Genome is not supported by built-in annotations. ",
            "Please supply custom ones of use one of hg38, mm10"
        )
    } else {
        annotations <- switch(genome,
            "hg38" = RLHub::annots_primary_hg38(quiet = TRUE),
            "mm10" = RLHub::annots_primary_mm10(quiet = TRUE)
        )
    }

    # Get TX annotations
    txfeat <- annotations[vapply(X = names(annotations), FUN = grepl, pattern = "Transcript_Features", FUN.VALUE = logical(1))]
    txfeat <- txfeat %>%
        dplyr::bind_rows(.id = "feature") %>%
        dplyr::mutate(
            feature = gsub(.data$feature, pattern = ".+__", replacement = "")
        )

    # Convert to valr bed format
    rlbed <- object %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(seqnames = as.character(.data$seqnames)) %>%
        dplyr::rename(
            chrom = .data$seqnames,
            name = 6
        ) %>%
        dplyr::select(
            .data$chrom, .data$start, .data$end, .data$name
        )

    # Overlap with features
    olfeats <- valr::bed_intersect(rlbed, txfeat)

    # Get summarization
    olsum <- olfeats %>%
        dplyr::group_by(.data$name.x) %>%
        dplyr::summarise(
            feats = paste0(.data$feature.y, collapse = ", ")
        ) %>%
        dplyr::mutate(
            feature = dplyr::case_when(
                grepl(.data$feats, pattern = "TSS") ~ "TSS",
                grepl(.data$feats, pattern = "TTS") ~ "TTS",
                grepl(.data$feats, pattern = "fiveUTR") ~ "fiveUTR",
                grepl(.data$feats, pattern = "threeUTR") ~ "threeUTR",
                grepl(.data$feats, pattern = "Exon") ~ "Exon",
                grepl(.data$feats, pattern = "Intron") ~ "Intron"
            )
        ) %>%
        dplyr::select(name = .data$name.x, .data$feature)

    # Add the intergenic
    olsum_final <- rlbed %>%
        dplyr::filter(!.data$name %in% olsum$name) %>%
        dplyr::mutate(
            feature = "Intergenic"
        ) %>%
        dplyr::select(
            .data$name, .data$feature
        ) %>%
        dplyr::bind_rows(olsum)

    if (!quiet) message(" - Done")

    # Add to object
    methods::slot(object@metadata$results,
        name = "txFeatureOverlap"
    ) <- olsum_final

    return(object)
}
