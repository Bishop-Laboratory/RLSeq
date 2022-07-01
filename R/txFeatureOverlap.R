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
#' are in `tbl` format. For a full description of all columns in the output
#' table see [RLHub::???????????????????].
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
        annotations <- switch(paste0(genome),
            "hg38" = RLHub::annots_primary_hg38(quiet = TRUE),
            "mm10" = RLHub::annots_primary_mm10(quiet = TRUE)
        )
    }

    # Get TX annotations
    txfeat <- annotations[vapply(X = names(annotations), FUN = grepl, pattern = "Transcript_Features", FUN.VALUE = logical(1))]
    txfeat <- txfeat %>%
        bind_rows(.id = "feature") %>%
        mutate(
            feature = gsub(feature, pattern = ".+__", replacement = "")
        )

    # Convert to valr bed format
    rlbed <- object %>%
        as_tibble() %>%
        dplyr::mutate(seqnames = as.character(seqnames)) %>%
        dplyr::rename(
            chrom = seqnames,
            name = 6
        ) %>%
        dplyr::select(
            chrom, start, end, name
        )

    # Overlap with features
    olfeats <- valr::bed_intersect(rlbed, feats)

    # Get summarization
    olsum <- olfeats %>%
        dplyr::group_by(name.x) %>%
        dplyr::summarise(
            feats = paste0(feature.y, collapse = ", ")
        ) %>%
        dplyr::mutate(
            feature = case_when(
                grepl(feats, pattern = "TSS") ~ "TSS",
                grepl(feats, pattern = "fiveUTR") ~ "fiveUTR",
                grepl(feats, pattern = "Exon") ~ "Exon",
                grepl(feats, pattern = "Intron") ~ "Intron",
                grepl(feats, pattern = "threeUTR") ~ "threeUTR",
                grepl(feats, pattern = "TTS") ~ "TTS"
            )
        ) %>%
        dplyr::select(name = name.x, feature)

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
