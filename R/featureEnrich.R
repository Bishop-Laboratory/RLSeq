#' Test Genomic Feature Enrichment
#'
#' Tests the enrichment of genomic features in supplied peaks. See *details*.
#'
#' @param object An RLRanges object.
#' @param annotype The type of annotations to use. Can be one of "primary"
#' or "full". Default: "primary". See
#' [RLHub::annotations] for greater detail.
#' @param annotations A custom annotation list of the same structure
#' described in [RLHub::annotations].
#' @param downsample If a numeric, data will be down sampled to the requested
#' number of peaks. This improves the speed of genomic shuffling and
#' helps prevent p-value inflation.
#' If FALSE, then downsampling will not be performed. Default: 10000.
#' @param quiet If TRUE, messages will be suppressed. Default: FALSE
#' @details
#'
#' ## Method
#'
#' Annotations relevant to R-loops were curated as part of the
#' [RLBase-data](https://github.com/Bishop-Laboratory/RLBase-data) workflow
#' and are provided via [RLHub::annotations].
#'
#' In `featureEnrich`, each annotation "type" (e.g., "Exons", "Introns", etc)
#' is compared to the supplied RLRanges, yielding enrichment statistics with
#' the following procedure:
#'
#' 1. For each annotation type, the peaks are overlapped with the annotations.
#' 2. Then, [valr::bed_reldist] is used to find the relative distance
#' distribution between the peaks and the annotations for both the supplied
#' RLRanges and shuffled RLRanges (via [valr::bed_shuffle]).
#' Significance of the relative distance is calculated via [stats::ks.test].
#' 3. Then, Fisher’s exact test is implemented via [valr::bed_fisher]
#' to obtain the significance of the overlap and the odds ratio.
#'
#' @return An RLRanges object containing the results of the enrichment test
#' accessed via `rlresult(object, "featureEnrichment")`. The results
#' are in `tbl` format. For a full description of all columns in the output
#' table see [RLHub::feat_enrich_samples].
#'
#' @examples
#'
#' # Example RLRanges dataset
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # RL Region Test
#' featureEnrich(rlr)
#'
#' # With custom annotations
#' small_anno <- list(
#'     "Centromeres" = readr::read_csv(
#'         system.file("extdata", "Centromeres.csv.gz", package = "RLSeq"),
#'         show_col_types = FALSE
#'     )
#' )
#' featureEnrich(rlr, annotations = small_anno)
#' @export
featureEnrich <- function(object,
    annotype = c("primary", "full"),
    annotations = NULL,
    downsample = 10000,
    quiet = FALSE) {

    # Cutoff for stats tests
    MIN_ROWS <- 200

    # Check genome
    genome <- GenomeInfoDb::genome(object)[1]
    if (!genome %in% c("hg38", "mm10") & is.null(annotations)) {
        stop(
            "Genome is not supported by built-in annotations. ",
            "Please supply custom ones of use one of hg38, mm10"
        )
    } else if (is.null(annotations)) {
        annotations <- switch(paste0(annotype[1], "_", genome),
            "primary_hg38" = RLHub::annots_primary_hg38(quiet = TRUE),
            "primary_mm10" = RLHub::annots_primary_mm10(quiet = TRUE),
            "full_hg38" = RLHub::annots_full_hg38(quiet = TRUE),
            "full_mm10" = RLHub::annots_full_mm10(quiet = TRUE)
        )
    }

    # Get annotations
    if (is.null(annotations)) {
        stop("Annotations must be supplied!")
    }

    # Get the genome
    chromSizes <- getChromSizes(object)

    # Wrangle the peaks
    toTest <- object %>%
        dplyr::as_tibble() %>%
        dplyr::select(chrom = .data$seqnames, .data$start, .data$end)

    # Downsample
    if (is.numeric(downsample) && nrow(toTest) > downsample) {
        toTest <- dplyr::sample_n(toTest, size = downsample)
    }

    # Get shuffle
    toTestShuff <- regioneR::randomizeRegions(
        A = as.data.frame(toTest), genome = as.data.frame(chromSizes)
    ) %>%
        dplyr::as_tibble() %>%
        dplyr::select(chrom = .data$seqnames, .data$start, .data$end)

    # Only keep annotations above cutoff
    annots <- annotations[vapply(annotations, nrow, numeric(1)) > MIN_ROWS]

    # Test on all annotations
    if (!quiet) message(" - Calculating enrichment...")

    annoRes <- lapply(
        seq(annots),
        function(j) {
            annoSubNow <- annots[[j]]
            typeNow <- names(annots)[j]
            # Get shared seqnames -- important to avoid NAs
            shared_seqnames <- intersect(
                unique(annoSubNow$chrom),
                unique(toTest$chrom)
            )

            # Get the peaks to test
            x <- toTest %>%
                dplyr::filter(.data$chrom %in% shared_seqnames) %>%
                dplyr::mutate(chrom = as.character(.data$chrom))

            # Get shuffled peaks
            xshuff <- toTestShuff %>%
                dplyr::filter(.data$chrom %in% shared_seqnames) %>%
                dplyr::mutate(chrom = as.character(.data$chrom))

            # Get the annotations
            y <- annoSubNow %>%
                dplyr::filter(.data$chrom %in% shared_seqnames)

            # Get shared chromSizes
            chromSizesNow <- chromSizes %>%
                dplyr::filter(.data$chrom %in% shared_seqnames)

            # Use the peak_stats test to get the enrichment
            pat <- "(.+)__(.+)"
            y <- valr::bed_merge(y)
            pkstats <- peak_stats(x, xshuff, y, chromSizesNow, quiet = quiet)
            dplyr::bind_cols(
                dplyr::tibble(
                    db = gsub(typeNow, pattern = pat, replacement = "\\1"),
                    type = gsub(typeNow, pattern = pat, replacement = "\\2"),
                    num_tested_peaks = nrow(y),
                    num_total_peaks = nrow(toTest),
                    num_tested_anno_ranges = nrow(x),
                    num_total_anno_ranges = nrow(annoSubNow)
                ),
                pkstats
            )
        }
    ) %>%
        dplyr::bind_rows()

    if (!quiet) message(" - Done")

    # Add to object
    methods::slot(object@metadata$results,
        name = "featureEnrichment"
    ) <- annoRes

    return(object)
}


#' Build peak statistics tibble
#'
#' A helper function for building the peak statistics tibble
#'
#' @param x The R-loop peaks to test.
#' @param xshuff x, but shuffled around the genome to build a control peakset.
#' @param y The annotations against which to test x.
#' @param chromSizeTbl A tibble containing the sizes of each
#' chromosome in x and y.
#' @param quiet If TRUE, messages will be suppressed. Default: FALSE
#' @return A tibble containing the test results.
peak_stats <- function(x, xshuff, y, chromSizeTbl, quiet = FALSE) {

    # Cutoff for stats tests
    MIN_ROWS <- 200

    # Obtain distance test results (rel and abs). Based upon:
    # https://rnabioco.github.io/valr/articles/interval-stats.html
    reldist_rl <- valr::bed_reldist(x, y, detail = TRUE)
    if (!length(reldist_rl$chrom) | nrow(x) < MIN_ROWS | nrow(y) < MIN_ROWS) {
        if (!quiet) warning("Not enough observations for interval tests...")

        # Return results
        dplyr::tibble(
            avg_reldist_rl = NA,
            avg_reldist_shuf = NA,
            pval_reldist = NA,
            stat_fisher_rl = NA,
            pval_fisher_rl = NA,
            stat_fisher_shuf = NA,
            pval_fisher_shuf = NA
        )
    } else {
        reldist_shuf <- valr::bed_reldist(xshuff, y, detail = TRUE)
        ks <- stats::ks.test(
            x = reldist_rl$.reldist,
            y = reldist_shuf$.reldist,
            exact = FALSE
        )


        pval_reldist <- ks$p.value

        # Obtain fisher test results
        fshres_rl <- valr::bed_fisher(x, y, chromSizeTbl)
        fshres_shuf <- valr::bed_fisher(xshuff, y, chromSizeTbl)

        # Return results
        dplyr::tibble(
            avg_reldist_rl = mean(reldist_rl$.reldist),
            avg_reldist_shuf = mean(reldist_shuf$.reldist),
            pval_reldist = ifelse(pval_reldist == 0, 2.2E-16, pval_reldist),
            stat_fisher_rl = fshres_rl$estimate,
            stat_fisher_shuf = fshres_shuf$estimate,
            pval_fisher_rl = ifelse(fshres_rl$p.value == 0,
                .Machine$double.xmin, fshres_rl$p.value
            ),
            pval_fisher_shuf = ifelse(fshres_shuf$p.value == 0,
                .Machine$double.xmin, fshres_shuf$p.value
            )
        )
    }
}
