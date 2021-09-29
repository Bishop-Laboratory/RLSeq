#' Tests Genomic Feature Enrichment
#'
#' @description Tests the enrichment of genomic features in supplied peaks.
#'
#' @param object An RLRanges object.
#' @param annotype The type of annotations to use. Options include "primary"
#'  (missing ChIP and eCLiP data) and "full" (all annotations). Default: "primary".
#' @param annotations Annotation list. See details.
#' @param downsample If a numeric, data will be down sampled to the requested
#'  number of peaks.
#' This improves the speed of genomic shuffling and 
#' helps prevent p-value inflation.
#' If FALSE, then downsampling will not be performed. Default: 10000.
#' @param quiet If TRUE, messages will be suppressed. Default: FALSE
#' the same format as RLSeq::annotationLst.
#' @details
#'
#' ## annotations 
#' 
#' \strong{RLHub}: \code{RLHub::annots_full_hg38()}, 
#' for example, will return the full suite of annotations for the hg38 genome.
#' 
#' \strong{Custom}: For custom annotations, the only requirement is that \code{annotations} is a
#' named list of \code{tbl}s in which each \code{tbl} follows the structure:
#'
#' \tabular{lllll}{
#'   chrom \tab start \tab end \tab strand \tab name\cr
#'   chr1 \tab 10015 \tab 10498 \tab + \tab skewr__C_SKEW__1\cr
#'   chr1 \tab 10614 \tab 11380 \tab + \tab skewr__G_SKEW__1\cr
#'   ...
#' }
#'
#' Such a list can be generated directly from RLBase S3 bucket files:
#'
#' \preformatted{
#'   small_anno <- list(
#'     "Centromeres" = readr::read_csv(
#'       paste0(
#'           "https://rlbase-data.s3.amazonaws.com",
#'           "/annotations/hg38/Centromeres.csv.gz"
#'       )
#'     ),
#'     "SkewR" = readr::read_csv(
#'       paste0(
#'           "https://rlbase-data.s3.amazonaws.com/",
#'           "annotations/hg38/SkewR.csv.gz")
#'     )
#'   )
#' }
#'
#' @return An RLRanges object containing the results of the enrichment test.
#' @examples
#'
#' # Example RLRanges dataset
#' rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))
#'
#' # RL Region Test
#' featureEnrich(rlr)
#'
#' # With custom annotations
#' small_anno <- list(
#'     "Centromeres" = readr::read_csv(
#'         paste0(
#'             "https://rlbase-data.s3.amazonaws.com/",
#'             "annotations/hg38/Centromeres.csv.gz"
#'         ),
#'         show_col_types = FALSE
#'     ),
#'     "SkewR" = readr::read_csv(
#'         paste0(
#'             "https://rlbase-data.s3.amazonaws.com/",
#'             "annotations/hg38/SkewR.csv.gz"
#'         ),
#'         show_col_types = FALSE
#'     )
#' )
#' featureEnrich(rlr, annotations = small_anno)
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
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
        if (! quiet) {
            annotations <- utils::getFromNamespace(
                paste0("annots_", annotype[1], "_", genome), "RLHub"
            )() 
        } else {
            annotations <- suppressMessages(utils::getFromNamespace(
                paste0("annots_", annotype[1], "_", genome), "RLHub"
            )())
        }
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
    # there's an issue in the shuffle for valr that necessitates this pattern...
    # Returns "maximum iterations exceeded in bed_shuffle" error
    # However, increasing number of attempts does nothing to fix this.
    # Only changing the seed seems to provide any workaround.
    seeds <- seq(1000)
    while (TRUE) {
        seed <- seeds[1]
        toTestShuff <- try(
            valr::bed_shuffle(
                x = toTest,
                genome = chromSizes,
                seed = seed
            ),
            silent = TRUE
        )
        if (!methods::is(toTestShuff, "try-error")) break
        seeds <- seeds[-1]
    }

    # Only keep annotations above cutoff
    annots <- annotations[vapply(annotations, nrow, numeric(1)) > MIN_ROWS]

    # Test on all annotations
    if (!quiet) {
        message(" - Calculating enrichment...")
    }
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

    if (!quiet) {
        message(" - Done")
    }

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
    if (! length(reldist_rl$chrom) | nrow(x) < MIN_ROWS | nrow(y) < MIN_ROWS) {
        if (!quiet) {
            warning("Not enough observations for interval tests...")
        }

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
        if (quiet) {
            ks <- suppressWarnings(
                stats::ks.test(
                    x = reldist_rl$.reldist,
                    y = reldist_shuf$.reldist,
                    exact = FALSE
                )
            )
        } else {
            ks <- stats::ks.test(
                x = reldist_rl$.reldist,
                y = reldist_shuf$.reldist,
                exact = FALSE
            )
        }


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
