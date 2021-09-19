#' Tests Genomic Feature Enrichment
#'
#' @description Tests the enrichment of genomic features in supplied R-loop peaks.
#'
#' @param peaks A GRanges object containing R-loop peaks
#' @param genome UCSC genome identifier indicating genome of supplied peaks.
#' @param annotations Annotation list. See details.
#' @param downsample If a numeric, data will be down sampled to the requested number of peaks.
#' This improves the speed of genomic shuffling and helps prevent p-value inflation.
#' If FALSE, then downsampling will not be performed. Default: 10000.
#' @param quiet If TRUE, messages will be suppressed. Default: False
#' @param cores Cores for use in parallel operations. Default: 1
#' the same format as RLSeq::annotationLst.
#' @details
#'   \strong{annotations}:
#'     A list which can be generated from the RLHub package or via local sources.
#'
#'     \emph{RLHub}: \code{RLHub['annots_hg38']} will return the full suite of annotations
#'       for the hg38 genome.
#'
#'     \emph{custom}: The only requirement is that \code{annotations} is a named
#'       list of \code{tbl}s in which each \code{tbl} follows the structure:
#'
#'       \tabular{lllll}{
#'         chrom \tab start \tab end \tab strand \tab name\cr
#'         chr1 \tab 10015 \tab 10498 \tab + \tab skewr__C_SKEW__1\cr
#'         chr1 \tab 10614 \tab 11380 \tab + \tab skewr__G_SKEW__1\cr
#'         ...
#'       }
#'
#'       Such a list can be generated directly from RLBase S3 bucket files:
#'
#'       \preformatted{
#'         small_anno <- list(
#'           "Centromeres" = readr::read_csv("https://rlbase-data.s3.amazonaws.com/annotations/hg38/Centromeres.csv.gz"),
#'           "SkewR" = readr::read_csv("https://rlbase-data.s3.amazonaws.com/annotations/hg38/SkewR.csv.gz")
#'           )
#'       }
#'
#' @return A tibble containing the results of the enrichment test.
#' @examples
#' small_anno <- list(
#'   "Centromeres" = readr::read_csv("https://rlbase-data.s3.amazonaws.com/annotations/hg38/Centromeres.csv.gz"),
#'   "SkewR" = readr::read_csv("https://rlbase-data.s3.amazonaws.com/annotations/hg38/SkewR.csv.gz")
#' )
#' RLSeq::featureEnrich(RLSeq::SRX1025890_peaks, annotations = small_anno, genome = "hg38")
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
featureEnrich <- function(peaks,
                          genome = c("hg38", "mm10"),
                          annotations,
                          downsample = 10000,
                          quiet = FALSE,
                          cores = 1) {

  # Cutoff for stats tests
  MIN_ROWS <- 200

  # Choose apply function
  applyFun <- pbapply::pblapply
  if (quiet) {
    applyFun <- lapply
  }
  if (cores > 1) {
    applyFun <- parallel::mclapply
    options(mc.cores = cores)
  }

  # Get annotations
  if (is.null(annotations)) {
    stop("Annotations must be supplied!")
  }

  # Get the genome
  chromSizes <- getChromSizes(genome) %>%
    dplyr::rename(chrom = .data$X1, size = .data$X2)

  # Wrangle the peaks
  toTest <- peaks %>%
    tibble::as_tibble() %>%
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
  seeds <- 1:1000
  while (TRUE) {
    seed <- seeds[1]
    toTestShuff <- try(valr::bed_shuffle(
      x = toTest,
      genome = chromSizes,
      seed = seed
    ), silent = TRUE)
    if (class(toTestShuff)[1] != "try-error") break
    seeds <- seeds[-1]
  }

  # Only keep annotations above cutoff
  annots <- annotations[sapply(annotations, nrow) > MIN_ROWS]

  # Test on all annotations
  if (!quiet) {
    message(" - Calculating enrichment...")
  }
  annoRes <- applyFun(
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
        tibble::tibble(
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

  return(annoRes)
}


#' Build peak statistics tibble
#'
#' A helper function for building the peak statistics tibble
#'
#' @param x The R-loop peaks to test.
#' @param xshuff x, but shuffled around the genome to build a control peakset.
#' @param y The annotations against which to test x.
#' @param chromSizeTbl A tibble containing the sizes of each chromosome in x and y.
#'
#' Columns should be "chrom" (chromosome names) and "size" (number of base pairs).
peak_stats <- function(x, xshuff, y, chromSizeTbl, quiet = FALSE) {

  # Cutoff for stats tests
  MIN_ROWS <- 200

  # Obtain distance test results (rel and abs). Based upon:
  # https://rnabioco.github.io/valr/articles/interval-stats.html
  reldist_rl <- valr::bed_reldist(x, y, detail = TRUE)
  if (!length(reldist_rl$chrom) | nrow(x) < MIN_ROWS | nrow(y) < MIN_ROWS) {
    if (!quiet) {
      warning("Not enough observations for interval tests...")
    }

    # Return results
    tibble::tibble(
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
        ks.test(reldist_rl$.reldist,
          reldist_shuf$.reldist,
          exact = FALSE
        )
      )
    } else {
      ks <- ks.test(reldist_rl$.reldist,
        reldist_shuf$.reldist,
        exact = FALSE
      )
    }

    pval_reldist <- ks %>%
      broom::tidy() %>%
      dplyr::pull("p.value")

    # Obtain fisher test results
    fshres_rl <- valr::bed_fisher(x, y, chromSizeTbl)
    fshres_shuf <- valr::bed_fisher(xshuff, y, chromSizeTbl)

    # Return results
    tibble::tibble(
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
