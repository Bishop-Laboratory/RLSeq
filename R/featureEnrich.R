#' Tests Genomic Feature Enrichment
#'
#' Tests the enrichment of genomic features in supplied R-loop peaks.
#'
#' @param peaks A GRanges object containing R-loop peaks
#' @param genome UCSC genome identifier to use. Can be "hg38" or "mm10" currently. 
#' @param annotations Custom annotation list to use instead of built-in. It must follow
#' @param quiet If TRUE, messages will be suppressed. Default: False
#' @param cores Cores for use in parallel operations. Default: 1
#' the same format as RLSeq::annotationLst.
#' @return A tibble containing each annotation, the Log2 fold change of observed
#'  overlap vs expected, and the p.adjusted value.
#' @examples
#' 
#' RLSeq::featureEnrich(RLSeq::SRX1025890_peaks, genome="hg38")
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
featureEnrich <- function(peaks,
                          genome=c("hg38", "mm10"), 
                          annotations=NULL,
                          quiet = FALSE,
                          cores = 1) {
  
  # Choose apply function
  applyFun <- ifelse(quiet, lapply, pbapply::pblapply)
  if (cores > 1) {
    applyFun <- parallel::mclapply
    options(mc.cores = cores)
  }
  
  # Get the genome 
  chromSizes <- getChromSizes(genome) %>%
    dplyr::rename(chrom = .data$X1, size = .data$X2)
  
  # Wrangle the peaks
  toTest <- peaks %>%
    tibble::as_tibble() %>%
    dplyr::select(chrom = .data$seqnames, .data$start, .data$end)
  
  # Get shuffle
  toTestShuff <- valr::bed_shuffle(toTest, chromSizes)
  
  # Get the annotation list
  if (! is.null(annotations)) {
    annotationGen <- annotations[[genome]]
  } else {
    annotationGen <- RLSeq::annotationLst[[genome]]
  }
  
  # Flatten annotations into tbl list
  subpat <- "(.+)__(.+)__(.+)"
  if (! quiet) {
    message(" - Preparing annotations...")
  }
  annots <- lapply(annotationGen, FUN = function(x) {
    if ("tbl" %in% class(x)) {
      x
    } else {
      dplyr::bind_rows(x)
    }
  }) %>% 
    dplyr::bind_rows() %>%
    dplyr::mutate(
      comb = gsub(.data$name, pattern = subpat, 
                  replacement = "\\1__\\2", perl=TRUE)
    ) %>% dplyr::group_by(.data$comb) %>%
    {setNames(dplyr::group_split(.), dplyr::group_keys(.)[[1]])}
  annots <- annots[sapply(annots, nrow) > 200]
  
  # Test on all annotations
  # For each element in the database, calculate the projection statistics
  if (! quiet) {
    message(" - Calculating enrichment...")
  }
  annoRes <- applyFun(
    seq(annots), 
    function(j) {
      annoSubNow <- annots[[j]]
      typeNow <- names(annots)[j]
      # Get shared seqnames -- important to avoid NAs
      shared_seqnames <- intersect(unique(annoSubNow$chrom), 
                                   unique(toTest$chrom))
      
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
      
      # Use the bed_projection test to get the overlap enrichment
      pat <- "(.+)__(.+)"
      pkstats <- peak_stats(x, xshuff, y, chromSizesNow)
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
    }) %>%
    dplyr::bind_rows()
  
  if (! quiet) {
    message(" - Done")
  }
  
  return(annoRes)
}


#' Build peak statistics tibble
#' 
#' A helper function for building the peak statistics tibble
#' 
#' @param x The R-loop peaks to test.
#' @param y The annotations against which to test x.
#' @param chromSizeTbl A tibble containing the sizes of each chromosome in x and y.
#' Columns should be "chrom" (chromosome names) and "size" (number of base pairs).
peak_stats <- function(x, xshuff, y, chromSizeTbl) {
  
  # Obtain distance test results (rel and abs). Based upon:
  # https://rnabioco.github.io/valr/articles/interval-stats.html
  reldist_rl <- valr::bed_reldist(x, y, detail = TRUE)
  if (! length(reldist_rl$chrom) | nrow(x) < 200 | nrow(y) < 200) {
    warning("Not enough observations for interval tests...")
    # Return results
    tibble::tibble(
      avg_reldist_rl = NA,
      avg_reldist_shuf = NA,
      pval_reldist = NA,
      stat_fisher_rl = NA,
      pval_fisher_rl = NA,
      stat_fisher_shuf = NA,
      pval_fisher_shuf = NA,
      stat_projection_rl = NA,
      pval_projection_rl = NA,
      stat_projection_shuf = NA,
      pval_projection_shuf = NA
    )
    
  } else {
    reldist_shuf <- valr::bed_reldist(xshuff, y, detail = TRUE)
    pval_reldist <- ks.test(reldist_rl$.reldist, 
                            reldist_shuf$.reldist) %>%
      broom::tidy() %>% 
      dplyr::pull("p.value")
    
    # Obtain fisher test results
    fshres_rl <- valr::bed_fisher(x, y, chromSizeTbl) 
    fshres_shuf <- valr::bed_fisher(xshuff, y, chromSizeTbl) 
    # Obtain projection test results
    projres_rl <- valr::bed_projection(x, y, chromSizeTbl)
    projres_shuf <- valr::bed_projection(xshuff, y, chromSizeTbl)
    
    # Return results
    tibble::tibble(
      avg_reldist_rl = mean(reldist_rl$.reldist),
      avg_reldist_shuf = mean(reldist_shuf$.reldist),
      pval_reldist = ifelse(pval_reldist == 0, 2.2E-16, pval_reldist),
      stat_fisher_rl = fshres_rl$estimate,
      stat_fisher_shuf = fshres_shuf$estimate,
      pval_fisher_rl = ifelse(fshres_rl$p.value == 0, 
                              .Machine$double.xmin, pval_reldist),
      pval_fisher_shuf = ifelse(fshres_shuf$p.value == 0, 
                                .Machine$double.xmin, pval_reldist)
    )
  }
} 


