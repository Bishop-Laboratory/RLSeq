#' Tests Genomic Feature Enrichment
#'
#' Tests the enrichment of genomic features in supplied R-loop peaks.
#'
#' @param peaks A GRanges object containing R-loop peaks
#' @param genome UCSC genome identifier to use. Can be "hg38" or "mm10" currently. 
#' @return A tibble containing each annotation, the Log2 fold change of observed
#'  overlap vs expected, and the p.adjusted value.
#' @examples
#' 
#' RSeqR::featureEnrich(RSeqR::SRX1025890_peaks, genome="hg38")
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
featureEnrich <- function(peaks, genome=c("hg38", "mm10")) {
  
  # Wrangle the peaks
  toTest <- peaks %>%
    tibble::as_tibble() %>%
    dplyr::select(.data$seqnames, .data$start, .data$end)
  
  # Get the genome 
  chromSizes <- getChromSizes(genome) %>%
    dplyr::rename(chrom = .data$X1, size = .data$X2) 
  
  # Test on all annotations
  annoRes <- lapply(
    RSeqR::annotationLst[[genome]], 
    function(annoNow) {
      type <- annoNow$type[1]
      message("- - ", type)
      
      # Get shared seqnames -- important to avoid NAs
      shared_seqnames <- intersect(unique(annoNow$seqnames), 
                                   unique(toTest$seqnames))
      
      # Get the annotations
      x <- annoNow %>%
        dplyr::filter(.data$seqnames %in% shared_seqnames) %>%
        dplyr::rename(chrom = .data$seqnames) 
      
      # Get the peaks to test
      y <- toTest %>%
        dplyr::filter(.data$seqnames %in% shared_seqnames) %>%
        dplyr::mutate(chrom = as.character(.data$seqnames))
      
      # Use the bed_projection test to get the overlap enrichment
      valr::bed_projection(x, y, chromSizes) %>%
        dplyr::mutate(type = !! type,
                      # So that we don't get -inf
                      obs_exp_ratio = log2(.data$obs_exp_ratio + .001),
                      pmod = dplyr::case_when(
                        .data$p.value == 0 ~ .Machine$double.xmin,
                        TRUE ~ .data$p.value
                      ),
                      padj = -log10(
                        p.adjust(
                          .data$pmod, method = "bonferroni"
                        )
                      ) * sign(.data$obs_exp_ratio)) %>%
        dplyr::select(annotation=.data$type, .data$obs_exp_ratio, .data$padj)
    }
  ) %>%
    dplyr::bind_rows() 
  
  return(annoRes)
  
}
