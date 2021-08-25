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
    dplyr::select(seqnames, start, end)
  
  # Get the genome 
  chromSizes <- RSeqR::getChromSizes(genome) %>%
    dplyr::rename(chrom = X1, size = X2) 
  
  # Test on all annotations
  annoRes <- purrr::map(
    RSeqR::annotationLst[[genome]], 
    function(annoNow) {
      type <- annoNow$type[1]
      message("Testing - ", type)
      
      # Get shared seqnames -- important to avoid NAs
      shared_seqnames <- intersect(unique(annoNow$seqnames), 
                                   unique(toTest$seqnames))
      
      # Get the annotations
      x <- annoNow %>%
        dplyr::filter(seqnames %in% shared_seqnames) %>%
        dplyr::rename(chrom = seqnames) 
      
      # Get the peaks to test
      y <- toTest %>%
        dplyr::filter(seqnames %in% shared_seqnames) %>%
        dplyr::rename(chrom = seqnames) %>%
        dplyr::mutate(chrom = as.character(chrom))
      
      # Use the bed_projection test to get the overlap enrichment
      valr::bed_projection(x, y, chromSizes) %>%
        dplyr::mutate(type = !! type,
                      # So that we don't get -inf
                      obs_exp_ratio = log2(obs_exp_ratio + .001),
                      pmod = dplyr::case_when(
                        p.value == 0 ~ .Machine$double.xmin,
                        TRUE ~ p.value
                      ),
                      padj = -log10(
                        p.adjust(
                          pmod, method = "bonferroni"
                        )
                      ) * sign(obs_exp_ratio)) %>%
        dplyr::select(annotation=type, obs_exp_ratio, padj)
    }
  ) %>%
    dplyr::bind_rows() 
  
  return(annoRes)
  
}
