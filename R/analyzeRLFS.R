#' Analyze RLFS
#'
#' Analyzes the enrichment of peaks within R-loop forming sequences.
#'
#' @param peaks A GRanges object containing the R-loop ranges to check.
#' @param genome UCSC genome which peaks were generated from. Ignored if "chrom_sizes" or "RLFS" is specified. 
#' @param chrom_sizes A data.frame containing two columns with chromosome names
#' and chromosome sizes. 
#' @param RLFS A GRanges object containing the R-loop forming sequences to 
#' compare against. 
#' @param ... Arguments passed to `RegioneR::permTest()`
#' @return A named list containing the results of RegioneR, 
#' peaks annotated with the RLFS they overlap with, 
#' and the Z score results within 3000 BP of an RLFS. 
#' @examples
#' 
#' result <- RSeqR::analyzeRLFS(RSeqR::SRX1025890_peaks, genome="hg38")
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
analyzeRLFS <- function(peaks, 
                        genome="hg38", 
                        chrom_sizes=NULL,
                        RLFS=NULL, 
                        ...) {
  
  # Check whether genome supplied or chrom_sizes and RLFS supplied
  if (! is.character(genome) & (is.null(chrom_sizes) | is.null(RLFS))) {
    stop("User must supply genome or chrom_sizes and RLFS!")
  } else if (! is.null(genome)) {
    checkRes <- checkGenome(genome)
    RLFS <- getRLFS(genome)
    chrom_sizes <- getChromSizes(genome)
  } 
  
  # Check 
  
  
}





