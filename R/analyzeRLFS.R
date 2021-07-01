#' Analyze RLFS
#'
#' Analyzes the enrichment of peaks within R-loop forming sequences.
#'
#' @param peaks A GRanges object containing the R-loop ranges to check.
#' @param genome UCSC genome which peaks were generated from. Ignored if "chrom_sizes" or "RLFS" is specified. 
#' @param chrom_sizes A data.frame containing two columns with chromosome names
#' and chromosome sizes. 
#' @param mask A GRanges object containing the regions of the genome to mask. Can be generated from 
#' regioneR::getMask().
#' @param RLFS A GRanges object containing the R-loop forming sequences to 
#' compare against. 
#' @param ... Arguments passed to `regioneR::permTest()`
#' @return A named list containing the results of RegioneR, 
#' peaks annotated with the RLFS they overlap with, 
#' and the Z score results within 3000 BP of an RLFS. 
#' @examples
#' 
#' result <- RSeqR::analyzeRLFS(RSeqR::SRX1025890_peaks, genome="hg38")
#' 
#' @export
analyzeRLFS <- function(peaks, 
                        genome="hg38", 
                        chrom_sizes=NULL,
                        mask=NULL,
                        RLFS=NULL, 
                        ...) {
  
  # Check RLFS, chrom_sizes, and mask if not provided
  if (is.null(RLFS)) {RLFS <- RSeqR:::getRLFS(genome)}
  if (is.null(chrom_sizes)) {chrom_sizes <- RSeqR:::getChromSizes(genome)}
  if (is.null(mask)) {mask <- regioneR::getMask("hg38")}
    
  # Prevent stranded assignment
  GenomicRanges::strand(RLFS) <- "*"
  RLFS <- GenomicRanges::reduce(RLFS)
  
  # Final check
  stopifnot("data.frame" %in% class(chrom_sizes) &
              "GRanges" %in% class(mask) &
              "GRanges" %in% class(RLFS))
  
  # Run RLFS perm test
  pt <- regioneR::permTest(A=peaks, B=rlfs, 
                           genome=as.data.frame(chrom_sizes), 
                           mask=mask, 
                           allow.overlaps = FALSE,
                           randomize.function=regioneR::circularRandomizeRegions, 
                           evaluate.function=regioneR::numOverlaps, 
                           alternative = "greater",
                           ...)
  
  return(pt)
}





