#' Annotate R-Loops with Genes
#'
#' Annotates R-loop peaks as a GRanges object with gene-level information
#'
#' @param peaks A GRanges object containing R-loop peaks
#' @param genome UCSC genome which peaks were generated from. 
#' Only "hg38" and "mm10" currently available. 
#' Use RSeqR::liftUtil() to convert to the correct format, if needed. Additionally,
#' this can be an EnsDb or TxDb object.
#' @return A tibble object containing annotated R-loop peaks
#' @examples
#' 
#' geneAnnotation(RSeqR::SRX1025890_peaks)
#' 
#' @export
geneAnnotation <- function(peaks, genome) {
  
  # Available genomes
  if (genome == "hg38") {
    if( ! requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
      stop("EnsDb.Hsapiens.v86 is required.",
           " Please install it with BiocManager::install('EnsDb.Hsapiens.v86')")
    }
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  } else if (genome == "mm10") {
    if( ! requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE)) {
      stop("EnsDb.Mmusculus.v79 is required.",
           " Please install it with BiocManager::install('EnsDb.Mmusculus.v79')")
    }
    edb <- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
  } else if (class(genome) %in% c("EnsDb", "TxDb")) {
    edb <- genome
  } else {
    stop("genome must be 'hg38' or 'mm10' or an object of class 'EnsDb' or 'TxDb'")
  }
  
  # Get the ensembl genes and conver to UCSC style
  edb <-  GenomicFeatures::genes(edb)
  GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
  
  # Wrangle EnsDb to tibble
  annoData <- edb %>%
    as.data.frame() %>%
    dplyr::select(chrom=seqnames, start, end, strand, gene_name) %>%
    tibble::as_tibble() %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    dplyr::mutate(chrom=as.character(chrom))
  
  # Wrangle peaks to tibble
  peaksIntersect <- peaks %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::select(chrom=seqnames, start, end, width) %>%
    dplyr::mutate(chrom=as.character(chrom))
  
  # Intersect
  anno <- valr::bed_intersect(peaksIntersect, annoData, suffix = c("__userPeaks", "__EnsGenes"))
  
  return(anno)
}
