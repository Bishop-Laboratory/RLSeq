#' Annotate R-Loops with Genes
#'
#' Annotates R-loop peaks as a GRanges object with gene-level information
#'
#' @param peaks A GRanges object containing R-loop peaks
#' @return A GRanges object containing annotated R-loop peaks
#' @examples
#' 
#' geneAnnotation(RSeqR::SRX1025890_peaks)
#' 
#' @export
geneAnnotation <- function(peaks) {
  annoData <- ChIPpeakAnno::toGRanges(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
  anno <- ChIPpeakAnno::annotatePeakInBatch(peaks, AnnotationData = annoData,
                                            output = 'overlapping',
                                            select = 'all')
  mapping <- AnnotationDbi::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                   keys = AnnotationDbi::keys(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86),
                                   columns = "SYMBOL")
  annodf <- dplyr::left_join(as.data.frame(anno), mapping, by = c('feature' = "GENEID"))
  return(ChIPpeakAnno::toGRanges(annodf))
}
