#' Annotate R-Loops
#'
#' Annotates DRIP-Seq peaks as a GRanges object with gene-level information
#'
#' @param peaks A GRanges object containing DRIP-Seq peaks
#' @return A GRanges object containing annotated DRIP-Seq peaks
#' @export

# TO DO: 1. Example in doc at top of file.

#bed <- "tests/testthat/ERX2277510_E-MTAB-6318DRIP_mOHT_hg38.unstranded.bed"
#peaks <- ChIPpeakAnno::toGRanges(bed, format = "BED", header = FALSE)

annotateRLoops <- function(peaks) {
  annoData <- ChIPpeakAnno::toGRanges(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
  anno <- ChIPpeakAnno::annotatePeakInBatch(peaks, AnnotationData = annoData,
                                            output = 'overlapping',
                                            select = 'all')
  mapping <- AnnotationDbi::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                   keys = AnnotationDbi::keys(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86),
                                   columns = "SYMBOL")
  annodf <- dplyr::left_join(as.data.frame(anno), mapping, by = c('feature' = "GENEID"))
  return(plyranges::as_granges(annodf))
}

#res <- annotateRLoops(peaks)
#res
