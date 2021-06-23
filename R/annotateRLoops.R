#' Annotate R-Loops
#'
#' Annotates DRIP-Seq peaks as a GRanges object with gene-level information
#'
#' @param peaks A GRanges object containing DRIP-Seq peaks
#' @return A GRanges object containing annotated DRIP-Seq peaks
#' @export

#bed <- "tests/testthat/ERX2277510_E-MTAB-6318DRIP_mOHT_hg38.unstranded.bed"
#peaks <- ChIPpeakAnno::toGRanges(bed, format = "BED", header = FALSE)

annotateRLoops <- function(peaks) {
  annoData <- ChIPpeakAnno::toGRanges(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86) # Use this database for gene IDs as well using select from annotationDbi
  anno <- ChIPpeakAnno::annotatePeakInBatch(peaks, AnnotationData = annoData,
                                            output = 'overlapping',
                                            select = 'all')
  anno <- ChIPpeakAnno::addGeneIDs(anno, orgAnn = "org.Hs.eg.db",
                                   feature_id_type = "ensembl_gene_id",
                                   IDs2Add = c("symbol"))
  return(anno)
}

#res <- annotateRLoops(peaks)
