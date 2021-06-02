#' Annotate Peaks
#'
#' Annotates DRIP-Seq peaks with gene-level information
#'
#' @param peaks A .bed file containing DRIP-Seq peaks
#' @return A .bed file containing annotated DRIP-Seq peaks
#' @export

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("ChIPpeakAnno")

#library(tidyverse)
#library(ChIPpeakAnno)
#library(EnsDb.Hsapiens.v86)
#library(org.Hs.eg.db)

#bed <- "R/ERX2277510_E-MTAB-6318DRIP_mOHT_hg38.unstranded.bed"
#peaks <- ChIPpeakAnno::toGRanges(bed, format = "BED", header = FALSE)

annotatePeaks <- function(peaks) {
  annoData <- toGRanges(EnsDb.Hsapiens.v86)

  anno <- ChIPpeakAnno::annotatePeakInBatch(peaks, AnnotationData = annoData, output = 'overlapping', select = 'all')

  anno <- ChIPpeakAnno::addGeneIDs(anno, orgAnn = "org.Hs.eg.db", feature_id_type = "ensembl_gene_id", IDs2Add = c("symbol"))

  anno_df <- as.data.frame(anno)
  names(anno_df)[1] <- '#Seqnames'
  anno_df['peak'] <- make.unique(anno_df[['peak']], sep = '_')

  return(anno_df)
}

#res <- annotatePeaks(peaks)
