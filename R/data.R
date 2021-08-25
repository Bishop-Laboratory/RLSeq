#' List. Results of running RSeq on SRX1025890
"SRX1025890"

#' GRanges. Peaks of running RSeq on SRX1025890
"SRX1025890_peaks"

#' data.frame. Contains information about the genes available in UCSC database.
#' Result of running RSeqR:::buildAvailableGenomes()
"available_genomes"

#' List of GRanges. Contains genomic masks for several commons species.
#' See data-raw/masked_genomes.R for details.
"genomeMasks"

#' List of list of tibbles. Containing genomic annotations for use with RSeqR.
#' See data-raw/anotationLst.R for details.
"annotationLst"

#' A Caret model used for preprocessing datasets in predictCondition().
"prepFeatures"

#' A Caret model used for predicting whether a peakset resembles 
#' "case" (normal) or "control" (RNaseH1-like) conditions. 
#' Used in predictionCondition().
"fftModel"
