#' GRanges. Peaks from running RSeq-CLI on SRX1025890 - in "hg38" format.
"SRX1025890_peaks"

#' GRanges. Peaks from running liftUtil(SRX1035890_peaks, "hg38", "hg19") in "hg19" format.
"SRX1025890_peaks_hg19"

#' data.frame. Contains information about the genes available in UCSC database.
#' Result of running RLSeq:::buildAvailableGenomes()
"available_genomes"

#' List of GRanges. Contains genomic masks for several commons species.
#' See data-raw/masked_genomes.R for details.
"genomeMasks"
