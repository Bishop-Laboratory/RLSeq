#' This script creates the \code{inst/ext-data/rlrsmall.rds} data object

# Example dataset
rlbase <- "https://rlbase-data.s3.amazonaws.com"
pks <- file.path(rlbase, "peaks", "SRX7671349_hg38.broadPeak")
cvg <- file.path(rlbase, "coverage", "SRX7671349_hg38.bw")

# Get RLRanges object
rlr <- RLRanges(pks,
    coverage = cvg, genome = "hg38", label="NEG",
    mode = "RDIP", sampleName = "RDIP-Seq +RNH1"
)

# Run RLSeq
rlr <- RLSeq::RLSeq(rlr)

# Save
saveRDS(rlr, file = "inst/ext-data/rlrsmall.rds", compress = "xz")
