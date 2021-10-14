#' This script creates the objects in `inst/extdata`

# Download hg38 centromere annotations
download.file(
    "https://rlbase-data.s3.amazonaws.com/annotations/hg38/Centromeres.csv.gz",
    destfile = "inst/extdata/Centromeres.csv.gz"
)

# Example dataset
rlbase <- "https://rlbase-data.s3.amazonaws.com"
pks <- file.path(rlbase, "peaks", "SRX7671349_hg38.broadPeak")
cvg <- file.path(rlbase, "coverage", "SRX7671349_hg38.bw")

# Download
download.file(pks, destfile = "inst/extdata/SRX7671349_hg38.broadPeak")

# Get RLRanges object
rlr <- RLSeq::RLRanges(pks,
    coverage = cvg, genome = "hg38", label="NEG",
    mode = "RDIP", sampleName = "RDIP-Seq +RNH1"
)

# Run RLSeq
rlr <- RLSeq::RLSeq(rlr)

# Save
dir.create("inst/extdata/", showWarnings = FALSE)
saveRDS(rlr, file = "inst/extdata/rlrsmall.rds", compress = "xz")
