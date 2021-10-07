test_that(desc = "Test that RLRanges() works", {
    rlbase <- "https://rlbase-data.s3.amazonaws.com"
    pks <- file.path(rlbase, "peaks", "SRX7671349_hg38.broadPeak")
    cvg <- file.path(rlbase, "coverage", "SRX7671349_hg38.bw")
    expect_s4_class(
        RLSeq::RLRanges(
            peaks = pks,
            coverage = cvg,
            mode = "RDIP",
            genome = "hg38",
            label = "NEG",
            sampleName = "RDIP +RNH"
        ),
        "RLRanges"
    )
})
