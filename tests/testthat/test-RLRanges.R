test_that(desc = "Test that RLRanges() works", {
  pks <- "https://rlbase-data.s3.amazonaws.com/peaks/SRX1025890_hg38.broadPeak"
  cvg <- "https://rlbase-data.s3.amazonaws.com/coverage/SRX1025890_hg38.bw"
  expect_s4_class(
    RLRanges(
      peaks = pks,
      coverage = cvg,
      mode = "DRIP",
      genome = "hg38",
      condType = "POS",
      sampleName = "TC32 DRIP-Seq"
    ),
    "tbl_df"
  )
})
