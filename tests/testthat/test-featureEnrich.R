
test_that(desc = "Test that annotatePeaks returns a GRanges object", {
  expect_s3_class(
    RSeqR::featureEnrich(RSeqR::SRX1025890_peaks, genome="hg38"),
    "tbl_df"
  )
})


