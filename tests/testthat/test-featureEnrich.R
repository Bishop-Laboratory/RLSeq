test_that(desc = "Feature Enrich returns a tibble", {
  expect_s3_class(
    RSeqR::featureEnrich(RSeqR::SRX1025890_peaks, genome="hg38"),
    "tbl_df"
  )
})


