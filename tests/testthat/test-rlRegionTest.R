test_that(desc = "Feature Enrich returns a tibble", {
  expect_type (
    RSeqR::rlRegionTest(RSeqR::SRX1025890_peaks, genome="hg38"),
    "list"
  )
})


