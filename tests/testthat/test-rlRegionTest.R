test_that(desc = "RL Region Test returns as list", {
  expect_type (
    RSeqR::rlRegionTest(RSeqR::SRX1025890_peaks, genome="hg38"),
    "list"
  )
})


