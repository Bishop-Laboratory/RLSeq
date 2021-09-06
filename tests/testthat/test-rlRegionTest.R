test_that(desc = "RL Region Test returns as list", {
  expect_type (
    RLSeq::rlRegionTest(RLSeq::SRX1025890_peaks, genome="hg38"),
    "list"
  )
})


