test_that(desc = "Predict Condition returns a list", {
  expect_type(
    RLSeq::predictCondition(RLSeq::analyzeRLFS(RLSeq::SRX1025890_peaks, genome = "hg38")),
    "list"
  )
})
