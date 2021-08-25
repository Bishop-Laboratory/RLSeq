test_that(desc = "Predict Condition returns a list", {
  expect_is(
    RSeqR::predictCondition(RSeqR::analyzeRLFS(RSeqR::SRX1025890_peaks, genome="hg38")),
    "list"
  )
})

