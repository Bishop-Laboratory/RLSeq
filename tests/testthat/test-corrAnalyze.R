test_that(desc = "Corr Analyze returns a numeric matrix", {
  expect_type( {
    BW_FILE <- paste0(RLSeq:::RLBASE_BW_URL, "SRX1025890_hg38.bw")
    RLSeq::corrAnalyze(coverage = BW_FILE, genome="hg38")
  },
    "double"
  )
})
