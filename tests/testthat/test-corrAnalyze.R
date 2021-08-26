test_that(desc = "Corr Analyze returns a numeric matrix", {
  expect_type( {
    BW_FILE <- "https://rmapdb-data.s3.us-east-2.amazonaws.com/bigwigs/rseq-coverage-unstranded/SRX1025890_hg38.bw"
    RSeqR::corrAnalyze(coverage = BW_FILE, genome="hg38")
  },
    "double"
  )
})
