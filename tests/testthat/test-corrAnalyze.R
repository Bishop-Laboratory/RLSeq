test_that(desc = "Corr Analyze returns a numeric matrix", {
  expect_type( {
    URL <- "https://rmapdb-data.s3.us-east-2.amazonaws.com/bigwigs/rseq-coverage-unstranded/SRX1025890_TC32_NT_DRIP.hg38.bw"
    BW_FILE <- "SRX1025890.bw"
    if (! file.exists(BW_FILE)) {
      download.file(URL, destfile=BW_FILE)
    }
    RSeqR::corrAnalyze(coverage = BW_FILE, genome="hg38")
  },
    "double"
  )
})
