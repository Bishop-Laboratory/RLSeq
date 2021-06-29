test_that(desc = "Test analyzeRLFS", {
  expect_error(
    RSeqR::analyzeRLFS(RSeqR::SRX1025890_peaks,
                       genome=NULL, 
                       chrom_sizes=NULL,
                       RLFS=NULL)
  )
  file.remove("report.html")
})

