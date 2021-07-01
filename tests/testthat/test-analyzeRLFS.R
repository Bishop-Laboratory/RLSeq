test_that(desc = "analyzeRLFS Arg Fail", {
  expect_error(
    RSeqR::analyzeRLFS(RSeqR::SRX1025890_peaks,
                       genome=NULL, 
                       chrom_sizes=NULL,
                       RLFS=NULL)
  )
})

test_that(desc = "analyzeRLFS Pass", {
  expect_s3_class(
    RSeqR::analyzeRLFS(RSeqR::SRX1025890_peaks,
                       genome=NULL, 
                       chrom_sizes=NULL,
                       RLFS=NULL),
    class = "permTestResults"
  )
})
