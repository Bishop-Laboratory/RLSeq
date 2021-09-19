test_that(desc = "analyzeRLFS Arg Fail", {
  expect_error(
    RLSeq::analyzeRLFS(RLSeq::SRX1025890_peaks,
      genome = NULL,
      chrom_sizes = NULL,
      RLFS = NULL
    )
  )
})

test_that(desc = "analyzeRLFS Pass", {
  expect_type(
    RLSeq::analyzeRLFS(RLSeq::SRX1025890_peaks,
      genome = "hg38",
      chrom_sizes = NULL,
      RLFS = NULL
    ),
    type = "list"
  )
})
