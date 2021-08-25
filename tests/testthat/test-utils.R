test_that("getRLFSAnno works", {
  expect_s4_class(RSeqR:::getRLFSAnno("hg38"), "GRanges")
})

test_that("urlExists works", {
  expect_true(RSeqR:::urlExists("www.google.com"))
})

test_that("getChromSizes works", {
  expect_s3_class(RSeqR:::getChromSizes("hg38"), class="tbl")
})

test_that("checkRLFSAnno works", {
  expect_true(RSeqR:::checkRLFSAnno("hg38"))
})

test_that("getChain works", {
  expect_s4_class(RSeqR:::getChain("hg19", "hg38"), "Chain")
})

test_that("liftUtil works", {
  expect_s4_class(RSeqR:::liftUtil(RSeqR::SRX1025890_peaks_hg19, "hg19", "hg38"),
                  "GRanges")
})
