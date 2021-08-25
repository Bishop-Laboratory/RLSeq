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

test_that(desc = "Test that toBed returns a list object and writes an output file", {
  # Create temporary file for output
  file <- tempfile()
  # Check that the file does not exist before running the funciton
  expect_false(file.exists(paste0(file, ".bed")))
  # Write .bed file and ascertain that toBed returns a list type object
  expect_type(
    RSeqR::grangesToBed(RSeqR::SRX1025890_peaks, write = TRUE, filename = file),
    "list"
  )
  # Check that the output file now exists
  expect_true(file.exists(paste0(file, ".bed")))
  # Delete the output file
  file.remove(paste0(file, ".bed"))
})
