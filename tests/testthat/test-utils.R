test_that("getRLFSAnno works", {
  expect_s4_class(RLSeq:::getRLFSAnno("hg38"), "GRanges")
})

test_that("urlExists works", {
  expect_true(RLSeq:::urlExists("www.google.com"))
})

test_that("getChromSizes works", {
  expect_s3_class(RLSeq:::getChromSizes("hg38"), class = "tbl")
})

test_that("checkRLFSAnno works", {
  expect_true(RLSeq:::checkRLFSAnno("hg38"))
})

test_that("getChain works", {
  expect_s4_class(RLSeq:::getChain("hg19", "hg38"), "Chain")
})

test_that("liftUtil works", {
  expect_s4_class(
    RLSeq:::liftUtil(RLSeq::SRX1025890_peaks_hg19, "hg19", "hg38"),
    "GRanges"
  )
})

test_that(desc = "Test that toBed returns a list object and writes an output file", {
  # Create temporary file for output
  file <- tempfile()
  # Check that the file does not exist before running the funciton
  expect_false(file.exists(paste0(file, ".bed")))
  # Write .bed file and ascertain that toBed returns a list type object
  expect_type(
    RLSeq::grangesToBed(RLSeq::SRX1025890_peaks, write = TRUE, filename = file),
    "list"
  )
  # Check that the output file now exists
  expect_true(file.exists(paste0(file, ".bed")))
  # Delete the output file
  file.remove(paste0(file, ".bed"))
})

test_that("getGSSignal works", {
  coverage <- paste0(RLSeq:::RLBASE_BW_URL, "SRX1025890_hg38.bw")
  expect_s4_class(RLSeq:::getGSSignal(coverage = coverage, "hg38"), class = "GRanges")
})
