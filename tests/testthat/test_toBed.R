
test_that(desc = "Test that toBed returns a list object and writes an output file", {
  # Create temporary file for output
  file <- tempfile()
  # Check that the file does not exist before running the funciton
  expect_false(file.exists(paste0(file, ".bed")))
  # Write .bed file and ascertain that toBed returns a list type object
  expect_type(
    RSeqR::toBed(ChIPpeakAnno::toGRanges("ERX2277510_E-MTAB-6318DRIP_mOHT_hg38.unstranded.bed",
                                                 format = "BED", header = FALSE), write = TRUE, filename = file),
    "list"
  )
  # Check that the output file now exists
  expect_true(file.exists(paste0(file, ".bed")))
  # Delete the output file
  file.remove(paste0(file, ".bed"))
})


