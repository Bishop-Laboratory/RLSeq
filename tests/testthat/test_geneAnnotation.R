
test_that(desc = "Test that annotatePeaks returns a GRanges object", {
  expect_s4_class(
    RSeqR::geneAnnotation(ChIPpeakAnno::toGRanges("SRX1025890_TC32_NT_DRIP_hg38.unstranded.broadPeak",
                                                  format = "BED", header = FALSE)),
    "GRanges"
  )
})


