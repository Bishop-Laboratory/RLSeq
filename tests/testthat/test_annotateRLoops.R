
test_that(desc = "Test that annotatePeaks returns a GRanges object", {
  expect_s4_class(
    RSeqR::annotateRLoops(ChIPpeakAnno::toGRanges("ERX2277510_E-MTAB-6318DRIP_mOHT_hg38.unstranded.bed",
                                                  format = "BED", header = FALSE)),
    "GRanges"
  )
})


