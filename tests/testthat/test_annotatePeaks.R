
test_that(desc = "Test annotate_peaks", {
  expect_s4_class(
    RSeqR::annotatePeaks(ChIPpeakAnno::toGRanges("ERX2277510_E-MTAB-6318DRIP_mOHT_hg38.unstranded.bed",
                                                 format = "BED", header = FALSE)),
    "GRanges"
  )
})


