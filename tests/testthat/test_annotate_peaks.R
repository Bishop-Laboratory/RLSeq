
test_that(desc = "Test annotate_peaks", {
  expect_type(
    RSeqR::annotatePeaks("ERX2277510_E-MTAB-6318DRIP_mOHT_hg38.unstranded.bed"),
    "character"
  )
})


