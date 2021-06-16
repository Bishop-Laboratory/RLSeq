
test_that(desc = "Test that toBed returns a list object", {
  expect_type(
    RSeqR::toBed(ChIPpeakAnno::toGRanges("ERX2277510_E-MTAB-6318DRIP_mOHT_hg38.unstranded.bed",
                                                 format = "BED", header = FALSE)),
    "list"
  )
})


