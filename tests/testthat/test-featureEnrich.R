test_that(desc = "Feature Enrich returns a tibble", {
  small_anno <- list(
    "Centromeres" = readr::read_csv("https://rlbase-data.s3.amazonaws.com/annotations/hg38/Centromeres.csv.gz"),
    "SkewR" = readr::read_csv("https://rlbase-data.s3.amazonaws.com/annotations/hg38/SkewR.csv.gz")
  )
  expect_s3_class(
    RLSeq::featureEnrich(RLSeq::SRX1025890_peaks, annotations = small_anno, genome = "hg38"),
    "tbl_df"
  )
})
