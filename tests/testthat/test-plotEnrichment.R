test_that(desc = "Test plotEnrichment", {



  # Load the RLBase enrichment results to compare with
  rlbase_enrich <- file.path(rlbase, "RLHub", "feature_enrichment_per_sample.rda")
  tmp <- tempfile()
  download.file(rlbase_enrich, destfile = tmp, quiet = TRUE)
  load(tmp)

  # Annotations can be found in RLHub
  annot <- file.path(rlbase, "RLHub", "annotations_primary_hg38.rda")
  tmp <- tempfile()
  download.file(annot, destfile = tmp, quiet = TRUE)
  load(tmp)

  # Perform test
  sampleRes <- featureEnrich(
    peaks = RLSeq::SRX1025890_peaks,
    genome = "hg38",
    annotations = annotations_primary,
    quiet = TRUE
  )

  # Plot the enrich
  plotEnrichment(
    sampleRes = annoResS96,
    rlbaseRes = feature_enrichment_per_sample
  )

  expect_type(RLSeq::plotEnrichment(sampleRes,
    rlbaseRes = rlbaseRes,
    splitby = "verdict", facetby = "none"
  ), "list")
  expect_type(RLSeq::plotEnrichment(sampleRes,
    rlbaseRes = rlbaseRes,
    splitby = "condType", facetby = "none"
  ), "list")
  expect_type(RLSeq::plotEnrichment(sampleRes,
    rlbaseRes = rlbaseRes
  ), "list")
  expect_error(RLSeq::plotEnrichment(sampleRes,
    rlbaseRes = rlbaseRes,
    yval = "ASDASDSA"
  ))
  file.remove("report.html")
})
