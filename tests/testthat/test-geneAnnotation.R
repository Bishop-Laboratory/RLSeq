test_that(desc = "Test that gene annotation works", {
  
  # Load 
  load("rlr.rda")
  
  # Plot the enrich
  expect_s4_class(
    RLSeq::geneAnnotation(
      object = rlr, 
      quiet = TRUE
    ),
    "RLRanges"
  )
  
})
