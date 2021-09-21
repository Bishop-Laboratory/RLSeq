test_that(desc = "Test that RLSeq works", {
  
  # Load 
  load("rlr.rda")
  
  # Plot the enrich
  expect_s4_class(
    RLSeq::RLSeq(
      object = rlr, 
      quiet = TRUE
    ),
    "RLRanges"
  )
  
})
