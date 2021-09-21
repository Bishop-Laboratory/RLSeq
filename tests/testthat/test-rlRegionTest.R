test_that(desc = "RL Region Test works", {
  
  # Load 
  load("rlr.rda")
  
  # Plot the enrich
  expect_s4_class(
    RLSeq::rlRegionTest(
      object = rlr
    ),
    "RLRanges"
  )
  
})
