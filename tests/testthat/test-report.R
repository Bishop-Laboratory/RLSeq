test_that(desc = "RLSeq report", {
  # Load 
  load("rlr.rda")
  
  # Plot the enrich
  expect_true(
    RLSeq::report(
      object = rlr, 
      quiet = TRUE
    )
  )
})
