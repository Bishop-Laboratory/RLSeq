test_that(desc = "Predict Condition works", {
  # Load 
  load("rlr.rda")
  
  # Predict condition
  expect_s4_class(
    RLSeq::predictCondition(
      object = rlr
    ),
    "RLRanges"
  )
  
  # Predict condition attempt without RLFS (should fail)
  slot(rlr@metadata$results, "rlfsRes", check=FALSE) <- NULL
  expect_error(
    RLSeq::predictCondition(
      object = rlr
    )
  )
  
})
