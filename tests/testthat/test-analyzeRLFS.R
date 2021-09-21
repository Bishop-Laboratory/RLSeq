test_that(desc = "RLFS analysis", {
    # Load 
    load("rlr.rda")
    
    # Plot the enrich
    expect_s4_class(
        RLSeq::analyzeRLFS(
            object = rlr,
            quiet = TRUE
        ),
        "RLRanges"
    )
})
