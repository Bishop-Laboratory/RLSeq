test_that(desc = "Corr Analyze returns a numeric matrix", {
    # Load 
    load("rlr.rda")
    
    # Plot the enrich
    expect_s4_class(
        RLSeq::corrAnalyze(
            object = rlr
        ),
        "RLRanges"
    )
})
