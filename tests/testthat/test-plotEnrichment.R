test_that(desc = "Test plotEnrichment", {

    # Load 
    load("rlr.rda")
    
    # Plot the enrich
    expect_type(
        RLSeq::plotEnrichment(
            object = rlr
        ),
        "list"
    )
    expect_type(
        RLSeq::plotEnrichment(
            object = rlr,
            splitby = "verdict"
        ),
        "list"
    )
    expect_type(
        RLSeq::plotEnrichment(
            object = rlr,
            splitby = "condType"
        ),
        "list"
    )
    expect_error(
        RLSeq::plotEnrichment(
            object = rlr,
            splitby = "ASDASD"
        )
    )
})
