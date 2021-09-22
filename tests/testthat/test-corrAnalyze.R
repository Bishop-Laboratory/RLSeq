test_that(desc = "Corr Analyze returns a numeric matrix", {
    # Load
    rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))

    # Plot the enrich
    expect_s4_class(
        RLSeq::corrAnalyze(
            object = rlr
        ),
        "RLRanges"
    )
})
