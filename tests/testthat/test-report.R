test_that(desc = "RLSeq report", {
    # Load
    rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))

    # Plot the enrich
    expect_true(
        RLSeq::report(
            object = rlr,
            quiet = TRUE
        )
    )
})
