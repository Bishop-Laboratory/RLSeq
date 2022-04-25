test_that(desc = "Test that RLSeq works", {

    # Load
    rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))

    # Plot the enrich
    expect_s4_class(
        RLSeq::RLSeq(
            object = rlr,
            quiet = TRUE,
            useMask = FALSE,
            ntimes = 30
        ),
        "RLRanges"
    )
})
