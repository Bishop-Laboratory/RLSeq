test_that(desc = "RL Region Test works", {

    # Load
    rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))

    # Plot the enrich
    expect_s4_class(
        RLSeq::rlRegionTest(
            object = rlr
        ),
        "RLRanges"
    )
})
