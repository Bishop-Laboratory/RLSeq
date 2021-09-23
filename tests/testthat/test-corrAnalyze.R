test_that(desc = "Corr Analyze returns a numeric matrix", {
    # Load
    rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))

    # Get correlation. Should error if os is Windows.
    if (.Platform$OS.type != "windows") {
        expect_s4_class(
            RLSeq::corrAnalyze(
                object = rlr
            ),
            "RLRanges"
        )
    } else {
        expect_error(
            RLSeq::corrAnalyze(
                object = rlr
            )
        )
    }
    
})
