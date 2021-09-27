test_that(desc = "Analyze RLFS works", {

    # Load
    rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))

    # With no Z
    expect_s4_class(
        RLSeq::analyzeRLFS(
            object = rlr,
            noZ=TRUE,
            ntimes=10,
            quiet = TRUE
        ),
        "RLRanges"
    )
    
    # With no mask
    expect_s4_class(
        RLSeq::analyzeRLFS(
            object = rlr,
            noZ=TRUE,
            useMask = FALSE,
            ntimes=10,
            quiet = TRUE
        ),
        "RLRanges"
    )
    
})
