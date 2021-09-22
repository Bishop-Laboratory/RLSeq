test_that("Access functions work", {

    # Load
    rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))

    # Get RLFS anno
    expect_s4_class(
        RLSeq:::getRLFSAnno(
            object = rlr
        ),
        "GRanges"
    )

    # Get RLFS anno
    expect_s3_class(
        RLSeq:::getChromSizes(
            object = rlr
        ),
        class = "tbl"
    )
})
