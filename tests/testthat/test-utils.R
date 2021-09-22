test_that("Access functions work", {
    
    # Load 
    load("rlr.rda")
    
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
