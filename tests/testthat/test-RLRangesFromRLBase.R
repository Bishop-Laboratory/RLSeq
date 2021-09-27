test_that(desc = "RLRangesFromRLBase works", {
    expect_s4_class(RLRangesFromRLBase("SRX1070676"), "RLRanges")
})
