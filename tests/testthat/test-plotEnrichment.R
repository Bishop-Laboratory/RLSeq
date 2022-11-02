test_that(desc = "Test plotEnrichment", {

    # Load
    rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))

    # Test conditions
    warns <- capture_warnings({
        res <- RLSeq::plotEnrichment(
            object = rlr,
            modes = c("DRIP", "DRIPc", "qDRIP", "sDRIP", "ssDRIP"),
            pred_POS_only = TRUE,
            splitby = "label"
        )
    })
    expect_error(
        RLSeq::plotEnrichment(
            object = rlr,
            splitby = "ASDASD"
        )
    )
    expect_type(
        res,
        "list"
    )
    expect_true(
        length(res) > 0
    )
})
