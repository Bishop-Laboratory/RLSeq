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
    expect_match(warns[1], regexp = "User-supplied sample test value is NA for .*")
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
