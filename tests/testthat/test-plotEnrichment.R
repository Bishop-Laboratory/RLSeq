test_that(desc = "Test plotEnrichment", {

    # Load
    rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))

    # Test conditions
    expect_type(
        RLSeq::plotEnrichment(
            object = rlr,
            modes = c("DRIP", "DRIPc", "qDRIP", "sDRIP", "ssDRIP"),
            onlyCase = TRUE,
            splitby = "condType"
        ),
        "list"
    )
    expect_type(
        RLSeq::plotEnrichment(
            object = rlr,
            modes = c("DRIP", "DRIPc", "qDRIP", "sDRIP", "ssDRIP"),
            onlyCase = FALSE,
            splitby = "verdict",
            returnData = TRUE
        ),
        "list"
    )
    expect_error(
        RLSeq::plotEnrichment(
            object = rlr,
            splitby = "ASDASD"
        )
    )
})
