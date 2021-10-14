test_that(desc = "RLSeq report", {
    # Load
    rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))

    # Plot the enrich
    expect_true(
        RLSeq::report(
            object = rlr,
            quiet = TRUE
        )
    )

    # Download coverage locally and try with local file
    coverage <- tempfile(tmpdir = ".", fileext = ".bw")
    download.file(rlr@metadata$coverage, destfile = coverage, quiet = TRUE)
    rlrloc <- rlr
    coverage <- file.path(
        normalizePath(dirname(coverage)),
        basename(coverage)
    )
    rlrloc@metadata$coverage <- coverage
    expect_true(
        RLSeq::report(
            object = rlrloc,
            quiet = TRUE
        )
    )
})
