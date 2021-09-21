test_that(desc = "Feature Enrich returns a tibble", {
    
    # Get small annotation set
    small_anno <- list(
        "Centromeres" = readr::read_csv("https://rlbase-data.s3.amazonaws.com/annotations/hg38/Centromeres.csv.gz", show_col_types = FALSE),
        "SkewR" = readr::read_csv("https://rlbase-data.s3.amazonaws.com/annotations/hg38/SkewR.csv.gz", show_col_types = FALSE)
    )
    
    # Load 
    load("rlr.rda")
    
    # Plot the enrich
    expect_s4_class(
        RLSeq::featureEnrich(
            object = rlr,
            annotations = small_anno, 
            quiet = TRUE
        ),
        "RLRanges"
    )
    
})
