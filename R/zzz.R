.onLoad <- function(libname, pkgname) {
    # Links
    BASE_UCSC <- "http://hgdownload.soe.ucsc.edu/goldenPath"
    RLBASE_URL <- "https://rlbase-data.s3.amazonaws.com"
    RLBASE_S3 <- "s3://rlbase-data"

    # Data
    # This is my workaround for the no lazyload NOTE in bioccheck
    auxdata <- readRDS(
        system.file("int-data", "auxdata.rds", package = "RLSeq")
    )
    genomeMasks <- readRDS(
        system.file("int-data", "genomeMasks.rds", package = "RLSeq")
    )
    available_genomes <- readRDS(
        system.file("int-data", "available_genomes.rds", package = "RLSeq")
    )

    # Assign
    assign("BASE_UCSC", BASE_UCSC, envir = parent.env(environment()))
    assign("RLBASE_URL", RLBASE_URL, envir = parent.env(environment()))
    assign("RLBASE_S3", RLBASE_S3, envir = parent.env(environment()))
    assign("auxdata", auxdata, envir = parent.env(environment()))
    assign("genomeMasks", genomeMasks, envir = parent.env(environment()))
    assign(
        "available_genomes", available_genomes,
        envir = parent.env(environment())
    )

    # Globals
    utils::globalVariables(".")
}
