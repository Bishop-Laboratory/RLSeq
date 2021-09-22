.onLoad <- function(libname, pkgname) {
    # Links
    BASE_UCSC <- "http://hgdownload.soe.ucsc.edu/goldenPath"
    RLBASE_URL <- "https://rlbase-data.s3.amazonaws.com"
    RLBASE_S3 <- "s3://rlbase-data"

    # Assign
    assign("BASE_UCSC", BASE_UCSC, envir = parent.env(environment()))
    assign("RLBASE_URL", RLBASE_URL, envir = parent.env(environment()))
    assign("RLBASE_S3", RLBASE_S3, envir = parent.env(environment()))
    
    # Globals
    utils::globalVariables(".")
}
