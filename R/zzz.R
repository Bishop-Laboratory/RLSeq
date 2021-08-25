.onLoad <- function(libname, pkgname) {
  # Add constants
  # From https://stackoverflow.com/questions/49056642/how-to-make-variable-available-to-namespace-at-loading-time/49094961
  BASE_UCSC <- 'http://hgdownload.soe.ucsc.edu/goldenPath/'
  RLFS_BED_URL <- "https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/"
  assign("BASE_UCSC", BASE_UCSC, envir = parent.env(environment()))
  assign("RLFS_BED_URL", RLFS_BED_URL, envir = parent.env(environment()))
}
