.onLoad <- function(libname, pkgname) {
  # Add constants
  # From https://stackoverflow.com/questions/49056642/how-to-make-variable-available-to-namespace-at-loading-time/49094961
  BASE_UCSC <- 'http://hgdownload.soe.ucsc.edu/goldenPath/'
  RLFS_BED_URL <- "https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/"
  RMAP_BW_URL <- "https://rmapdb-data.s3.us-east-2.amazonaws.com/bigwigs/rseq-coverage-unstranded/"
  RMAP_PEAKS_URL <- "https://rmapdb-data.s3.us-east-2.amazonaws.com/macs2-peaks-unstranded/"
  RMAP_PEAKS_S3 <- "s3://rmapdb-data/macs2-peaks-unstranded/"
  assign("BASE_UCSC", BASE_UCSC, envir = parent.env(environment()))
  assign("RLFS_BED_URL", RLFS_BED_URL, envir = parent.env(environment()))
  assign("RMAP_BW_URL", RMAP_BW_URL, envir = parent.env(environment()))
  assign("RMAP_PEAKS_URL", RMAP_PEAKS_URL, envir = parent.env(environment()))
  assign("RMAP_PEAKS_S3", RMAP_PEAKS_S3, envir = parent.env(environment()))
}
