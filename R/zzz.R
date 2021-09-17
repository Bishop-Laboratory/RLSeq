.onLoad <- function(libname, pkgname) {
  # Add constants
  # From https://stackoverflow.com/questions/49056642/how-to-make-variable-available-to-namespace-at-loading-time/49094961
  BASE_UCSC <- 'http://hgdownload.soe.ucsc.edu/goldenPath/'
  
  RLFS_BED_URL <- "https://rlbase-data.s3.amazonaws.com/rlfs-beds/"
  RLFS_BED_S3 <- "s3://rlbase-data/rlfs-beds/"
  RLBASE_BW_URL <- "https://rlbase-data.s3.amazonaws.com/coverage/"
  RLBASE_BW_S3 <- "s3://rlbase-data/coverage/"
  RLBASE_PEAKS_URL <- "https://rlbase-data.s3.amazonaws.com/peaks/"
  RLBASE_PEAKS_S3 <- "s3://rlbase-data/peaks/"
  RLBASE_QUANT_URL <- "https://rlbase-data.s3.amazonaws.com/quant/"
  RLBASE_QUANT_S3 <- "s3://rlbase-data/quant/"
  RLBASE_BAM_STATS_S3 <- "s3://rlbase-data/bam_stats/"
  RLBASE_BAM_STATS_URL <- "https://rlbase-data.s3.amazonaws.com/bam_stats/"
  RLBASE_FASTQ_STATS_S3 <- "s3://rlbase-data/fastq_stats/"
  RLBASE_FASTQ_STATS_URL <- "https://rlbase-data.s3.amazonaws.com/fastq_stats/"
  RLBASE_REPORT_S3 <- "s3://rlbase-data/report/"
  RLBASE_REPORT_URL <- "https://rlbase-data.s3.amazonaws.com/report/"
  RLBASE_DATADUMPS_S3 <- "s3://rlbase-data/datadump/"
  RLBASE_DATADUMPS_URL <- "https://rlbase-data.s3.amazonaws.com/datadump/"
  
  assign("BASE_UCSC", BASE_UCSC, envir = parent.env(environment()))
  assign("RLFS_BED_URL", RLFS_BED_URL, envir = parent.env(environment()))
  assign("RLFS_BED_S3", RLFS_BED_S3, envir = parent.env(environment()))
  assign("RLBASE_BW_URL", RLBASE_BW_URL, envir = parent.env(environment()))
  assign("RLBASE_BW_S3", RLBASE_BW_S3, envir = parent.env(environment()))
  assign("RLBASE_PEAKS_URL", RLBASE_PEAKS_URL, envir = parent.env(environment()))
  assign("RLBASE_PEAKS_S3", RLBASE_PEAKS_S3, envir = parent.env(environment()))
  assign("RLBASE_QUANT_URL", RLBASE_QUANT_URL, envir = parent.env(environment()))
  assign("RLBASE_QUANT_S3", RLBASE_QUANT_S3, envir = parent.env(environment()))
  assign("RLBASE_BAM_STATS_S3", RLBASE_BAM_STATS_S3, envir = parent.env(environment()))
  assign("RLBASE_BAM_STATS_URL", RLBASE_BAM_STATS_URL, envir = parent.env(environment()))
  assign("RLBASE_FASTQ_STATS_S3", RLBASE_FASTQ_STATS_S3, envir = parent.env(environment()))
  assign("RLBASE_FASTQ_STATS_URL", RLBASE_FASTQ_STATS_URL, envir = parent.env(environment()))
  assign("RLBASE_REPORT_S3", RLBASE_REPORT_S3, envir = parent.env(environment()))
  assign("RLBASE_REPORT_URL", RLBASE_REPORT_URL, envir = parent.env(environment()))
  assign("RLBASE_DATADUMPS_S3", RLBASE_DATADUMPS_S3, envir = parent.env(environment()))
  assign("RLBASE_DATADUMPS_URL", RLBASE_DATADUMPS_URL, envir = parent.env(environment()))
  
}
