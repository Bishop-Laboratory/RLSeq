#' Access RLBase samples as RLRanges
#'
#' Accessor function which returns any sample in
#' [RLBase](https://gccri.bishop-lab.uthscsa.edu/rlbase/) as an 
#' [RLRanges] object for use with RLSeq. For a full list of available 
#' samples, see [RLHub::rlbase_samples].
#'
#' @param acc The sample ID of the RLBase object. See the `rlsample` column
#' in [RLHub::rlbase_samples].
#' @param rlsamples The tibble provided by [RLHub::rlbase_samples].
#' Providing these data ahead of time speeds up this operation. Default: NULL.
#' @return An RLRanges object with all results available.
#' @examples
#' 
#' rlr <- RLRangesFromRLBase("SRX1070676")
#' 
#' @export
RLRangesFromRLBase <- function(acc,
    rlsamples = NULL) {
    if (is.null(rlsamples)) rlsamples <- RLHub::rlbase_samples(quiet = TRUE)
    i <- which(rlsamples$rlsample == acc)[1]
    aws.s3::s3readRDS(object = rlsamples$rlranges_rds_s3[i], bucket = RLBASE_S3)
}
