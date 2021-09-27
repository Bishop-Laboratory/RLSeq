#' Access RLBase samples as RLRanges
#' 
#' Accessor function which returns any sample in RLBase as an RLRanges object for
#' use with RLSeq.
#' 
#' @param acc The sample ID of the RLBase object. \code{rlsample} column in \code{RLHub::rlbase_samples()}.
#' @param rlsamples The tibble provided by \code{RLHub::rlbase_samples()}.
#' Providing these data ahead of time speeds up this operation. Default: NULL.
#' @return An RLRanges object with all results available.
#' @examples
#'
#' rlr <- RLRangesFromRLBase("SRX1070676")
#' 
#' @export
RLRangesFromRLBase <- function(
  acc,
  rlsamples=NULL
) {
  if (is.null(rlsamples)) rlsamples <- suppressMessages(RLHub::rlbase_samples())
  i <- which(rlsamples$rlsample == acc)[1]
  aws.s3::s3readRDS(object = rlsamples$rlranges_rds_s3[i], bucket = RLBASE_S3)
}
