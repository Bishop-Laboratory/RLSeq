#' Access RLBase samples as RLRanges
#' 
#' Accessor function which returns any sample in RLBase as an RLRanges object for
#' use with RLSeq.
#' 
#' @param acc The sample ID of the RLBase object. \code{rlsample} column in \code{RLHub::rlbase_samples()}.
#' @param quiet If TRUE, messages will be suppressed. Default: FALSE
#' @return An RLRanges object. 
#' @examples
#'
#' rlr <- RLRangesFromRLBase("SRX1070676")
#' 
#' @export
RLRangesFromRLBase <- function(
  acc,
  quiet=FALSE
) {
  rlsamples <- suppressMessages(RLHub::rlbase_samples())
  cvg <- file.path(RLBASE_URL, rlsamples$coverage_s3[rlsamples$rlsample == acc])
  pks <- file.path(RLBASE_URL, rlsamples$peaks_s3[rlsamples$rlsample == acc])
  gen <- rlsamples$genome[rlsamples$rlsample == acc]
  md <- rlsamples$mode[rlsamples$rlsample == acc]
  
  RLRanges(
    peaks = pks,
    coverage = cvg,
    genome = rlsamples$genome[rlsamples$rlsample == acc],
    mode = rlsamples$mode[rlsamples$rlsample == acc],
    condType = rlsamples$condType[rlsamples$rlsample == acc],
    sampleName = rlsamples$name[rlsamples$rlsample == acc],
    quiet = quiet
  )
}
