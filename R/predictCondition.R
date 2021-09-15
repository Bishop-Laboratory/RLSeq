#' Predict Condition
#'
#' Uses the results of analyzeRLFS() to predict whether a sample is "Case" 
#' (Normal) or "Control" (RNaseH1-like). This is a useful quality metric when 
#' interpreting R-loop mapping results.
#'
#' @param rlfsRes The results list from running analyzeRLFS().
#' @param ... Internal use only.
#' @return A list containing the results of the fourier analysis and the model prediction.
#' @examples
#' 
#' rlfsRes <- RLSeq::analyzeRLFS(RLSeq::SRX1025890_peaks, genome="hg38")
#' RLSeq::predictCondition(rlfsRes)
#' 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @importFrom stats fft acf
#' @export
predictCondition <- function(rlfsRes, ...) {
  
  # Dots are used to supply custom models for testing purposes only.
  dots <- list(...)
  if (length(dots) == 2) {
    prepFeatures <- dots$prepFeatures
    fftModel <- dots$fftModel
  } else if (length(dots) != 0) {
    stop(
      "Inappropriate arguments supplied: ",
      paste0(
      names(dots)[which(! names(dots) %in% c("prepFeatures", "fftModel"))],
      collapse = ", "
      )
    )
  } else {
    prepFeatures <- RLSeq::prepFeatures
    fftModel <- RLSeq::fftModel
  }
  
  # Get pval
  pval <- rlfsRes$perTestResults$`regioneR::numOverlaps`$pval
  
  # Get Z
  Z <- rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifted.z.scores
  
  # Get edge and center Z values
  shifts <- rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifts
  Zcenter <- Z[which(shifts == 0)]
  Zleft <- Z[which.min(shifts)]
  Zright <- Z[which.max(shifts)]
  
  # compute Fourier transform of Z
  W <- fft(Z)
  
  # compute autocorrelation on Z
  Zacf <- drop(acf(Z, lag.max = length(Z)/1, plot = FALSE, type = "covariance")$acf)
  
  # compute Fourier transform of the autocorrelation 
  Wacf <- fft(Zacf)
  
  # compute first and second moments for Z, Zacf, Re(W), Im(W), Re(Wacf), Im(Wacf) 
  featuresRaw <- data.frame(
    Z1 = mean(Z),
    Z2 = sqrt(sum(Z^2)),
    Zacf1 = mean(Zacf),
    Zacf2 = sqrt(sum(Zacf^2)),
    ReW1 = mean(Re(W)),
    ReW2 = sqrt(sum(Re(W)^2)),
    ImW1 = mean(Im(W)),
    ImW2 = sqrt(sum(Im(W)^2)),
    ReWacf1 = mean(Re(Wacf)),
    ReWacf2 = sqrt(sum(Re(Wacf)^2)),
    ImWacf1 = mean(Im(Wacf)),
    ImWacf2 = sqrt(sum(Im(Wacf)^2))
  )
  
  # Standardize features
  predict.prp <- utils::getFromNamespace("predict.preProcess", "caret")
  features <- predict.prp(prepFeatures, featuresRaw)
  
  # Predict using stacked model
  predict.cs <- utils::getFromNamespace("predict.caretStack", "caretEnsemble")
  pred <- predict.cs(fftModel, features)
  
  # Test each criteria for labeling "Control"
  criteriaOne <- pval < .05
  criteriaTwo <- Zcenter > 0
  criteriaThree <- Zcenter > Zleft & Zcenter > Zright
  criteriaFour <- toupper(pred) == "CASE"
  
  # return results
  list(
    Features = featuresRaw %>%
      tidyr::pivot_longer(dplyr::everything(),
                          names_to = "feature", 
                          values_to = "raw_value") %>%
      dplyr::inner_join(
        tidyr::pivot_longer(
          features, dplyr::everything(),
          names_to = "feature", 
          values_to = "processed_value"
        ), by = "feature"
      ),
    Criteria = list(
      "PVal Significant" = criteriaOne,
      "ZApex > 0" = criteriaTwo,
      "ZApex > ZEdges" = criteriaThree,
      "Predicted 'Case'" = criteriaFour
    ),
    Verdict = ifelse(criteriaOne & criteriaTwo & 
                       criteriaThree & criteriaFour, 
                     "Case", "Control")
  ) %>%
    return()
}
