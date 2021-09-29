#' Predict Condition
#'
#' Uses the results of analyzeRLFS() to predict whether a sample is "POS"
#' (robust R-loop mapping) or "NEG" (poor R-loop mapping).
#'
#' @param object An RLRanges object with \code{analyzeRLFS()} already run. 
#' Ignored if rlfsRes provided.
#' @param rlfsRes If object not supplied, provide the rlfsRes list which is
#' obtained from \code{rlresult(object, "rlfsRes")}. 
#' @param ... Internal use only.
#' @return An RLRanges object with predictions included.
#' @examples
#'
#' # Example data with analyzeRLFS already run
#' rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))
#'
#' # predict condition
#' rlr <- predictCondition(rlr)
#' 
#' # With rlfsRes
#' predRes <- predictCondition(rlfsRes=rlresult(rlr, "rlfsRes"))
#' 
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
#' @importFrom stats fft acf
#' @import caretEnsemble
#' @export
predictCondition <- function(object, rlfsRes=NULL, ...) {

    # Obtain RLFS-Res if not supplied
    if (is.null(rlfsRes)) {
        rlfsRes <- rlresult(object, resultName = "rlfsRes")
        rlfsGiven <- FALSE
    } else {
        rlfsGiven <- TRUE
    }
    

    # Check for missing packages and stop if found
    pkgs <- vapply(
        X = c("kernlab", "randomForest", "rpart", "MASS"),
        FUN = requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
    if (any(!pkgs)) {
        stop(
            'Packages needed for predictCondition() but not installed: "',
            paste0(names(pkgs)[which(!pkgs)], collapse = '", "'), '"'
        )
    }

    # Dots are used to supply custom models for testing purposes only.
    dots <- list(...)
    if (length(dots) == 2) {
        prepFeatures <- dots$prepFeatures
        fftModel <- dots$fftModel
    } else {
        prepFeatures <- suppressMessages(RLHub::prep_features())
        fftModel <- suppressMessages(RLHub::fft_model())
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
    Zacf <- drop(
        acf(
            Z, lag.max = length(Z) / 1, plot = FALSE,
            type = "covariance"
        )$acf
    )

    # compute Fourier transform of the autocorrelation
    Wacf <- fft(Zacf)

    # compute 1st and 2nd moments for Z, Zacf, Re(W), Im(W), Re(Wacf), Im(Wacf)
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

    # Test each criteria for labeling "POS"
    criteriaOne <- pval < .05
    criteriaTwo <- Zcenter > 0
    criteriaThree <- Zcenter > Zleft & Zcenter > Zright
    criteriaFour <- toupper(pred) == "POS" || toupper(pred) == "CASE"

    # Wrangle features raw
    finalFeatures <- dplyr::bind_rows(features, featuresRaw) %>%
        t() %>%
        as.data.frame()
    colnames(finalFeatures) <- c("raw_value", "processed_value")
    finalFeatures$feature <- rownames(finalFeatures)
    finalFeatures <- finalFeatures %>%
        dplyr::as_tibble() %>%
        dplyr::relocate(.data$feature, .before = .data$raw_value)

    # return results
    reslst <- list(
        Features = finalFeatures,
        Criteria = list(
            "PVal Significant" = criteriaOne,
            "ZApex > 0" = criteriaTwo,
            "ZApex > ZEdges" = criteriaThree,
            "Predicted 'POS'" = criteriaFour
        ),
        prediction = ifelse(criteriaOne & criteriaTwo &
            criteriaThree & criteriaFour,
        "POS", "NEG"
        )
    )

    # Add to RLRanges and result
    if (! rlfsGiven) {
        methods::slot(object@metadata$results, "predictRes") <- reslst
        return(object)
    } else {
        return(reslst)
    }
}
