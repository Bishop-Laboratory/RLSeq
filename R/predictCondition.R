#' Predict Condition
#'
#' Uses the results of [analyzeRLFS] to predict whether a sample is "POS"
#' (robust R-loop mapping) or "NEG" (poor R-loop mapping). See *details*.
#'
#' @param object An RLRanges object with [analyzeRLFS] already run.
#' Ignored if `rlfsRes` is provided.
#' @param rlfsRes If object not supplied, provide the rlfsRes list which is
#' obtained from `rlresult(object, "rlfsRes")`.
#' @param ... Internal use only.
#' @details
#'
#' Following R-loop forming sequences (RLFS) analysis, the quality model
#' (see [RLHub::models]) is implemented for predicting the sample condition in
#' coordination with other results from [analyzeRLFS].
#' A prediction of “POS” indicates robust R-loop mapping, whereas “NEG”
#' indicates poor R-loop mapping. The succeeding sections describe this
#' process in greater detail.
#'
#' ### Application of binary classification model
#'
#' First, the binary classifier is applied, yielding a preliminary
#' prediction of quality. This is accomplished via the following
#' steps:
#'
#' 1. Calculate the Fourier transform of the Z-score distribution
#' (see [analyzeRLFS]).
#' 2. Reduce the dimensions to the engineered feature set (see table below).
#' 3. Apply the preprocessing model (see [RLHub::models]) to normalize these
#' features
#' 4. Apply the classifier (see [RLHub::models]) to render a quality prediction.
#'
#' #### Engineered feature set
#'
#' Abbreviations: Z, Z-score distribution; ACF, autocorrelation function;
#' FT, Fourier Transform.
#'
#' |feature |description                              |
#' |:-------|:----------------------------------------|
#' |Z1      |mean of Z                                |
#' |Z2      |variance of Z                            |
#' |Zacf1   |mean of Z ACF                            |
#' |Zacf2   |variance of Z ACF                        |
#' |ReW1    |mean of FT of Z (real part)              |
#' |ReW2    |variance of FT of Z (real part)          |
#' |ImW1    |mean of FT of Z (imaginary part)         |
#' |ImW2    |variance of FT of Z (imaginary part)     |
#' |ReWacf1 |mean of FT of Z ACF (real part)          |
#' |ReWacf2 |variance of FT of Z ACF (real part)      |
#' |ImWacf1 |mean of FT of Z ACF (imaginary part)     |
#' |ImWacf2 |variance of FT of Z ACF (imaginary part) |
#'
#' ### Final quality prediction
#'
#' The results from the binary classifier are combined with other results from
#' [analyzeRLFS] to yield a final prediction. To yield a prediction of “POS”
#' all the following must be `TRUE`:
#'
#' 1. The RLFS Permutation test P value is significant (p < .05). Stored as
#' `PVal Significant` in the results object.
#' 2. The Z-score distribution at 0bp is > 0. Stored as `ZApex > 0` in
#' the results object.
#' 3. The Z-score distribution at 0bp is > the start and the end. Sored as
#' `ZApex > ZEdges` in the results object.
#' 4. binary The classifier predicts a label of “POS”. Stored as
#' `Predicted 'POS'` in the results object.
#'
#' @return An RLRanges object with predictions accessible via
#' `rlresult(object, "predictRes")`.
#'
#' ### Structure
#'
#' The results object is a named `list` of the structure:
#'
#' * `Features`
#'   - A `tbl` with three columns that describe the engineered features used
#'   for prediction:
#'     * `feature`: the name of the feature (see *details*).
#'     * `raw_value`: The raw value of that feature in the supplied object.
#'     * `processed_value`: The normalized value of that feature after
#'     preprocessing (see *details*).
#' * Criteria
#'   - The four criteria which must all be `TRUE` to render a
#'   prediction of "POS" (see *details*).
#' * prediction
#'   - The final prediction. "POS" indicates robust R-loop mapping, "NEG"
#'   indicates poor R-loop mapping.
#'
#' @examples
#'
#' # Example data with analyzeRLFS already run
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # predict condition
#' rlr <- predictCondition(rlr)
#'
#' # With rlfsRes
#' predRes <- predictCondition(rlfsRes = rlresult(rlr, "rlfsRes"))
#' @export
predictCondition <- function(object, rlfsRes = NULL, ...) {

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
        pks <- paste0(names(pkgs)[which(!pkgs)], collapse = '", "')
        stop(
            'Packages needed for predictCondition() but not installed: "',
            pks, '"'
        )
    }

    # Dots are used to supply custom models for testing purposes only.
    dots <- list(...)
    if (length(dots) == 2) {
        prepFeatures <- dots$prepFeatures
        fftModel <- dots$fftModel
    } else {
        prepFeatures <- RLHub::prep_features(quiet = TRUE)
        fftModel <- RLHub::fft_model(quiet = TRUE)
    }

    # Get pval
    pval <- rlfsRes$perTestResults$`regioneR::numOverlaps`$pval

    # Get Z
    Z <- rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifted.z.scores
    if (any(is.infinite(Z))) {
        stop(
            "Found Inf values in Z score distribution.",
            " Set 'ntimes' higher in analyzeRLFS"
        )
    }

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
            Z,
            lag.max = length(Z) / 1, plot = FALSE,
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
    features <- stats::predict(prepFeatures, featuresRaw)

    # Predict using stacked model
    pred <- stats::predict(fftModel, features)

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
    if (!rlfsGiven) {
        methods::slot(object@metadata$results, "predictRes") <- reslst
        return(object)
    } else {
        return(reslst)
    }
}
