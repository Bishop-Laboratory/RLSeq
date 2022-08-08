#' RLSeq
#'
#' Executes the RLSeq analysis workflow.
#'
#' @param object An [RLRanges] object.
#' @param quiet If `TRUE`, messages are suppressed. Default: `FALSE`.
#' @param skip Analysis steps to skip.
#' Default: `NULL`. See *details* for options.
#' @param ... Arguments passed to [analyzeRLFS].
#' @return An [RLRanges] object with results available (see [rlresult]).
#' @details
#'
#' The `RLSeq()` function does all of the following by default:
#'
#' 1. **RLFS Perm Test**. Runs the [analyzeRLFS] function to test the
#'   enrichment of user-supplied ranges within R-loop-forming sequences.
#'   *Cannot be skipped.*
#' 2. **Predict Condition**. Runs the [predictCondition] function to
#'   predict whether the user-supplied sample robustly maps R-loops or not.
#'   *Cannot be skipped.*
#' 3. **Noise analysis**. Runs the [noiseAnalyze] function to analyze the
#'   signal-noise distribution within the supplied dataset. Skip with
#'   `skip="noiseAnalyze"`.
#' 4. **Feature enrichment test**. Runs the [featureEnrich] function to
#'   test the enrichment of user-supplied ranges within R-loop-relevant
#'   genomic features. Skip with `skip="featureEnrich"`.
#' 5. **Transcript Feature Overlap**. Runs the [txFeatureOverlap] function to
#'   get the overlap of transcript features and user-supplied peaks.
#' 6. **Correlation Analysis**. Runs the [corrAnalyze] function to test
#'   the correlation of user-supplied R-loop signal with other samples in
#'   RLBase around "gold-standard" R-loop regions.
#'   Skip with `skip="corrAnalyze"`.
#' 7. **Gene annotation**. Runs the [geneAnnotation] function to find overlap
#'   of genes with the user-supplied ranges. Skip with `skip="geneAnnotation"`.
#' 8. **R-loop Region Analysis**. Runs the [rlRegionTest] function to find
#'   the overlap of user-supplied ranges with consensus R-loop sites
#'   (RL-Regions). Skip with `skip="rlRegionTest"`.
#'
#' @examples
#'
#' # Example RLRanges
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Run RLSeq
#' # `useMask=FALSE`, `ntime=10`, and `skip=` for demonstration purposes here.
#' rlr <- RLSeq(
#'     rlr,
#'     useMask = FALSE, ntimes = 10,
#'     skip = c(
#'         "featureEnrich", "corrAnalyze", "geneAnnotation", "rlRegionTest"
#'     )
#' )
#'
#' @importFrom dplyr %>% .data bind_rows tibble relocate as_tibble mutate select
#' @importFrom dplyr filter distinct sample_n bind_cols row_number desc arrange
#' @importFrom stats fft acf predict ks.test
#' @importFrom RLHub prep_features fft_model rlregions_meta
#' @importFrom RLHub annots_primary_hg38 annots_primary_mm10
#' @importFrom RLHub annots_full_hg38 annots_full_mm10
#' @import caretEnsemble
#' @importClassesFrom GenomicRanges GenomicRanges
#' @importFrom GenomeInfoDb genome seqinfo seqlevelsStyle getChromInfoFromUCSC
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges trim show strand reduce GRanges
#' @importFrom methods slotNames slot new is prototype
#' @importFrom regioneR toGRanges randomizeRegions
#' @importFrom aws.s3 s3readRDS
#' @importFrom callr r
#' @importFrom valr bed_intersect bed_fisher bed_merge bed_reldist
#' @importFrom valr gr_to_bed bed_map
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom GenomicFeatures genes
#' @export
RLSeq <- function(object, quiet = FALSE, skip = NULL, ...) {
    if (!quiet) message("[1/8] RLFS Perm Test")
    object <- analyzeRLFS(object, quiet = TRUE, ...)

    if (!quiet) message("[2/8] Predict Condition")
    object <- predictCondition(object)

    if (!quiet) message("[3/8] Noise Analysis")
    if (object@metadata$coverage != "" & !"noiseAnalyze" %in% skip) {
        if (GenomeInfoDb::genome(object)[1] %in% c("hg38", "hg19", "mm10")) {
            if (.Platform$OS.type != "windows") {
                object <- noiseAnalyze(object)
            } else {
                if (!quiet) {
                    warning(
                        "noiseAnalyze does not work on windows OS. Please",
                        " run noiseAnalyze() directly with ",
                        "`force=TRUE` to override."
                    )
                }
            }
        } else {
            warning(
                "Only 'hg38', 'hg19', and 'mm10' genome ranges are available",
                " for noise analysis. Specify custom ranges using the
                'windows' option in noiseAnalyze(). Skipping."
            )
        }
    } else {
        message("No coverage provided... skipping.")
    }

    if (!quiet) message("[4/8] Feature Enrichment Test")
    if (!"featureEnrich" %in% skip) {
        if (GenomeInfoDb::genome(object)[1] %in% c("hg38", "mm10")) {
            object <- featureEnrich(object, quiet = TRUE)
        } else {
            warning(
                "RLSeq only contains built-in annotations for 'hg38' and",
                "'mm10' genomes. Please lift-over or run featureEnrich() ",
                "separately with custom annotations. Skipping."
            )
        }
    }

    if (!quiet) message("[5/8] Transcript Feature Overlap Analysis")
    if (!"txFeatureOverlap" %in% skip) {
        if (GenomeInfoDb::genome(object)[1] %in% c("hg38", "mm10")) {
            object <- txFeatureOverlap(object, quiet = TRUE)
        } else {
            warning(
                "RLSeq only contains built-in annotations for 'hg38' and",
                "'mm10' genomes. Please lift-over or run featureEnrich() ",
                "separately with custom annotations. Skipping."
            )
        }
    }

    if (!quiet) message("[6/8] Correlation Analysis")
    if (object@metadata$coverage != "" & !"corrAnalyze" %in% skip) {
        if (GenomeInfoDb::genome(object)[1] == "hg38") {
            if (.Platform$OS.type != "windows") {
                object <- corrAnalyze(object)
            } else {
                if (!quiet) {
                    warning(
                        "corrAnalyze does not work on windows OS. Please",
                        " run corrAnalyze() directly with ",
                        "`force=TRUE` to override."
                    )
                }
            }
        } else {
            warning(
                "Only 'hg38' genome ranges are available",
                " for correlation analysis. Skipping."
            )
        }
    } else {
        message("No coverage provided... skipping.")
    }

    if (!quiet) message("[7/8] Gene Annotation")
    if (!"geneAnnotation" %in% skip) {
        objectanno <- try(geneAnnotation(object), silent = TRUE)
        if ("try-error" %in% class(objectanno)) {
            warning(objectanno)
        } else {
            object <- objectanno
        }
    }

    if (!quiet) message("[8/8] R-loop Region Analysis")
    if (!"rlRegionTest" %in% skip) {
        if (GenomeInfoDb::genome(object)[1] == "hg38") {
            object <- rlRegionTest(object)
        } else {
            warning(
                "Only 'hg38' genome ranges",
                " are available for RL-region analysis. Skipping."
            )
        }
    }
    return(object)
}
