#' @rdname RLRanges
setClass(
    "RLRanges",
    contains = "GRanges",
    representation = representation(
        peaks = "GRanges"
    )
)
setValidity(
    "RLRanges",
    function(object) {
        # Verification of types
        stopifnot(is.character(object@metadata$mode))
        stopifnot(is.character(object@metadata$coverage))
        stopifnot(is.character(object@metadata$sampleName))

        ## Verification of other options ##

        # Mode
        if (
            !object@metadata$mode %in% auxdata$available_modes$mode &&
                length(object@metadata$mode) > 1
        ) {
            stop("'mode' must be one of auxdata$available_modes$mode or empty.")
        }

        # Genome
        if (any(!GenomeInfoDb::genome(object) %in% auxdata$available_genomes)) {
            stop("'genome' must be one of auxdata$available_genomes.")
        }

        # Seqinfo
        if (all(is.na(GenomeInfoDb::seqinfo(object)@seqlengths))) {
            stop(
                "Problem with genome, seqlengths not found. Please supply pea",
                "ks as GRanges with seqinfo (with seqlengths) already included."
            )
        }

        # label
        if (
            !object@metadata$label %in% c("POS", "NEG") &&
                object@metadata$label != ""
        ) {
            stop(
                "'label' must be one of 'POS'",
                " or 'NEG'-- or be unspecified."
            )
        }

        # coverage
        coverage <- object@metadata$coverage
        if (
            object@metadata$coverage != "" &&
                (!file.exists(coverage) && !urlExists(coverage))
        ) {
            stop(
                "Coverage could not be found. Content of 'coverage' slot: ",
                coverage, ". Set coverage with coverage(object) <-"
            )
        }
    }
)

setMethod(
    f = "show",
    signature = "RLRanges",
    definition = function(object) {
        # Which results are available?
        res <- object@metadata$results
        sn <- methods::slotNames(res)
        fld <- sn[
            vapply(
                X = sn,
                FUN = function(x) {
                    length(methods::slot(res, x)) > 1
                },
                FUN.VALUE = logical(1)
            )
        ]
        if (!length(fld)) fld <- "None"

        GenomicRanges::show(GenomicRanges::GRanges(object))
        cat(
            paste0("\n", object@metadata$sampleName, ":"),
            "\n  Mode:", object@metadata$mode,
            "\n  Genome:", GenomeInfoDb::genome(object)[1],
            "\n  Label:", object@metadata$label
        )
        cat(
            "\n\nRLSeq Results Available:",
            "\n ", paste0(fld, collapse = ", "), "\n"
        )
        cat(
            ifelse(
                "predictRes" %in% fld,
                paste0("\nprediction: ", res@predictRes$prediction, "\n\n"),
                "\n"
            )
        )
    }
)


#' Construct RLRanges Dataset
#'
#' `RLRanges` is a subclass of `GRanges`, which stores R-loop peaks
#' and metadata about the R-loop-mapping experiment, along with results from
#' the analyses in [RLSeq].
#'
#' @param peaks Path/URL to peak file or a GRanges object. This file should be
#' in "broadPeak" format if possible. If not, then `qcol` should be specified.
#' @param coverage Path/URL to the coresponding coverage file
#' (in "bigWig" format). If not supplied, correlation tests will be skipped.
#' @param genome UCSC genome ID. Acceptable types are listed in
#' [auxdata] (`available_genomes` entry).
#' @param mode Type of R-loop mapping from which peaks and coverage were
#' derived. Acceptable types are listed in
#' [auxdata] (`available_modes` entry). Can be unspecified.
#' @param label "POS" (positive R-loop-mapping sample; e.g., DRIP-Seq S9.6 -RNH1)
#' or "NEG" (negative control sample; e.g., DRIP-Seq S9.6 +RNH1 or Input).
#' Can be unspecified.
#' @param sampleName A unique name for identifying this sample.
#' Can be unspecified.
#' @param qcol The name of the metadata column which contains the score or
#' significance of each peak. For broadPeak (preferred), this is the
#' qvalue (column 11 after accounting for extra columns created during
#' peakset building).
#' If not specified, the last column will be chosen by default. **NOTE**:
#' if supplying narrowPeak form peaks, the last column will NOT be appropriate
#' and QCol should be specified as `11`.
#' If FALSE or if no metadata columns exist, it will be left blank and some
#' operations in `report()` will not fully run.
#' @return An object of class `RLRanges`. These objects are an extension of
#' `GRanges` with the addition of sample metadata entries and [RLResults].
#' @aliases RLRanges RLRanges-class
#' @rdname RLRanges
#' @examples
#'
#' # Example dataset
#' rlbase <- "https://rlbase-data.s3.amazonaws.com"
#' cvg <- file.path(rlbase, "coverage", "SRX7671349_hg38.bw")
#' pks <- system.file("extdata", "SRX7671349_hg38.broadPeak", package = "RLSeq")
#'
#' # Get RLRanges object
#' rlr <- RLRanges(pks,
#'     coverage = cvg, genome = "hg38", label = "NEG",
#'     mode = "RDIP", sampleName = "RDIP-Seq +RNH1", qcol = 9
#' )
#'
#' @export
RLRanges <- function(peaks = GenomicRanges::GRanges(),
    coverage = character(1),
    genome = character(1),
    mode = character(1),
    label = character(1),
    sampleName = "User-selected sample",
    qcol = NULL) {

    # Obtain GRanges -- works even if file-path given
    peaks <- regioneR::toGRanges(peaks)

    # Set the qval column
    md <- as.data.frame(methods::slot(peaks, "elementMetadata"))
    if (ncol(md) != 0 & !isFALSE(qcol)) {
        if (is.null(qcol)) {
            qcol <- colnames(md)[ncol(md)]
        } else if (is.numeric(qcol)) {
            qcol <- colnames(as.data.frame(peaks))[qcol]
            message(
                "Column ", qcol, " out of ",
                length(colnames(as.data.frame(peaks))),
                " chosen as 'qcol'. Are you sure this is correct? View",
                " available columns with `colnames(as.data.frame(rlranges))`"
            )
        }
        colnames(
            methods::slot(peaks, "elementMetadata")
        )[colnames(md) == qcol] <- "qval"
    } else {
        warning(
            "No qcol (padjusted) column available...",
            " This will limit some operations in report()."
        )
    }

    # Add in genome info
    GenomeInfoDb::seqlevelsStyle(peaks) <- "UCSC"
    GenomeInfoDb::genome(peaks) <- genome

    # Add in seq info if not already available
    if (any(is.na(GenomeInfoDb::seqinfo(peaks)@seqlengths))) {
        si <- GenomeInfoDb::getChromInfoFromUCSC(genome, as.Seqinfo = TRUE)
        GenomeInfoDb::seqlevels(
            peaks,
            pruning.mode = "coarse"
        ) <- GenomeInfoDb::seqlevels(si)
        GenomeInfoDb::seqinfo(peaks) <- si

        # Trim out-of-bounds ranges
        peaks <- GenomicRanges::trim(peaks)
    }

    # Normalize coverage path
    if (coverage != "") {
        if (!urlExists(coverage) && file.exists(coverage)) {
            # POS: It is a file, which exists. Absolute path.
            coverage <- file.path(
                normalizePath(dirname(coverage)),
                basename(coverage)
            )
        } else if (!urlExists(coverage)) {
            stop(
                "Coverage could not be found. Content of 'coverage': ",
                coverage
            )
        }
    }

    # Add new data to the metadata
    peaks@metadata <- c(
        peaks@metadata,
        list(
            mode = mode,
            label = label,
            coverage = coverage,
            sampleName = sampleName,
            results = methods::new("RLResults")
        )
    )

    # Build object
    methods::new(
        "RLRanges",
        peaks
    )
}


#' RLSeq Results
#'
#' Functions for creating and accessing the R-loop results (RL Results). These
#' are a type of object used for holding the results of the tests implemented
#' in RLSeq. They can be accessed using the `rlresult` function.
#'
#' @aliases rlresult RLResults RLResults-class
#' @rdname RLResults
#' @param object [RLRanges] object.
#' @param resultName Name of the result slot to access. See *details*.
#' @return The contents of the requested slot.
#' @details
#'
#' ## Slot descriptions
#'
#' * `featureEnrichment`
#'   - The `tbl` generated from running [featureEnrich].
#'   - The structure and column descriptions are provided in detail within
#'   [RLHub::feat_enrich_samples].
#' * `correlationMat`
#'   - The `matrix` generated from running [corrAnalyze].
#'   - Contains pairwise pearson correlations between all samples in
#'   [RLBase](https://gccri.bishop-lab.uthscsa.edu/rlbase/)
#'   and the supplied RLRanges object.
#' * `rlfsRes`
#'   - The `list` generated from running [analyzeRLFS].
#'   - See [analyzeRLFS] for description of structure.
#' * `noiseAnalysis`
#'   - The `tbl` generated from running [noiseAnalyze].
#' * `txFeatureOverlap`
#'   - The `tbl` generated from running [txFeatureOverlap].
#' * `geneAnnoRes`
#'   - The `tbl` generated from running [geneAnnotation].
#' * `predictRes`
#'   - The `list` generated from running [predictCondition].
#' * `rlRegionRes`
#'   - The `list` generated from running [rlRegionTest].
#'
#' @examples
#'
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' rlresult(rlr, "predictRes")
#' @export
rlresult <- function(object, resultName) {

    # Check class
    stopifnot(methods::is(object, "RLRanges"))

    # Obtain data from slot
    lst <- methods::slot(object@metadata$results, name = resultName)
    if (!length(lst) > 1) {
        stop(resultName, " not found in object.")
    } else {
        return(lst)
    }
}


#' RLResults-class
#'
#' @rdname RLResults
setClass(
    "RLResults",
    slots = c(
        featureEnrichment = "tbl",
        txFeatureOverlap = "tbl",
        correlationMat = "matrix",
        rlfsRes = "list",
        noiseAnalysis = "tbl",
        geneAnnoRes = "tbl",
        predictRes = "list",
        rlRegionRes = "list"
    ),
    prototype = methods::prototype(
        featureEnrichment = dplyr::tibble(),
        txFeatureOverlap = dplyr::tibble(),
        correlationMat = matrix(),
        rlfsRes = list(),
        noiseAnalysis = dplyr::tibble(),
        geneAnnoRes = dplyr::tibble(),
        predictRes = list(),
        rlRegionRes = list()
    )
)
