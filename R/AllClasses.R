#' @rdname RLRanges
#' @importClassesFrom GenomicRanges GenomicRanges
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
        if (!object@metadata$mode %in% auxdata$available_modes$mode &&
            length(object@metadata$mode) > 1) {
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
        if (!object@metadata$label %in% c("POS", "NEG") &&
            object@metadata$label != "") {
            stop(
                "'label' must be one of 'POS'",
                " or 'NEG'-- or be unspecified."
            )
        }

        # coverage
        coverage <- object@metadata$coverage
        if (object@metadata$coverage != "" &&
            (!file.exists(coverage) && !urlExists(coverage))) {
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

        sgr <- utils::getFromNamespace("show_GenomicRanges", "GenomicRanges")
        sgr(object)
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
#' \code{RLRanges} is a subclass of \code{GRanges}, which stores R-loop peaks
#' and metadata about the R-loop-mapping experiment.
#'
#' @param peaks Path/URL to peak file or a GRanges object.
#' @param coverage Path/URL to bigWig file. If not supplied, correlation tests
#'  will be skipped.
#' @param genome UCSC genome ID. Acceptable types are listed in 
#' auxdata$available_genomes.
#' @param mode Type of R-loop mapping from which peaks and coverage were
#' derived. Acceptable types are listed in RLSeq::auxdata$available_modes$mode.
#'  Can
#' be unspecified.
#' @param label "POS" (e.g., S9.6 -RNH1) or "NEG" (e.g., S9.6 +RNH1 or Input).
#'  Can be unspecified.
#' @param sampleName A unique name for identifying this sample. 
#' Can be unspecified.
#' @param qcol The name of the metadata column which contains the adjusted p
#'  value. 
#' If not specified, the last column will be chosen (standard for .broadPeak 
#' files).
#' If FALSE or if no metadata columns exist, it will be left blank and some 
#' operations in report() will not fully run.
#' @param quiet If TRUE, messages and warnings are suppressed. Default: FALSE.
#' @return An object of class "RLRanges".
#' @examples
#'
#' # Example dataset
#' rlbase <- "https://rlbase-data.s3.amazonaws.com"
#' pks <- file.path(rlbase, "peaks", "SRX7671349_hg38.broadPeak")
#' cvg <- file.path(rlbase, "coverage", "SRX7671349_hg38.bw")
#'
#' # Get RLRanges object
#' rlr <- RLRanges(pks,
#'     coverage = cvg, genome = "hg38", label="NEG",
#'     mode = "RDIP", sampleName = "RDIP-Seq +RNH1"
#' )
#' @export
RLRanges <- function(peaks = GenomicRanges::GRanges(),
    coverage = character(1),
    genome = character(1),
    mode = character(1),
    label = character(1),
    sampleName = "User-selected sample",
    qcol = NULL,
    quiet = FALSE) {

    # Obtain GRanges -- works even if file-path given
    peaks <- regioneR::toGRanges(peaks)
    
    # Set the qval column
    md <- as.data.frame(methods::slot(peaks, "elementMetadata"))
    if (ncol(md) != 0 & ! isFALSE(qcol)) {
        if (is.null(qcol)) {
            qcol <- colnames(md)[ncol(md)]
        } 
        colnames(
            methods::slot(peaks, "elementMetadata")
        )[colnames(md) == qcol] <- "qval"
    } else {
        if (! quiet) {
            warning("No qcol (padjusted) column available...",
                    " This will limit some operations in report().")
        }
    }

    # Add in genome info
    GenomeInfoDb::seqlevelsStyle(peaks) <- "UCSC"
    GenomeInfoDb::genome(peaks) <- genome

    # Add in seq info if not already available
    if (any(is.na(GenomeInfoDb::seqinfo(peaks)@seqlengths))) {
        si <- GenomeInfoDb::getChromInfoFromUCSC(genome, as.Seqinfo = TRUE)
        if (quiet) {
            suppressWarnings(
                GenomeInfoDb::seqlevels(peaks) <- GenomeInfoDb::seqlevels(si)
            )
            suppressWarnings(GenomeInfoDb::seqinfo(peaks) <- si)
        } else {
            GenomeInfoDb::seqlevels(peaks) <- GenomeInfoDb::seqlevels(si)
            GenomeInfoDb::seqinfo(peaks) <- si
        }

        # Trim out-of-bounds ranges
        peaks <- GenomicRanges::trim(peaks)
    }
    
    # Normalize coverage path
    if (coverage != "") {
        if (! urlExists(coverage) && file.exists(coverage)) {
            # POS: It is a file, which exists. Absolute path.
            coverage <- file.path(normalizePath(dirname(coverage)), 
                                  basename(coverage))
        } else if (! urlExists(coverage)) {
            stop("Coverage could not be found. Content of 'coverage': ",
                coverage)
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


#' Result accessor function
#'
#' @param object RLRanges object.
#' @param resultName Name of the result slot to access. See details.
#' @return The contents of the requested slot.
#' @details 
#' 
#' \strong{"featureEnrichment"} The \code{tbl} generated from running the 
#' \code{featureEnrich()} function.
#' 
#' \strong{"correlationMat"}. The \code{matrix} generated from running the 
#' \code{corrAnalyze()} function.
#' 
#' \strong{"rlfsRes"}. The \code{list} generated from running the 
#' \code{analyzeRLFS()} function.
#' 
#' \strong{"geneAnnoRes"}. The \code{tbl} generated from running the 
#' \code{geneAnnotation()} function.
#' 
#' \strong{"predictRes"}. The \code{list} generated from running the 
#' \code{predictCondition()} function.
#' 
#' \strong{"rlRegionRes"}. The \code{list} generated from running the 
#' \code{rlRegionTest()} function.
#' 
#' @examples 
#' 
#' rlr <- readRDS(system.file("ext-data", "rlrsmall.rds", package = "RLSeq"))
#' 
#' rlresult(rlr, "predictRes")
#' 
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


#' Construct RLResults Object
#'
#' \code{RLResults} is a class used internally for storing the 
#' results of running
#' \code{RLSeq} operations on an \code{RLRanges} object.
#' @slot featureEnrichment The \code{tbl} generated from running the 
#' \code{featureEnrich()} function.
#' @slot correlationMat The \code{matrix} generated from running the 
#' \code{corrAnalyze()} function.
#' @slot rlfsRes The \code{list} generated from running the 
#' \code{analyzeRLFS()} function.
#' @slot geneAnnoRes The \code{tbl} generated from running the 
#' \code{geneAnnotation()} function.
#' @slot predictRes The \code{list} generated from running the 
#' \code{predictCondition()} function.
#' @slot rlRegionRes The \code{list} generated from running the 
#' \code{rlRegionTest()} function.
setClass(
    "RLResults",
    slots = c(
        featureEnrichment = "tbl",
        correlationMat = "matrix",
        rlfsRes = "list",
        geneAnnoRes = "tbl",
        predictRes = "list",
        rlRegionRes = "list"
    ),
    prototype = methods::prototype(
        featureEnrichment = dplyr::tibble(),
        correlationMat = matrix(),
        rlfsRes = list(),
        geneAnnoRes = dplyr::tibble(),
        predictRes = list(),
        rlRegionRes = list()
    )
)
