#' Annotate R-Loops with Genes
#'
#' Annotates RLRanges with entrez ids for overlapping genes. See *details*.
#'
#' @param object An RLRanges object.
#' @param txdb The TxDb or EnsDb object containing gene annotations. If not
#' supplied, annotations will be automatically downloaded from AnnotationHub.
#' See also [GenomicFeatures::TxDb].
#' @return An RLRanges object with gene overlaps included. The results
#' are available via `rlresult(object, "geneAnnoRes")`. The result object
#' is a `tbl` with a mapping of `peak_name` (peak names from `names(object)`)
#' to `gene_id` (entrez gene IDs).
#' @details
#'
#' The `geneAnnotation` function provides a simple procedure for annotating
#' RLRanges with gene IDs by overlap.
#'
#' ### Annotations
#'
#' First. gene annotations are
#' automatically downloaded using [AnnotationHub::query] with the following
#' pattern:
#'
#' ```r
#' AnnotationHub::query(
#'     x = ah,
#'     pattern = c("TxDb", "UCSC", "knownGene", genome)
#' )
#' ```
#'
#' Where `genome` is the UCSC genome id for the RLRanges object. If these
#' annotations are unavailable, they should be provded using the `txdb`
#' parameter. See also [GenomicFeatures::TxDb].
#'
#' ### Overlaps
#'
#' The annotations are subsequently overlapped with the ranges in the
#' supplied RLRanges object using [valr::bed_intersect] and saved in the
#' [RLResults] object as a `tbl` with a mapping of peak names to `gene_id`
#' (entrez gene IDs).
#'
#' @examples
#'
#' # Example RLRanges data
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Perform gene annotation
#' rlr <- geneAnnotation(rlr)
#'
#' # Supply custom TxDb if needed
#' if (GenomeInfoDb::genome(rlr)[1] == "hg19") {
#'     library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'     rlr <- geneAnnotation(rlr, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)
#' }
#' @export
geneAnnotation <- function(object, txdb = NULL) {

    # Get genome from object
    genome <- GenomeInfoDb::genome(object)[1]
    if (!genome %in% c("hg38", "mm10") && is.null(txdb)) {
        stop(
            "No gene annotations available for object genome, ", genome,
            ". Provide TxDb or EnsDb object."
        )
    }

    # If no TxDb provided, obtain from annotationhub
    if (is.null(txdb)) {
        ah <- AnnotationHub::AnnotationHub()
        ahDb <- AnnotationHub::query(
            x = ah,
            pattern = c("TxDb", "UCSC", "knownGene", genome)
        )
        txdb <- ah[[names(which.max(ahDb@.db_uid)), verbose = FALSE]]
    }

    # Get the ensembl genes and conver to UCSC style
    edb <- GenomicFeatures::genes(txdb)
    GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"

    # Wrangle EnsDb to tibble
    annoData <- edb %>%
        as.data.frame() %>%
        dplyr::select(
            chrom = .data$seqnames, .data$start,
            .data$end, .data$strand, .data$gene_id
        ) %>%
        dplyr::as_tibble() %>%
        dplyr::distinct(.data$gene_id, .keep_all = TRUE) %>%
        dplyr::mutate(chrom = as.character(.data$chrom))

    # Wrangle peaks to tibble
    pkNames <- names(object)
    peaksIntersect <- object %>%
        as.data.frame() %>%
        dplyr::as_tibble() %>%
        dplyr::select(
            chrom = .data$seqnames, .data$start,
            .data$end, .data$width
        ) %>%
        dplyr::mutate(
            chrom = as.character(.data$chrom),
            pkName = {{ pkNames }}
        )

    # Intersect
    anno <- valr::bed_intersect(
        peaksIntersect,
        annoData,
        suffix = c("__userPeaks", "__Gene")
    )

    # Clean
    pk_to_gene <- dplyr::select(
        anno,
        peak_name = .data$pkName__userPeaks,
        gene_id = .data$gene_id__Gene
    ) %>%
        dplyr::distinct()

    # Return to object
    methods::slot(object@metadata$results, name = "geneAnnoRes") <- pk_to_gene

    return(object)
}
