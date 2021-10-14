#' Plot RLFS analysis results
#'
#' Plots the results of the R-loop-forming sequences (RLFS) analysis. The plot is a 
#' metaplot of the Z score distribution around RLFS with the p value from 
#' permutation testing annotated. See also [analyzeRLFS].
#' 
#' @param object An RLRanges object with [analyzeRLFS] already run.
#' Alternatively, can be the RLFS results object from an RLRanges (from
#' `rlresult(object, "rlfsRes")`).
#' @param plotName A Sample name used for plotting.
#' If blank, the RLRanges `sampleName` metadata entry is used (see [RLRanges]).
#' @param fft If TRUE, the Fourier transform of the Z-score is plotted instead.
#' Default: FALSE.
#' @param ... Additional parameters passed to [ggplot2::ggplot].
#' @return A ggplot object. See also [ggplot2::ggplot].
#' @examples
#'
#' # Example RLRanges dataset with analyzeRLFS() already run.
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Plot RLFS res
#' plotRLFSRes(rlr)
#'
#' # Plot the Fourier transform instead
#' plotRLFSRes(rlr, fft = TRUE)
#' @export
plotRLFSRes <- function(object,
    plotName = NULL,
    fft = FALSE,
    ...) {

    # Obtain RLFS-Res
    if (methods::is(object, "RLRanges")) {
        rlfsRes <- rlresult(object, resultName = "rlfsRes")
        if (is.null(plotName)) plotName <- object@metadata$sampleName
    } else if (
        methods::is(object, "list") &&
            "perTestResults" %in% names(object)
    ) {
        rlfsRes <- object
        if (is.null(plotName)) plotName <- "User-supplied sample"
    }

    # Obtain the plot data
    zs <- rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifted.z.scores
    pltdat <- dplyr::tibble(
        "zscore" = if (fft) Re(stats::fft(zs)) else zs,
        "shift" = rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifts
    )
    pval <- rlfsRes$perTestResults$`regioneR::numOverlaps`$pval

    # Control labels based on FFT
    plttitle <- ifelse(
        fft, "Fourier-transformed Z-Score around RLFS", "Z-score around RLFS"
    )
    pltylab <- ifelse(fft, "FFT of Z-score", "Peak Enrichment (Z-Score)")
    pltxlab <- ifelse(fft, "Frequency", "Distance to RLFS (bp)")
    # Make plot
    pltbase <- ggplot2::ggplot(
        data = pltdat,
        ggplot2::aes_string(
            y = "zscore",
            x = "shift"
        ),
        ...
    )
    if (!fft) {
        pltbase <- pltbase + ggplot2::geom_vline(
            color = "firebrick",
            xintercept = 0,
            linetype = "dashed"
        ) + ggplot2::scale_y_continuous(expand = c(0, 0))
    }
    pltbase +
        ggplot2::geom_line(size = 1) +
        ggplot2::ggtitle(plttitle) +
        ggplot2::labs(caption = paste0("p < ", round(pval, digits = 5))) +
        ggplot2::ylab(pltylab) +
        ggplot2::xlab(pltxlab) +
        ggprism::theme_prism(base_size = 14) +
        ggplot2::labs(subtitle = plotName)
}


#' Plot Correlation Results
#' 
#' Plots a heatmap to visualize the pairwise Pearson correlation matrix
#' generated via [corrAnalyze].
#'
#' @param object An RLRanges with [corrAnalyze] already run.
#' @param returnData If TRUE, plot data is returned instead of plotting.
#' Default: FALSE
#' @param complex If TRUE, [ComplexHeatmap::Heatmap] will be used for plotting.
#' Otherwise, [pheatmap::pheatmap] is used. Default: TRUE
#' @param ... For internal use.
#' @return A plot object or plotting data (if `returnData` is `TRUE`).
#' @examples
#'
#' # Example RLRanges data with corrAnalyze() already run.
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Corr heatmap
#' corrHeatmap(rlr)
#' @export
corrHeatmap <- function(object,
    returnData = FALSE,
    complex = TRUE,
    ...) {

    # Get dots -- get values used by RLBase
    # Provides RLBase with a speed boost by pre-supplying
    # data from memory when plotting. Not intended for regular use.
    dots <- list(...)
    selected <- NULL
    rlsamples <- NULL
    if (length(dots) > 0) {
        selected <- dots$selected
        rlsamples <- dots$rlsamples
    }

    # Get RLBase samples
    if (is.null(rlsamples)) rlsamples <- RLHub::rlbase_samples(quiet = TRUE)

    # Get the correlation matrix
    corrRes <- rlresult(object, resultName = "correlationMat")

    # Get the mode and prediction and label
    prediction <- rlresult(object, resultName = "predictRes")

    # Wrangle the annotation data
    rlsamples <- rlsamples[rlsamples$rlsample != object@metadata$sampleName, ]
    annoCorr <- rlsamples %>%
        dplyr::mutate(Selected = "RLBase") %>%
        dplyr::select(
            .data$rlsample, Mode = .data$mode,
            Label = .data$label,
            Prediction = .data$prediction, .data$Selected
        ) %>%
        dplyr::bind_rows(
            dplyr::tibble(
                rlsample = object@metadata$sampleName,
                Mode = object@metadata$mode,
                Label = object@metadata$label,
                Prediction = prediction$prediction,
                Selected = object@metadata$sampleName
            )
        ) %>%
        dplyr::mutate(
            Mode = dplyr::case_when(
                .data$Mode %in% auxdata$misc_modes ~ "misc",
                TRUE ~ .data$Mode
            )
        ) %>%
        dplyr::distinct(.data$rlsample, .keep_all = TRUE)
    annoCorr <- as.data.frame(annoCorr)
    rownames(annoCorr) <- annoCorr$rlsample
    annoCorr <- annoCorr[, -which(colnames(annoCorr) == "rlsample")]

    # Filter for available / desired samples
    toSelect <- colnames(corrRes)
    if (!is.null(selected)) {
        toSelect <- intersect(selected, toSelect)
    }
    corrRes <- corrRes[toSelect, toSelect]
    annoCorr <- annoCorr[toSelect, ]

    # Filter to remove any un-predicted samples
    keep <- which(!is.na(annoCorr$Prediction))
    corrRes <- corrRes[keep, keep]
    annoCorr <- annoCorr[keep, ]

    # Pallete
    paletteLength <- 250
    myColor <- grDevices::colorRampPalette(
        rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))
    )(paletteLength)
    myBreaks <- c(
        seq(-1, 0, length.out = ceiling(paletteLength / 2) + 1),
        seq(
            1 / paletteLength, 1,
            length.out = floor(paletteLength / 2)
        )
    )

    # Wrangle colors
    mode_cols <- auxdata$mode_cols$col
    names(mode_cols) <- auxdata$mode_cols$mode
    cond_cols <- auxdata$label_cols$col
    names(cond_cols) <- auxdata$label_cols$label
    verd_cols <- auxdata$prediction_cols$col
    names(verd_cols) <- auxdata$prediction_cols$prediction
    group_cols <- stats::setNames(c(
        auxdata$heat_cols$col[auxdata$heat_cols$selected == "user_selected"],
        auxdata$heat_cols$col[auxdata$heat_cols$selected == "RLBase"]
    ), nm = c(object@metadata$sampleName, "RLBase"))
    cat_cols <- list(
        "Mode" = mode_cols,
        "Label" = c(cond_cols, "grey"),
        "Prediction" = verd_cols,
        "Selected" = group_cols
    )
    cat_cols$Mode <- cat_cols$Mode[names(cat_cols$Mode) %in% annoCorr$Mode]
    continuous_pal <- circlize::colorRamp2(
        breaks = myBreaks[-1],
        colors = myColor
    )

    # Return data if requested
    if (returnData) {
        return(
            list(
                "corrRes" = corrRes,
                "continuous_pal" = continuous_pal,
                "cat_cols" = cat_cols,
                "annoCorr" = annoCorr,
                "pheatmap_breaks" = myBreaks,
                "pheatmap_color" = myColor
            )
        )
    }

    # Build heatmap
    if (complex) {
        hm <- ComplexHeatmap::Heatmap(
            corrRes,
            col = continuous_pal,
            row_dend_reorder = FALSE,
            column_dend_reorder = FALSE,
            show_column_names = FALSE,
            show_row_names = FALSE,
            top_annotation = ComplexHeatmap::HeatmapAnnotation(
                df = annoCorr, col = cat_cols
            ),
            name = "Corr (R)"
        )
    } else {
        hm <- pheatmap::pheatmap(
            corrRes,
            color = myColor, breaks = myBreaks,
            annotation_col = annoCorr[, c(3, 2, 1)],
            annotation_colors = cat_cols,
            show_colnames = FALSE,
            show_rownames = FALSE,
            silent = TRUE
        )
    }
    return(hm)
}


#' Plot Enrichment Test Results
#'
#' Creates a list of plots, one for each annotation database
#' (see [RLHub::annotations]).
#' These plots show the feature enrichment for the user-supplied sample in 
#' comparison to the samples in 
#' [RLBase](https://gccri.bishop-lab.uthscsa.edu/rlbase/). This
#' will only work if you did not use custom annotations with [featureEnrich].
#' 
#' @param object An RLRanges object with [featureEnrich] already run.
#' @param pred_POS_only If TRUE, only "POS" predicted samples included (see 
#' also [predictCondition]). Default: TRUE. 
#' @param label_POS_only If TRUE, only "POS" labeled samples included (samples
#' which are expected to robustly map R-loops,
#' e.g., "D210N" condition R-ChIP data). Default: FALSE.
#' @param splitby Metadata by which to split plots. Can be "none", "prediction",
#'  or "label".
#' @param limits Specify limits on data range. This is used for controlling 
#' the infinite estimation of odds ratio resulting from fisher's exact test. 
#' To remove limits, set c(-Inf, Inf). Default: c(-10, 15).
#' @param returnData If TRUE, plot data is returned instead of plot objects.
#'  Default: FALSE
#' @param ... For internal use.
#' @return A named list of [ggplot2::ggplot] objects. Names correspond to 
#' the annotations provided. See also [featureEnrich].
#' @examples
#'
#' # Example dataset with featureEnrich() already run.
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Make plots, split by prediction
#' plotEnrichment(rlr, pred_POS_only = FALSE, splitby = "prediction")
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
#' @export
plotEnrichment <- function(object,
    pred_POS_only = TRUE,
    label_POS_only = FALSE,
    splitby = c("none", "prediction", "label"),
    limits = c(-10, 15),
    returnData = FALSE,
    ...) {

    # TODO: Should there be an option to control this for the user?
    yval <- "stat_fisher_rl"

    # Get dots -- get values used by RLBase app
    # Provides RLBase with a speed boost by pre-supplying
    # data from memory when plotting. Not intended for regular use.
    dots <- list(...)
    rlbaseRes <- NULL
    rlsamples <- NULL
    if (length(dots) > 0) {
        rlbaseRes <- dots$rlbaseRes
        rlsamples <- dots$rlsamples
    }

    # Verify splitby
    stopifnot(splitby[1] %in% c("none", "prediction", "label"))

    # Get sample-level feature enrichment and RLBase samples
    if (is.null(rlbaseRes)) rlbaseRes <- RLHub::feat_enrich_samples(quiet = TRUE)
    if (is.null(rlsamples)) rlsamples <- RLHub::rlbase_samples(quiet = TRUE)

    # Get enrichment results from object
    sampleRes <- rlresult(object, resultName = "featureEnrichment")

    # Add metadata
    usamp <- object@metadata$sampleName
    sampleRes$experiment <- usamp
    sampleRes$label <- object@metadata$label
    sampleRes$mode <- object@metadata$mode
    predres <- rlresult(object, resultName = "predictRes")
    sampleRes$prediction <- predres$prediction

    # Filter RLBase
    if (pred_POS_only) {
        rlsamples <- dplyr::filter(rlsamples, .data$prediction == "POS")
    }
    if (label_POS_only) {
        rlsamples <- dplyr::filter(rlsamples, .data$label == "POS")
    }

    # Wrangle the RLBase data
    rlbaseRes <- rlbaseRes[rlbaseRes$db %in% sampleRes$db, ]
    rlbaseResFull <- dplyr::bind_rows(
        rlbaseRes
    ) %>%
        dplyr::inner_join(
            rlsamples,
            by = c("experiment" = "rlsample")
        )

    # Combine user sample and RLBase
    input_data <- sampleRes %>%
        dplyr::bind_rows(
            rlbaseResFull
        ) %>%
        dplyr::select(
            .data$db, .data$type, .data$experiment,
            .data$prediction, .data$label, .data$mode,
            dplyr::contains("stat_fisher_rl")
        ) %>%
        dplyr::mutate(
            # TODO: Find a better way to deal with Inf results.
            stat_fisher_rl = log2(.data$stat_fisher_rl),
            stat_fisher_rl = ifelse(
                .data$stat_fisher_rl > limits[2],
                limits[2],
                .data$stat_fisher_rl
            ),
            stat_fisher_rl = ifelse(
                .data$stat_fisher_rl < limits[1],
                limits[1],
                .data$stat_fisher_rl
            )
        ) %>%
        dplyr::filter(.data$experiment %in% c(rlsamples$rlsample, usamp)) %>%
        dplyr::filter(
            !is.na(.data$stat_fisher_rl)
        )

    # Wrap strings
    input_data$type <- gsub(input_data$type, pattern = "_", replacement = " ")

    # Make selected explicit
    input_data$selected <- ifelse(input_data$experiment == usamp, usamp, "")
    input_data <- dplyr::distinct(
        input_data, .data$db, .data$type, .data$experiment, .keep_all = TRUE
    )

    # Build plots
    datalst <- input_data %>%
        dplyr::group_by(.data$db) %>%
        {
            stats::setNames(dplyr::group_split(.), dplyr::group_keys(.)[[1]])
        }

    # Return data if requested
    if (returnData) {
        return(datalst)
    }

    # Get plots
    plts <- lapply(
        datalst,
        feature_ggplot,
        usamp = usamp,
        limits = limits,
        splitby = splitby
    )

    # Remove empty
    plts <- plts[!vapply(plts, is.null, logical(1))]

    # Return
    return(plts)
}


#' Feature ggplot
#'
#' The core plotting component of `plotEnrichment`
#'
#' @param x A `tbl` containing data for plotting.
#' @param usamp The name of the user-supplied sample
#' @param limits Specify limits on data. Useful for controlling infinite
#' estimation of odds ratio
#' resulting from fisher's exact test. To remove limits, set c(-Inf, Inf).
#'  Default: c(-10, 15).
#' @param splitby Metadata by which to split plots. Can be "none", "prediction",
#'  or "label".
#'
#' @returns A ggplot2 object.
feature_ggplot <- function(x, usamp, limits, splitby) {
    yval <- "stat_fisher_rl"
    db_now <- x$db[1]
    col <- auxdata$db_cols$col[auxdata$db_cols$db == db_now]

    if (!usamp %in% x$experiment) {
        warning(
            "User-supplied sample test value is NA for ",
            db_now, ". Cannot plot."
        )
    }

    # Define limits
    lmts <- c(
        max(min(x[, yval]), limits[1]) - 1,
        min(max(x[, yval]), limits[2]) + 1
    )


    # Make base plot
    pltbase <- ggplot2::ggplot(
        x,
        ggplot2::aes_string(
            x = "type",
            y = yval,
            label = "experiment",
            fill = if (splitby[1] != "none") splitby[1],
            color = if (splitby[1] != "none") splitby[1]
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed")


    if (splitby[1] == "none") {
        plt <- pltbase +
            ggplot2::geom_violin(
                color = "black",
                width = .8,
                fill = col,
                trim = FALSE,
                position = ggplot2::position_dodge(.8),
                alpha = .25,
                na.rm = TRUE
            ) +
            ggplot2::geom_boxplot(
                width = .2,
                color = "black",
                fill = col,
                position = ggplot2::position_dodge(.8),
                alpha = .5
            ) +
            ggplot2::geom_jitter(
                alpha = ifelse(x$experiment == usamp, 1, 0),
                size = 3,
                fill = col,
                width = 0,
                color = "black",
                shape = 23,
                stroke = 1.5
            )
    } else {
        # Get split cols
        cols <- auxdata[[paste0(tolower(splitby[1]), "_cols")]]
        colvec <- cols %>% dplyr::pull(.data$col)
        names(colvec) <- as.data.frame(cols)[, splitby[1]]

        # Add box and jitter
        plt <- pltbase +
            ggplot2::geom_violin(
                color = "black",
                width = .8,
                trim = FALSE,
                position = ggplot2::position_dodge(.8),
                alpha = .5,
                na.rm = TRUE
            ) +
            ggplot2::geom_boxplot(
                width = .2,
                color = "black",
                position = ggplot2::position_dodge(.8),
                alpha = 1
            ) +
            ggplot2::geom_jitter(
                alpha = ifelse(x$experiment == usamp, 1, 0),
                size = 3,
                color = "black",
                position = ggplot2::position_jitterdodge(
                    dodge.width = .8,
                    jitter.width = 0,
                    jitter.height = 0
                ),
                shape = 23,
                stroke = 1.5
            ) +
            ggplot2::scale_fill_manual(
                values = colvec, drop = TRUE
            )
    }

    # Final additions to plot and themeing
    # TODO: Want a more systematic ylab method
    plt +
        ggplot2::ylab(
            "Fisher's test odds ratio (log2)"
        ) +
        ggplot2::xlab(NULL) +
        ggplot2::labs(
            title = gsub(
                db_now,
                pattern = "_",
                replacement = " "
            ),
            subtitle = usamp,
            caption = paste0(
                "y-axis data min, max: ",
                paste0(limits, collapse = ", "),
                ". \u25C7 - User sample"
            )
        ) +
        ggprism::theme_prism(base_size = 14) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 45,
                vjust = 1,
                hjust = 1
            )
        ) +
        ggplot2::theme(
            legend.title = ggplot2::element_text(size = 18),
            legend.text = ggplot2::element_text(size = 14),
            plot.caption = ggplot2::element_text(size = 11)
        ) +
        ggplot2::guides(
            fill = ggplot2::guide_legend(
                override.aes = list(size = 0, stroke = 0)
            )
        ) +
        ggplot2::scale_y_continuous(limits = lmts)
}


#' Plot RL-Region overlap with RLRanges
#'
#' Convenience function for plotting the overlap between RLRanges and R-loop
#' regions (RL regions) as calculated by [rlRegionTest].
#'
#' @param object An RLRanges object with [rlRegionTest] already run.
#' @param returnData If TRUE, plot data is returned instead of plotting.
#'  Default: FALSE
#' @param rlregions_table The table of RLRegions to overlap sample ranges with.
#' Obtained from RLHUb using [RLHub::rlregions_meta]. Loaded from RLHub if
#' not supplied. Default: NULL.
#' @param ... Additional arguments passed to [VennDiagram::venn.diagram]
#' @return A [ggplot2::ggplot] object containing the venn diagram. Built 
#' using [ggplotify::as.ggplot].
#' @examples
#'
#' # Example dataset with rlRegionTest() already run.
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Plot RL-Region overlap
#' plotRLRegionOverlap(rlr)
#'
#' # Return data only
#' plotRLRegionOverlap(rlr, returnData = TRUE)
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
#' @export
plotRLRegionOverlap <- function(object, returnData = FALSE,
    rlregions_table = NULL, ...) {

    # Get RLRegions
    if (is.null(rlregions_table)
    ) {
        rlregions_table <- RLHub::rlregions_meta(quiet = TRUE)
    }

    # Retrieve results
    olres <- rlresult(object, "rlRegionRes")

    # Get abstracted sites
    pkolnms <- olres$Overlap$name__peaks
    rlolnms <- olres$Overlap$name__rlregion

    # Get overlapping ranges
    pkolgr <- object[names(object) %in% pkolnms, ]

    # Get the ranges for the rlregions
    rlReg <- tableToRegions(rlregions_table)
    rlolgr <- rlReg %>%
        dplyr::filter(.data$name %in% {{ rlolnms }}) %>%
        GenomicRanges::makeGRangesFromDataFrame()

    # Get the number of shared overlapping sites
    shared <- c(pkolgr, rlolgr) %>%
        GenomicRanges::reduce() %>%
        length()

    # Get number unique
    pkonly <- length(object[!names(object) %in% pkolnms, ])
    rlonly <- rlReg %>%
        dplyr::filter(!.data$name %in% {{ rlolnms }}) %>%
        nrow()

    # Build the tbl and plot
    sites <- seq(pkonly + rlonly + shared)

    # Get the plot data
    pltdata <- list(
        sites[c(rep(FALSE, pkonly), rep(TRUE, rlonly), rep(TRUE, shared))],
        sites[c(rep(TRUE, pkonly), rep(FALSE, rlonly), rep(TRUE, shared))]
    ) %>%
        stats::setNames(nm = c("RL-Regions", object@metadata$sampleName))
    subtitle <- paste0(
        "Fisher's Test: pval = ",
        signif(olres$Test_results$p.value),
        "; odds ratio = ",
        signif(olres$Test_results$estimate)
    )
    retData <- list("pltdata" = pltdata, "subtitle" = subtitle)
    if (returnData) {
        return(retData)
    }

    # Suppress output log from venn.diagram if futile logger available
    if (requireNamespace("futile.logger", quietly = TRUE)) {
        futile.logger::flog.threshold(
            futile.logger::ERROR, name = "VennDiagramLogger"
        )
    }
    
    # Make the plot
    gt <- pltdata %>%
        VennDiagram::venn.diagram(
            filename = NULL,
            main = "RL-Region Overlap",
            sub = subtitle,
            ...
        ) %>%
        grid::grid.draw() %>%
        ggplotify::grid2grob() %>%
        ggplotify::as.ggplot()

    # Return
    return(gt)
}
