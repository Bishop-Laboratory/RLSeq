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
    pltxlab <- ifelse(fft, "Frequency Domain", "Distance to RLFS (bp)")
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
    pltbase <- pltbase +
        ggplot2::geom_line(size = 1) +
        ggplot2::ggtitle(plttitle) +
        ggplot2::labs(caption = paste0("p < ", round(pval, digits = 5))) +
        ggplot2::ylab(pltylab) +
        ggplot2::xlab(pltxlab) +
        ggprism::theme_prism(base_size = 14) +
        ggplot2::labs(subtitle = plotName)
    if (fft) {
        # Remove xaxis ticks / numbers if plotting fourier frequency
        pltbase <- pltbase +
            ggplot2::theme(
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank()
            )
    }
    return(pltbase)
}


#' Plot noise analysis results as a fingerprint plot
#'
#' Plots the results of the noise analysis in [noiseAnalyze]. Creates a
#' fingerprint plot like those developed by
#' [Diaz et al, 2012](https://pubmed.ncbi.nlm.nih.gov/22499706/) and those
#' provided by
#' [deepTools](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html).
#'
#' The term "Fingerprint plot" comes from deepTools.
#'
#' @param object An RLRanges object with [noiseAnalyze] already run.
#' @return A ggplot object. See also [ggplot2::ggplot].
#'
#' 2. `noiseComparisonPlot`
#'   - A plot showing the noise analysis
#'   results from the user-supplied sample compared to similar samples from
#'   RLBase.
#'
#' @examples
#'
#' # Example RLRanges dataset with analyzeRLFS() already run.
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Plot RLFS res
#' plotFingerprint(rlr)
#'
#' @export
plotFingerprint <- function(object) {

    # Get the noise analysis results
    rlnoise <- rlresult(object, "noiseAnalysis")

    # Get sample name
    sname <- object@metadata$sampleName

    # Create a fingerprint plot
    fingerprint <- rlnoise %>%
        dplyr::mutate(
            rank = .data$rank / max(.data$rank)
        ) %>%
        ggplot2::ggplot(
            ggplot2::aes_string(
                x = "rank",
                y = "value"
            )
        ) +
        ggplot2::geom_abline(
            slope = 1, intercept = 0,
            color = "grey",
            linetype = "dashed",
            alpha = .5
        ) +
        ggplot2::geom_line() +
        ggplot2::labs(
            title = paste0("Fingerprint plot of ", sname)
        ) +
        ggplot2::xlab("Proportion of bins") +
        ggplot2::ylab("Proportion of signal") +
        ggplot2::theme_linedraw(
            base_size = 14
        ) +
        ggplot2::theme(
            plot.caption = ggplot2::element_text(size = 8)
        )

    return(fingerprint)
}


#' Creates a metaplot for comparing noise analysis results with RLBase
#'
#' Plots the average standardized signal from [noiseAnalyze] alongside the
#' samples in RLBase. For this plot, lower average signal indicates better
#' signal to noise ratio. **Note**: This plot may be misleading if you supplied
#' custom windows when running [noiseAnalyze].
#'
#' @param object An RLRanges object with [noiseAnalyze] already run.
#' @param mode A `character` containing the R-loop data mode to compare
#' against. See details for more information.
#' @param simple A `logical` which specifies whether the plot should only show
#' samples where the prediction and label are the same. Default: TRUE.
#' @param returnData If TRUE, plot data is returned instead of plotting.
#'  Default: FALSE
#' @return A [ggplot2::ggplot] object or a `tbl` if `returnData` is `TRUE`.
#' @details
#'
#' ## Mode
#'
#' The `mode` parameter specifies the R-loop modality to compare the
#' user-supplied sample against in the plot. The default, `"auto"`
#' specifies that the `mode` from the supplied `RLRanges` object will
#' be used. Only one mode can be specified. For a list of applicable modes,
#' see `auxdata$available_modes`.
#'
#' ## Plot
#'
#' The plot is a violin / jitter plot showing the distribution of average values
#' from the [noiseAnalyze] output across RLBase samples of the selected `mode`.
#' The user-supplied sample is annotated on the plot.
#'
#' @examples
#'
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Plot RL-Region overlap
#' noiseComparisonPlot(rlr)
#'
#' # Return data only
#' noiseComparisonPlot(rlr, returnData = TRUE)
#' @export
noiseComparisonPlot <- function(object, mode = "auto", simple = TRUE, returnData = FALSE) {

    # Retrieve results
    noiseres <- rlresult(object, "noiseAnalysis")

    # Get the noise index
    noiseindex <- mean(noiseres$value)

    # Wrangle along with sample info
    userdata <- dplyr::tibble(
        sample = object@metadata$sampleName,
        noise_index = noiseindex,
        label = object@metadata$label,
        prediction = object@metadata$results@predictRes$prediction,
        group = "User-supplied"
    )

    # Get genome
    genome <- GenomeInfoDb::genome(object)[1]

    # Get rlbase samps
    rlsamples <- RLHub::rlbase_samples()

    # Get the mode
    if (mode == "auto") mode <- object@metadata$mode

    # Get the RLBase noise analysis results and wrangle them
    toplt <- rlbaseNoiseAnalyze %>%
        dplyr::inner_join(rlsamples, by = "rlsample") %>%
        dplyr::filter(
            .data$mode == {{ mode }},
            .data$genome == {{ genome }}
        ) %>%
        dplyr::mutate(group = {{ mode }}) %>%
        dplyr::select(
            sample = .data$rlsample,
            noise_index = .data$value,
            .data$label, .data$prediction, .data$group
        ) %>%
        dplyr::bind_rows(userdata)

    if (returnData) {
        return(toplt)
    }

    # Plot
    if (!simple) {
        cond <- "cond"
        colvec <- auxdata$prediction_label_cols
        xlab <- "Sample QC prediction/label (prediction_label)"
        title <- "Noise index across prediction-label combinations"
        topltfinal <- toplt
    } else {
        cond <- "prediction"
        colvec <- stats::setNames(
            auxdata$prediction_cols$col,
            nm = auxdata$prediction_cols$prediction
        )
        xlab <- "Sample QC prediction"
        topltfinal <- toplt %>%
            dplyr::filter(
                .data$prediction == .data$label | .data$group == "User-supplied"
            )
        title <- "Noise comparison plot"
    }
    plt <- topltfinal %>%
        dplyr::mutate(
            cond = paste0(.data$prediction, "_", .data$label),
            cond = factor(.data$cond, levels = c("NEG_NEG", "NEG_POS", "POS_NEG", "POS_POS"))
        ) %>%
        ggplot2::ggplot(ggplot2::aes_string(x = cond, fill = cond, y = "noise_index")) +
        ggplot2::geom_boxplot(
            width = ifelse(simple, .8, 1),
            position = ggplot2::position_dodge(width = .9),
            outlier.shape = NA
        ) +
        ggplot2::geom_point(
            position = ggplot2::position_jitterdodge(
                dodge.width = .9,
                jitter.width = ifelse(simple, 0.25, 1),
                seed = 42
            )
        ) +
        ggplot2::geom_point(
            alpha = ifelse(topltfinal$group == "User-supplied", 1, 0),
            size = 4,
            color = "black",
            position = ggplot2::position_jitterdodge(
                dodge.width = .9,
                jitter.width = 0,
                jitter.height = 0
            ),
            shape = 23,
            stroke = 1.5
        ) +
        ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, .2))) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::scale_fill_manual(
            values = colvec
        ) +
        ggplot2::ylab("Noise index (log scale)") +
        ggplot2::xlab(xlab) +
        ggplot2::labs(
            title = title,
            subtitle = paste0("Data modality: ", mode, "-seq"),
            caption = paste0(
                "\u25C7 - User-supplied sample"
            )
        ) +
        ggplot2::theme(legend.position = "none")

    # Return plt
    return(plt)
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
            .data$rlsample,
            Mode = .data$mode,
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
        input_data, .data$db, .data$type, .data$experiment,
        .keep_all = TRUE
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


#' Plot Transcript Feature Overlap
#'
#' Plots the results of [txFeatureOverlap] alongside the average from
#' public R-loop datasets. This allows comparison of user-supplied
#' samples with data that is expected to be simialr.
#'
#' @param object An RLRanges object with [txFeatureOverlap] already run.
#' @param mode A `character` containing the R-loop data mode to compare
#' against. See details for more information.
#' @param returnData If TRUE, plot data is returned instead of plotting.
#'  Default: FALSE
#' @return A [ggplot2::ggplot] object or a `tbl` if `returnData` is `TRUE`.
#' @details
#'
#' ## Mode
#'
#' The `mode` parameter specifies the R-loop modality to compare the
#' user-supplied sample against in the plot. The default, `"auto"`
#' specifies that the `mode` from the supplied `RLRanges` object will
#' be used. Only one mode can be specified. For a list of applicable modes,
#' see `auxdata$available_modes`.
#'
#' ## Plot
#'
#' The plot is a stacked bar chart showing the proportion of peaks overlapping
#' transcript features for the supplied `RLRanges` object. Additionally, the
#' average of the [txFeatureOverlap] analysis for all samples within the
#' specified modes are also shown as a background comparison.
#'
#' This style of analysis enables a user to see the transcript features
#' overlapping their peaks and compare those results to the average within
#' relevant public datasets.
#'
#' @examples
#'
#' rlr <- readRDS(system.file("extdata", "rlrsmall.rds", package = "RLSeq"))
#'
#' # Plot RL-Region overlap
#' plotTxFeatureOverlap(rlr)
#'
#' # Return data only
#' plotRLRegionOverlap(rlr, returnData = TRUE)
#' @export
plotTxFeatureOverlap <- function(object, mode = "auto", returnData = FALSE) {

    # Retrieve results
    txfeatres <- rlresult(object, "txFeatureOverlap")

    # Get genome
    genome <- GenomeInfoDb::genome(object)[1]

    # Get rlbase samps
    rlsamples <- RLHub::rlbase_samples()

    # Get the mode
    if (mode == "auto") mode <- object@metadata$mode

    # Get the rlbase tx feat overlap
    rlspos <- rlsampleTxOl %>%
        dplyr::inner_join(rlsamples, by = "rlsample") %>%
        dplyr::filter(
            .data$label == "POS",
            .data$mode == {{ mode }},
            .data$numPeaks > 5000,
            .data$prediction == "POS",
            .data$genome == {{ genome }}
        ) %>%
        dplyr::group_by(.data$feature) %>%
        dplyr::summarise(
            pct = mean(.data$pct),
        ) %>%
        dplyr::mutate(group = paste0({{ mode }}, "-avg"))

    # Get stats for the user-supplied sample
    olstat <- txfeatres %>%
        dplyr::group_by(.data$feature) %>%
        dplyr::tally() %>%
        dplyr::mutate(
            pct = .data$n / sum(.data$n),
            group = object@metadata$sampleName
        ) %>%
        dplyr::select(-.data$n)
    toplt <- olstat %>%
        dplyr::bind_rows(rlspos) %>%
        dplyr::mutate(
            feature = factor(.data$feature, levels = rev(c(
                "TSS", "fiveUTR", "Exon", "Intron",
                "threeUTR", "TTS", "Intergenic"
            )))
        )
    if (returnData) {
        return(toplt)
    }

    # Plot
    plt <- toplt %>%
        ggplot2::ggplot(ggplot2::aes_string(x = "group", y = "pct", fill = "feature")) +
        ggplot2::geom_col(color = "black") +
        ggplot2::coord_flip() +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
            axis.title.y = ggplot2::element_blank()
        ) +
        ggplot2::ylab("Proportion of peaks") +
        ggplot2::ggtitle(
            "Transcript Feature Overlap",
            subtitle = paste0(
                object@metadata$sampleName,
                " (user-supplied data) vs RLBase ",
                paste0(mode, "-avg")
            )
        ) +
        ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(reverse = TRUE))

    # Return plt
    return(plt)
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
            futile.logger::ERROR,
            name = "VennDiagramLogger"
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
