#' Plot Perm Test results
#'
#' @param object An RLRanges with \code{analyzeRLFS()} already run.
#' @param ... Additional parameters passed to \code{ggplot}.
#' @return A ggplot object.
#' @export
plotRLFSRes <- function(object,
                        ...) {
  
  # Obtain RLFS-Res
  rlfsRes <- rlresult(object, resultName = "rlfsRes")
  
  # Obtain the plot data
  pltdat <- tibble::tibble(
    "zscore" = rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifted.z.scores,
    "shift" = rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifts
  )
  pval <- rlfsRes$perTestResults$`regioneR::numOverlaps`$pval

  # Make plot
  ggplot2::ggplot(
    data = pltdat,
    ggplot2::aes_string(
      y = "zscore",
      x = "shift"
    )
  ) +
    ggplot2::geom_vline(
      color = "firebrick",
      xintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::geom_line(size = 1) +
    ggplot2::ggtitle("Z-score around RLFS") +
    ggplot2::labs(caption = paste0("p < ", round(pval, digits = 5))) +
    ggplot2::ylab("Peak Enrichment (Z-Score)") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::xlab("Distance to RLFS (bp)") +
    ggprism::theme_prism(base_size = 14) +
    ggplot2::labs(subtitle = object@metadata$sampleName)
}


#' Plot Correlation Results
#'
#' @param object An RLRanges with \code{analyzeRLFS()} already run.
#' @param ... For internal use.
#' @return A Heatmap object.
#' @export
corrHeatmap <- function(object, ...) {

  # TODO: NEED RLHub for this
  rlsamples <- file.path(rlbase, "RLHub", "rlsamples.rda")
  tmp <- tempfile()
  download.file(rlsamples, destfile = tmp, quiet = TRUE)
  load(tmp)
  
  # Get the correlation matrix
  corrRes <- rlresult(object, resultName = "correlationMat")
  
  # Get the mode and prediction and condType
  prediction <- rlresult(object, resultName = "predictRes")
  
  # Get dots -- get values used by shiny
  dots <- list(...)
  selected <- NA
  if (length(dots) > 0) {
    selected <- dots$selected
  }
  
  # Wrangle the annotation data
  annoCorr <- rlsamples %>%
    dplyr::mutate(group = "RLBase") %>%
    dplyr::select(
      .data$rlsample, .data$mode, 
      .data$condType,
      .data$verdict, .data$group
    ) %>%
    dplyr::bind_rows(
      tibble::tibble(
        rlsample = object@metadata$sampleName,
        mode = object@metadata$mode,
        condType = object@metadata$condType,
        verdict = prediction$Verdict,
        group = "user_selected"
      )
    ) %>%
    dplyr::distinct(.data$rlsample, .keep_all = TRUE) %>%
    tibble::column_to_rownames(var = "rlsample")

  # Filter for available / desired samples
  toSelect <- colnames(corrRes)
  if (!is.na(selected)) {
    toSelect <- intersect(selected, toSelect)
  }
  corrNow <- corrRes[toSelect, toSelect]
  annoCorr <- annoCorr[toSelect, ]

  # Pallete
  paletteLength <- 100
  myColor <- colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))
  )(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(
    seq(min(corrNow), 0, length.out = ceiling(paletteLength / 2) + 1),
    seq(max(corrNow) / paletteLength, max(corrNow), 
        length.out = floor(paletteLength / 2))
  )
  
  # Wrangle colors
  mode_cols <- aux$mode_cols$POS
  names(mode_cols) <- aux$mode_cols$mode
  cond_cols <- aux$condtype_cols$col
  names(cond_cols) <- aux$condtype_cols$condType
  verd_cols <- aux$verdict_cols$col
  names(verd_cols) <- aux$verdict_cols$verdict
  group_cols <- c(
    "user_selected" = "firebrick",
    "RLBase" = "grey"
  )
  collst <- list(
    "mode" = mode_cols,
    "condType" = c(cond_cols, "grey"),
    "verdict" = verd_cols,
    "group" = group_cols
  )
  
  collst$mode <- collst$mode[which(names(collst$mode) %in% annoCorr$mode)]
  
  # Build heatmap
  chm <- ComplexHeatmap::Heatmap(
    corrRes, col = circlize::colorRamp2(breaks = myBreaks[-1], colors = myColor), 
    show_column_names = FALSE,
    show_row_names = FALSE,
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      df = annoCorr, col = collst
    ),
    name = "Corr (R)"
  )
  
  return(chm)
  
}


#' Plot Enrichment Test Results
#'
#' @param object The tibble obejct obtained from running \code{featureEnrich}.
#' @param rlbaseRes A tibble containing results from RLBase samples.
#' @param rlsamples A tibble containing the metadata of RLBase samples.
#' @param yval Column to use for obtaining Y-axis values. Default: "stat_fisher_rl".
#' If NULL, no comparison with RLBase will be made.
#' @return A named list of ggplot objects.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
plotEnrichment <- function(object,
                           rlbaseRes = NULL,
                           rlsamples = NULL,
                           splitby = c("none", "verdict", "condtype"),
                           facetby = c("none", "mode"),
                           yval = "stat_fisher_rl") {
  
  # TODO: NEEDS to be replaced with RLHub
  rlbase_enrich <- file.path(rlbase, "RLHub", "feature_enrichment_per_sample.rda")
  tmp <- tempfile()
  download.file(rlbase_enrich, destfile = tmp, quiet = TRUE)
  load(tmp)
  rlbaseRes <- feature_enrichment_per_sample
  
  # TODO: NEEDS to be replaced with RLHub
  rlbase_samples <- file.path(rlbase, "RLHub", "rlsamples.rda")
  tmp <- tempfile()
  download.file(rlbase_samples, destfile = tmp, quiet = TRUE)
  load(tmp)
  
  # Get enrichment results from object
  sampleRes <- rlresult(object, resultName = "featureEnrichment")
  
  # Get min, max for plotting
  limit <- ifelse(yval == "stat_fisher_rl", 10, Inf)

  # Wranlg the input data
  input_data <- sampleRes %>%
    dplyr::mutate(experiment = "User-supplied") %>%
    dplyr::bind_rows(
      dplyr::filter(
        rlbaseRes, .data$type %in% sampleRes$type
      )
    ) %>%
    dplyr::select(
      .data$db, .data$type, .data$experiment,
      dplyr::contains("stat_fisher_rl")
    ) %>%
    dplyr::mutate(
      stat_fisher_rl = log2(.data$stat_fisher_rl)
    ) %>%
    dplyr::filter(
      !is.na(.data$stat_fisher_rl),
      is.finite(.data$stat_fisher_rl)
    ) %>%
    dplyr::left_join(
      rlsamples,
      by = c("experiment" = "rlsample")
    )
  if (splitby[1] == "none") {
    input_data <- dplyr::filter(
      input_data,
      .data$verdict == "Case" | .data$experiment == "User-supplied"
    )
  }

  # Wrap strings
  input_data$type <- gsub(input_data$type, pattern = "_", replacement = " ")
  input_data$type[nchar(input_data$type) > 30] <-
    stringr::str_wrap(input_data$type[nchar(input_data$type) > 30], 30)

  # Build plots
  datalst <- input_data %>%
    dplyr::group_by(.data$db) %>%
    {setNames(dplyr::group_split(.), dplyr::group_keys(.)[[1]])}
  plts <- lapply(
    datalst,
    function(x) {
      db_now <- x$db[1]
      col <- aux$db_cols$col[aux$db_cols$db == db_now]
      
      if (!"User-supplied" %in% x$experiment) {
        return(NULL)
      }
      
      ggplot2::ggplot(
        x,
        ggplot2::aes_string(
          x = "type",
          y = yval,
          label = "experiment"
        )
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::geom_boxplot(
          alpha = .25,
          col = col,
          fill = col,
          outlier.shape = NA,
          outlier.size = 0
        ) +
        ggplot2::geom_jitter(
          color = col,
          fill = "#dbdbdb",
          alpha = 0.3,
          size = 0.3,
          width = .3
        ) +
        ggplot2::geom_point(
          data = x[x$experiment == "User-supplied", ],
          mapping = ggplot2::aes_string(
            col = "experiment"
          ),
          size = 2.5
        ) +
        ggplot2::scale_color_manual(
          values = c("User-supplied" = "black")
        ) +
        ggpubr::rremove("legend") +
        ggplot2::xlab(NULL) +
        ggplot2::ylab(
          stringr::str_to_title(
            gsub(yval,
                 pattern = "_",
                 replacement = " "
            )
          )
        ) +
        ggplot2::labs(
          title = gsub(db_now,
                       pattern = "_",
                       replacement = " "
          )
        ) +
        ggprism::theme_prism(base_size = 13) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = 45,
            vjust = 1,
            hjust = 1
          ),
          legend.title = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = 10),
          legend.position = "bottom"
        ) +
        ggplot2::scale_y_continuous(limits = c(-limit, limit))
    }
  )

  # Remove empty
  plts <- plts[!sapply(plts, is.null)]
  
  # Return
  return(plts)
}

