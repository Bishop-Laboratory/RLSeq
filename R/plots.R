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
#' @param returnData If TRUE, plot data is returned instead of plotting. Default: FALSE
#' @param ... For internal use.
#' @return A Heatmap object.
#' @export
corrHeatmap <- function(object,
                        returnData=FALSE,
                        ...) {

  # TODO: NEED RLHub for this
  rlbase <- "https://rlbase-data.s3.amazonaws.com"
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
      # .data$condType,
      .data$verdict, .data$group
    ) %>%
    dplyr::bind_rows(
      tibble::tibble(
        rlsample = object@metadata$sampleName,
        mode = object@metadata$mode,
        # condType = object@metadata$condType,
        verdict = prediction$Verdict,
        group = object@metadata$sampleName
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
  mode_cols <- aux$mode_cols$col
  names(mode_cols) <- aux$mode_cols$mode
  cond_cols <- aux$condtype_cols$col
  names(cond_cols) <- aux$condtype_cols$condType
  verd_cols <- aux$verdict_cols$col
  names(verd_cols) <- aux$verdict_cols$verdict
  group_cols <- setNames(c(
    aux$heat_cols$col[aux$heat_cols$selected == "user_selected"],
    aux$heat_cols$col[aux$heat_cols$selected == "RLBase"]
  ), nm = c(object@metadata$sampleName, "RLBase"))
  cat_cols <- list(
    "mode" = mode_cols,
    # "condType" = c(cond_cols, "grey"),
    "verdict" = verd_cols,
    "group" = group_cols
  )
  cat_cols$mode <- cat_cols$mode[which(names(cat_cols$mode) %in% annoCorr$mode)]
  continuous_pal <- circlize::colorRamp2(breaks = myBreaks[-1], colors = myColor)
  
  # Return data if requested
  if (returnData) {
    return(
      list(
        "corrRes" = corrRes,
        "continuous_pal" = continuous_pal,
        "cat_cols" = cat_cols,
        "annoCorr" = annoCorr
      )
    )
  }
  
  # Build heatmap
  chm <- ComplexHeatmap::Heatmap(
    corrRes, col = continuous_pal, 
    show_column_names = FALSE,
    show_row_names = FALSE,
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      df = annoCorr, col = cat_cols
    ),
    name = "Corr (R)"
  )
  
  return(chm)
  
}


#' Plot Enrichment Test Results
#'
#' @param object The tibble obejct obtained from running \code{featureEnrich}.
#' @param modes Which modes to include in plot? If empty, all modes included.
#' @param onlyCase If TRUE, only "case" predicted samples included. Default: TRUE. 
#' @param onlyPOS If TRUE, only "POS" labeled samples included. Default: FALSE. 
#' @param splitby Metadata by which to split plots. Can be "none", "verdict", or "condType".
#' @param returnData If TRUE, plot data is returned instead of plotting. Default: FALSE
#' @return A named list of ggplot objects.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
plotEnrichment <- function(object,
                           modes = NULL,
                           onlyCase = TRUE,
                           onlyPOS = FALSE,
                           splitby = c("none", "verdict", "condType"),
                           returnData = FALSE) {
  
  # TODO: Should there be an option to control this for the user?
  yval <- "stat_fisher_rl"
  
  # TODO: NEEDS to be replaced with RLHub
  rlbase <- "https://rlbase-data.s3.amazonaws.com"
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
  
  # Add metadata
  usamp <- object@metadata$sampleName
  sampleRes$experiment <- usamp
  sampleRes$condType <- object@metadata$condType
  sampleRes$mode <- object@metadata$mode
  predres <- rlresult(object, resultName = "predictRes")
  sampleRes$verdict <- predres$Verdict
  
  # Get min, max for plotting
  limit <- ifelse(yval == "stat_fisher_rl", 15, Inf)
  
  # Filter RLBase
  if (! is.null(modes)) {
    rlsamples <- dplyr::filter(rlsamples, .data$mode %in% {{ modes }})
  }
  if (onlyCase) {
    rlsamples <- dplyr::filter(rlsamples, .data$verdict == "Case")
  }
  if (onlyPOS) {
    rlsamples <- dplyr::filter(rlsamples, .data$condType == "POS")
  }

  # Wrangle the RLBase data
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
      .data$verdict, .data$condType, .data$mode,
      dplyr::contains("stat_fisher_rl")
    ) %>%
    dplyr::mutate(
      stat_fisher_rl = log2(.data$stat_fisher_rl)
    ) %>%
    dplyr::filter(
      !is.na(.data$stat_fisher_rl),
      is.finite(.data$stat_fisher_rl)
    ) %>%
    dplyr::filter(.data$experiment %in% c(rlsamples$rlsample, usamp))
  
  # Wrap strings
  input_data$type <- gsub(input_data$type, pattern = "_", replacement = " ")
  input_data$type[nchar(input_data$type) > 30] <-
    stringr::str_wrap(input_data$type[nchar(input_data$type) > 30], 30)

  # Make selected explicit
  input_data$selected <- ifelse(input_data$experiment == usamp, usamp, "")
  
  # Build plots
  datalst <- input_data %>%
    dplyr::group_by(.data$db) %>%
    {setNames(dplyr::group_split(.), dplyr::group_keys(.)[[1]])}
  
  # Return data if requested
  if (returnData) {
    return(datalst)
  }
  
  plts <- lapply(
    datalst,
    function(x) {
      db_now <- x$db[1]
      col <- aux$db_cols$col[aux$db_cols$db == db_now]
      
      if (! usamp %in% x$experiment) {
        return(NULL)
      }
      
      # Define limits
      lmts <- c(
        max(min(x[,yval]), -limit) - 1,
        min(max(x[,yval]), limit) + 1
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
        
      
      if (splitby == "none") {
        plt <- pltbase +
          ggplot2::geom_boxplot(
            color = "black",
            fill = col,
            position = ggplot2::position_dodge(.9),
            alpha = .5
          ) +
          ggplot2::geom_jitter(
            alpha=ifelse(x$experiment == usamp, 1, 0), 
            size = 3, 
            fill = col,
            width = 0,
            color = "black",
            shape=23, 
            stroke = 1.5
          )
      } else {
        # Get split cols
        cols <- aux[[paste0(tolower(splitby), "_cols")]]
        colvec <- cols %>% dplyr::pull(.data$col)
        names(colvec) <- cols %>% dplyr::pull({{ splitby }})
        
        # Add box and jitter
        plt <- pltbase + 
          ggplot2::geom_boxplot(
            color = "black",
            position = ggplot2::position_dodge(.9),
            alpha = 1
          ) +
          ggplot2::geom_jitter(
            alpha=ifelse(x$experiment == usamp, 1, 0), 
            size = 3, 
            color = "black",
            position = ggplot2::position_jitterdodge(
              dodge.width = .9,
              jitter.width = 0, 
              jitter.height = 0
            ),
            shape=23, 
            stroke = 1.5
          ) +
          ggplot2::scale_fill_manual(
            values = colvec, drop=TRUE
          )
      }
      
      # Final additions to plot and themeing
      plt +
        ggplot2::ylab(
          stringr::str_to_title(
            gsub(yval,
                 pattern = "_",
                 replacement = " "
            )
          )
        ) +
        ggplot2::xlab(NULL) +
        ggplot2::labs(
          title = gsub(
            db_now,
            pattern = "_",
            replacement = " "
          ),
          subtitle = usamp,
          caption = "\u25C7 - User sample"
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
          legend.title = ggplot2::element_text(size=18),
          legend.text = ggplot2::element_text(size=14),
          plot.caption = ggplot2::element_text(size=12)
        ) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(
            override.aes = list(size=0, stroke=0)
          )
        ) + 
        ggplot2::scale_y_continuous(limits = lmts)
      
    }
  )

  # Remove empty
  plts <- plts[!sapply(plts, is.null)]
  
  # Return
  return(plts)
}



#' Plot RL-Region overlap with RLRanges
#'
#' @param object An RLRanges object with \code{rlRegionTest()} already run.
#' @param returnData If TRUE, plot data is returned instead of plotting. Default: FALSE
#' @return A venn diagram ggplot object.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
plotRLRegionOverlap <- function(object, returnData=FALSE) {
  
  # Get RLRegions 
  # TODO: NEEDS to be in RLHub 
  rlbase <- "https://rlbase-data.s3.amazonaws.com"
  rlregions_table <- file.path(rlbase, "RLHub", "rlregions_table.rda")
  tmp <- tempfile()
  download.file(rlregions_table, destfile = tmp, quiet = TRUE)
  load(tmp)
  
  # Retrieve results
  olres <- rlresult(object, "rlRegionRes")
  
  # Get abstracted sites
  pkolnms <- olres$Overlap$name__peaks
  rlolnms <- olres$Overlap$name__rlregion
  
  # Get overlapping ranges
  pkolgr <- object[names(object) %in% pkolnms,]
  
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
  pkonly <- length(object[! names(object) %in% pkolnms,])
  rlonly <- rlReg %>% dplyr::filter(! .data$name %in% {{ rlolnms }}) %>% nrow()
  
  # Build the tbl and plot
  sites <- seq(pkonly + rlonly + shared)
  
  # Get the plot data
  pltdata <- list(
    sites[c(rep(FALSE, pkonly), rep(TRUE, rlonly), rep(TRUE, shared))],
    sites[c(rep(TRUE, pkonly), rep(FALSE, rlonly), rep(TRUE, shared))]
  ) %>%
    setNames(nm = c("RL-Regions", object@metadata$sampleName))
  if (returnData) {
    return(pltdata)
  }
  
  # Make the plot
  futile.logger::flog.threshold(futile.logger::ERROR,
                                name = "VennDiagramLogger")
  gt <- pltdata %>%
    VennDiagram::venn.diagram(filename = NULL,
                              fill = c("#9ad9ab", "#9aa0d9"),
                              main = "RL-Region Overlap",
                              sub = paste0(
                                "Fisher's Test: pval = ", 
                                signif(olres$Test_results$p.value),
                                "; odds ratio = ", 
                                signif(olres$Test_results$estimate)
                              ),
                              main.cex = 2,
                              cat.pos = c(200, 160),
                              cat.dist = c(0.05, 0.05),
                              margin = .05, sub.pos = c(0.5,0.15)) %>%
    grid::grid.draw() %>%
    ggplotify::grid2grob() %>% 
    ggplotify::as.ggplot()
  
  # Return
  return(gt)
  
}

