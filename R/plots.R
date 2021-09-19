#' Plot Perm Test results
#'
#' @description Simple function for plotting permTest
#' results using ggplot2.
#'
#' @param rlfsRes The list object generated from running \code{analyzeRLFS}.
#' @param ... Additional parameters passed to \code{ggplot}.
#' @return A ggplot object.
#' @export
plotRLFSRes <- function(rlfsRes,
                        ...) {

  # Validate rlfsRes
  stopifnot(is.list(rlfsRes))

  # Obtain the plot data
  pltdat <- tibble::tibble(
    "zscore" = rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifted.z.scores,
    "shift" = rlfsRes$`Z-scores`$`regioneR::numOverlaps`$shifts
  )
  pval <- rlfsRes$perTestResults$`regioneR::numOverlaps`$pval

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
    ggplot2::ggtitle("ZScore around RLFS") +
    ggplot2::labs(caption = paste0("p < ", round(pval, digits = 5))) +
    ggplot2::ylab("Peak Enrichment (Z-Score)") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::xlab("Distance to RLFS (bp)") +
    ggprism::theme_prism(base_size = 15)
}


#' RLBase Heatmap
rmapHeatmap <- function(corrRes, rmapSamples, prediction = NA, cleanAnno = FALSE, selected = NA) {

  # Wrangle the annotation data
  annoCorr <- RLSeq::rmapSamps %>%
    dplyr::mutate(group = "RMapDB") %>%
    dplyr::select(
      .data$id, .data$mode, .data$is_rnh_like,
      .data$prediction, .data$group
    ) %>%
    dplyr::bind_rows(
      tibble::tibble(
        id = "user_supplied",
        mode = NA,
        is_rnh_like = NA,
        prediction = prediction,
        group = "User-supplied"
      )
    ) %>%
    column_to_rownames(var = "id")

  # Filter for available / desired samples
  toSelect <- colnames(corrRes)
  if (!is.na(selected)) {
    toSelect <- intersect(selected, toSelect)
  }
  corrNow <- corrRes[toSelect, toSelect]
  annoNow <- annoCorr[toSelect, ]

  # Match up columns if desired (for Shiny implementation mostly)
  if (cleanAnno) {
    # Select isControl if there's a good reason to
    if (any(annoCorrNow$is_ctrl)) {
      annoCorrNow$is_ctrl <- as.factor(annoCorrNow$is_ctrl)
    } else {
      annoCorrNow <- annoCorrNow[, -which(colnames(annoCorrNow) == "is_ctrl")]
    }

    # Select pred_ctrl
    if (any(annoCorrNow$pred_ctrl)) {
      annoCorrNow$pred_ctrl <- as.factor(annoCorrNow$pred_ctrl)
    } else {
      annoCorrNow <- annoCorrNow[, -which(colnames(annoCorrNow) == "pred_ctrl")]
    }
  }

  # Pallete
  paletteLength <- 100
  myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(
    seq(min(corrNow), 0, length.out = ceiling(paletteLength / 2) + 1),
    seq(max(corrNow) / paletteLength, max(corrNow), length.out = floor(paletteLength / 2))
  )
  pheatColLst <- colList[which(names(colList) %in% colnames(annoCorrNow))]
  pheatColLst$Mode <- colList$mode[which(names(colList$mode) %in% annoCorrNow$Mode)]
  pheatmap::pheatmap(corrRes,
    show_rownames = FALSE,
    main = paste0(
      "RMapDB Corr Heatmap\n",
      current_samp()
    ), fontsize = 13.5,
    show_colnames = FALSE, silent = TRUE,
    annotation_colors = pheatColLst,
    color = myColor, breaks = myBreaks,
    annotation_col = annoCorrNow
  ) %>%
    pluck(4) %>%
    ggplotify::as.ggplot()
}



#' Plot Enrichment Test Results
#'
#' @description Plotting function for featureEnrich() results.
#'
#' @param sampleRes The tibble obejct obtained from running \code{featureEnrich}.
#' @param rlbaseRes A tibble containing results from RLBase samples.
#' @param rlsamples A tibble containing the metadata of RLBase samples.
#' @param yval Column to use for obtaining Y-axis values. Default: "stat_fisher_rl".
#' If NULL, no comparison with RLBase will be made.
#' @return A named list of ggplot objects.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
plotEnrichment <- function(sampleRes,
                           rlbaseRes = NULL,
                           rlsamples = NULL,
                           splitby = c("none", "verdict", "condtype"),
                           facetby = c("none", "mode"),
                           yval = "stat_fisher_rl") {
  limit <- ifelse(yval == "stat_fisher_rl", 10, Inf)

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
  plts <- input_data %>%
    dplyr::group_by(.data$db) %>%
    dplyr::group_split() %>%
    lapply(
      ., function(x) {
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
          xlab(NULL) +
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

  suppressWarnings(lapply(unique(anno_data$annotate_type), function(annoNow) {
    toPlt <- anno_data %>%
      rename(is_ctrl = is_rnh_like) %>%
      mutate(pred_ctrl = prediction == "Control") %>%
      filter(
        annotate_type == !!annoNow,
        id %in% rmapSampsRV()
      ) %>%
      mutate(Annotation = if (!!annoNow == "gene") {
        factor(Annotation, levels = annoPlot_genelvls)
      } else {
        Annotation
      }) %>%
      mutate(sample_now = id == !!current_samp()) %>%
      fill(everything(0)) %>%
      arrange(sample_now)

    if (input$splitAnnoBy == "None") {
      plt <- ggplot(toPlt, mapping = aes(
        x = Annotation,
        y = !!sym(opt)
      )) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_violin(
          # draw_quantiles = c(.50),
          trim = FALSE,  position = position_dodge(.9),
          fill = vCols[[annoNow]]
        ) +
        geom_boxplot(
          width = .12, color = "black", position = position_dodge(.9),
          fill = boxCols[[annoNow]],
          alpha = 1
        ) +
        geom_point(
          alpha = ifelse(toPlt$sample_now, 1, 0), size = 4, color = "black",
          shape = 23, stroke = 2, fill = "#32889c"
        )
    } else {
      if (input$splitAnnoBy != "mode") {
        plt <- ggplot(toPlt, mapping = aes(
          x = Annotation,
          fill = !!sym(input$splitAnnoBy),
          y = !!sym(opt)
        )) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_violin(
            trim = FALSE, position = position_dodge(.9)
          )
      } else {
        plt <- ggplot(toPlt, mapping = aes(
          x = Annotation,
          fill = !!sym(input$splitAnnoBy),
          y = !!sym(opt)
        )) +
          geom_hline(yintercept = 0, linetype = "dashed")
      }

      # So that only the colors in the factor level appear in the legend
      cols <- annoFillSplit[[input$splitAnnoBy]][unique(as.character(as.data.frame(toPlt)[, input$splitAnnoBy]))]

      # Add box and jitter
      plt <- plt + geom_boxplot(
        width = ifelse(input$splitAnnoBy != "mode", .12, .85),
        color = "black",
        position = position_dodge(.9),
        alpha = 1
      ) +
        geom_jitter(
          alpha = ifelse(toPlt$sample_now, 1, 0),
          size = 4, color = "black",
          position = position_jitterdodge(
            dodge.width = .9,
            jitter.width = 0,
            jitter.height = 0
          ),
          shape = 23, stroke = 2
        ) +
        scale_fill_manual(
          values = cols, drop = TRUE
        )
    }

    plt <- plt +
      ggpubr::rremove("legend") +
      ylab(opt) +
      xlab(NULL) +
      labs(title = annoPlot_titles[annoNow]) +
      theme_prism(base_size = 16) +
      theme(axis.text.x = element_text(
        angle = 45, vjust = 1,
        hjust = 1
      )) +
      theme(
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14)
      )
  })) %>%
    ggpubr::ggarrange(plotlist = ., nrow = 3, ncol = 1, align = "hv")
}
