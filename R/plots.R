
rmapHeatmap <- function(corrRes, rmapSamples, prediction=NA, cleanAnno=FALSE, selected=NA) {
  
  # Wrangle the annotation data
  annoCorr <- RLSeq::rmapSamps %>%
    dplyr::mutate(group = "RMapDB") %>%
    dplyr::select(.data$id, .data$mode, .data$is_rnh_like,
                  .data$prediction, .data$group) %>%
    dplyr::bind_rows(
      tibble::tibble(
        id="user_supplied",
        mode=NA,
        is_rnh_like=NA,
        prediction=prediction,
        group="User-supplied"
      )
    ) %>%
    column_to_rownames(var="id")
  
  # Filter for available / desired samples
  toSelect <- colnames(corrRes)
  if (! is.na(selected)) {
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
      annoCorrNow <- annoCorrNow[,-which(colnames(annoCorrNow) == "is_ctrl")]
    }
    
    # Select pred_ctrl 
    if (any(annoCorrNow$pred_ctrl)) {
      annoCorrNow$pred_ctrl <- as.factor(annoCorrNow$pred_ctrl)
    } else {
      annoCorrNow <- annoCorrNow[,-which(colnames(annoCorrNow) == "pred_ctrl")]
    }
  }
  
  # Pallete
  paletteLength <- 100
  myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(corrNow), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(corrNow)/paletteLength, max(corrNow), length.out=floor(paletteLength/2)))
  pheatColLst <- colList[which(names(colList) %in% colnames(annoCorrNow))]
  pheatColLst$Mode <- colList$mode[which(names(colList$mode) %in% annoCorrNow$Mode)]
  pheatmap::pheatmap(corrRes, show_rownames = FALSE, 
                     main = paste0("RMapDB Corr Heatmap\n",
                                   current_samp()
                     ), fontsize = 13.5,
                     show_colnames = FALSE, silent = TRUE,
                     annotation_colors = pheatColLst, 
                     color = myColor, breaks = myBreaks,
                     annotation_col = annoCorrNow) %>%
    pluck(4) %>%
    ggplotify::as.ggplot()
}


