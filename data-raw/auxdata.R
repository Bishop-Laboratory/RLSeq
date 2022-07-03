#' This script is used to generate auxdatailiary information used by
#' RLSeq and RLBase.
#' This includes color pallets, metadata, and other small items.
#' The result is a named list which RLSeq depends upon.

library(tidyverse)

### Metadata ###

## RLBase samples ##

# Load the RLBase samples
rlsamples <- RLHub::rlbase_samples()

## Annotations ##

# Load the annotations -- add family info
library(RLSeq)
library(RLHub)

rlbase_anno <- RLHub::annots_full_hg38()

# Load RLBP info
rlbps <- RLHub::rlbps()

# Get the available genomes
rlbase <- "https://rlbase-data.s3.amazonaws.com"
avgs <- file.path(rlbase, "misc", "available_genomes.rda")
tmp <- tempfile()
download.file(avgs, destfile = tmp, quiet = TRUE)
load(tmp)
genomes <- available_genomes %>%
    as_tibble() %>%
    dplyr::rename(rlfs_available = homer_anno_available) %>%
    dplyr::filter(rlfs_available, genes_available) %>%
    pull(UCSC_orgID)

# Get the available modes (no bisulfite currently supported)
available_modes <- rlsamples %>%
    dplyr::select(
        mode, family, ip_type, strand_specific, moeity,
        bisulfite_seq
    ) %>%
    distinct()


# Get the databases and add order/group information
pat <- "(.+)__(.+)"
annotypes <- tibble(
    db = gsub(names(rlbase_anno), pattern = pat, replacement = "\\1"),
    type = gsub(names(rlbase_anno), pattern = pat, replacement = "\\2")
)

### Color Palettes ###

## Annotation databases ##
pal_db <- c(
    "#BA3D5D", "#BD317E", "#BA2C9A", "#AD34B0", "#9546C0",
    "#6A5AC7", "#006CC4", "#007AB8", "#0082A6", "#00878E",
    "#008872", "#008650", "#00831F", "#007E00", "#4D7800",
    "#6F7100", "#876800", "#9A5E00", "#A95300", "#B44741"
)
set.seed(42)
db_cols <- tibble(
    db = unique(annotypes$db),
    col = sample(pal_db, length(unique(annotypes$db)))
)

## Prediction_Label colors ##
cols <- dplyr::tribble(
    ~cond, ~col,
    "NEG_NEG", "#8a2c2c",
    "NEG_POS", "#A76767",
    "POS_NEG", "#7aa4c4",
    "POS_POS", "#2270ab"
)
pred_lab_cols <- cols %>% dplyr::pull(.data$col)
names(pred_lab_cols) <- as.data.frame(cols)[, "cond"]

## S9.6 vs dRNaseH1 vs Other ##
ip_cols <- tribble(
    ~ip_type, ~POS,
    "S9.6", "#2E294E",
    "dRNH", "#C1292E",
    "Other", "#F3B700"
) %>%
    mutate(
        "NEG" = colorspace::desaturate(POS, amount = .6)
    )

## Modes within each
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 110)[1:n]
}
modes <- unique(rlsamples$mode)
misc <- names(which(table(rlsamples$mode) <= 12))
set.seed(6)
mode_cols <- tibble(
    mode = sample(modes[!modes %in% misc], size = length(modes[!modes %in% misc])),
    col = gg_color_hue(
        length(
            modes[!modes %in% misc]
        )
    )
)
misc_cols <- tibble(
    mode = "misc",
    col = "#8f8f8f"
)
mode_cols <- bind_rows(mode_cols, misc_cols)

## label
label_cols <- tribble(
    ~label, ~col,
    "POS", "#e0dede",
    "NEG", "#8a2c2c"
)

## prediction
prediction_cols <- tribble(
    ~prediction, ~col,
    "POS", "#e0dede",
    "NEG", "#8a2c2c"
)

## Heatcols
heatcols <- tribble(
    ~selected, ~col,
    "RLBase", "#e0dede",
    "user_selected", "#2c408a"
)

## Save
auxdata <- list(
    db_cols = db_cols,
    annotypes = annotypes,
    rlbps = rlbps,
    ip_cols = ip_cols,
    mode_cols = mode_cols,
    heat_cols = heatcols,
    label_cols = label_cols,
    prediction_cols = prediction_cols,
    prediction_label_cols = pred_lab_cols,
    available_modes = available_modes,
    available_genomes = genomes,
    misc_modes = misc
)
usethis::use_data(auxdata, overwrite = TRUE, compress = "xz")
