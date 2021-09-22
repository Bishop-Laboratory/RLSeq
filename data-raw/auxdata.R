#' This script is used to generate auxdatailiary information used by RLSeq and RLBase.
#' This includes color pallets, metadata, and other small items.
#' The result is a named list which RLSeq depends upon.

library(tidyverse)

### Metadata ###

## RLBase samples ##

rlbase <- "https://rlbase-data.s3.amazonaws.com"

# Load the RLBase samples
rlbase_enrich <- file.path(rlbase, "RLHub", "rlsamples.rda")
tmp <- tempfile()
download.file(rlbase_enrich, destfile = tmp, quiet = TRUE)
load(tmp)

## Annotations ##

# Load the annotations -- add family info
rlbase_anno <- file.path(rlbase, "RLHub", "annotations_all_hg38.rda")
tmp <- tempfile()
download.file(rlbase_anno, destfile = tmp, quiet = TRUE)
load(tmp)

# Load RLBP info
rlbps <- file.path(rlbase, "RLHub", "rlbps.rda")
tmp <- tempfile()
download.file(rlbps, destfile = tmp, quiet = TRUE)
load(tmp)

# Get the available genomes
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
    dplyr::select(mode, family, ip_type, strand_specific, moeity, bisulfite_seq) %>%
    distinct()


# Get the databases and add order/group information
pat <- "(.+)__(.+)"
annotypes <- tibble(
    db = gsub(names(annotations_all), pattern = pat, replacement = "\\1"),
    type = gsub(names(annotations_all), pattern = pat, replacement = "\\2")
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


## S9.6 vs dRNaseH1 vs Other ##
ip_cols <- tribble(
    ~ip_type, ~POS,
    "S9.6", "#2E294E",
    "dRNH", "#C1292E",
    "Other", "#F3B700"
) %>%
    mutate(
        "NEG" = colorspace::desaturate(POS, amount = .6),
        "NULL" = colorspace::desaturate(POS, amount = 1),
        "NULL" = case_when(
            is.na(`NULL`) ~ colorspace::desaturate(POS, amount = 0.75),
            TRUE ~ `NULL`
        )
    )

## Modes within each
# From https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 110)[1:n]
}
modes <- unique(rlsamples$mode)
set.seed(5)
mode_cols <- tibble(
    mode = sample(modes, size = length(modes)),
    col = gg_color_hue(
        length(
            modes
        )
    )
)

## CondType
condtype_cols <- tribble(
    ~condType, ~col,
    "POS", "#dedee0",
    "NEG", "#2e2e63",
    "NULL", "#05052b"
)

## Verdict
verdict_cols <- tribble(
    ~verdict, ~col,
    "Case", "#e0dede",
    "Control", "#8a2c2c"
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
    condtype_cols = condtype_cols,
    verdict_cols = verdict_cols,
    available_modes = available_modes,
    available_genomes = genomes
)
usethis::use_data(auxdata, compress = "xz", overwrite = TRUE)