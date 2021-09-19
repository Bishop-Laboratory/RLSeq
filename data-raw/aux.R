#' This script is used to generate auxiliary information used by RLSeq and RLBase.
#' This includes color pallets, metadata, and other small items.
#' The result is a named list which RLSeq depends upon.

library(tidyverse)

### Metadata ###

## RLBase samples ##

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
full_cols <- c(
  "#C12930", "#8D2840", "#9B416E", "#77518F", "#5D3C9D", "#354980",
  "#385D37", "#2E294E", "#501B35", "#8A5D54", "#413284", "#09360D",
  "#592566", "#F3B80C", "#B8934E", "#CB886D", "#BDB0D5", "#86C1B9",
  "#A94548", "#6C557F", "#463852", "#224431"
)

set.seed(42)
mode_cols <- tibble(
  mode = unique(rlsamples$mode),
  POS = sample(full_cols, length(unique(rlsamples$mode)))
) %>%
  mutate(
    "NEG" = colorspace::desaturate(POS, amount = .6),
    "NULL" = colorspace::desaturate(POS, amount = 1),
    "NULL" = case_when(
      is.na(`NULL`) ~ colorspace::desaturate(POS, amount = 0.75),
      TRUE ~ `NULL`
    )
  )

## CondType
condtype_cols <- tribble(
  ~condType, ~col,
  "POS", "#41919D",
  "NEG", "#6C8C91",
  "NULL", "#868686"
)

## Verdict
verdict_cols <- tribble(
  ~verdict, ~col,
  "Case", "#C69952",
  "Control", "#B59C7F"
)

## Save
aux <- list(
  db_cols = db_cols,
  annotypes = annotypes,
  rlbps = rlbps,
  ip_cols = ip_cols,
  mode_cols = mode_cols,
  condtype_cols = condtype_cols,
  verdict_cols = verdict_cols
)
usethis::use_data(aux, compress = "xz")
