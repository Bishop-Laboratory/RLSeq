# Script used on 8/25/2021 to get the rlRegions data object
RLTAB_LOC <- "~/projects/RMapDB-shiny/data/rltab.rda"
NEW_RLTAB <- "data-raw/rltab.rda"

# Move rlTab and load it
file.copy(RLTAB_LOC, to = NEW_RLTAB)
load(NEW_RLTAB)

# Clean it up for minimal file size
rltabShow %>%
  dplyr::select(
    `RL Region`, Location, `# of Studies`, `# of Samples`, Modes,
    `Mean Signal`, `Mean FDR`
  ) %>%
  object.size()/1E6

