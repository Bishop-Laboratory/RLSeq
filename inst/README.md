# `inst/` directory

The items in this folder are installed along with `RLSeq`. Briefly:

1. The `Rmd/` directory contains the RMarkdown templates necessary for running
the `report()` function. 
2. The `int-data/` contains objects which are made available to all `RLSeq()`
function via `.onload()` (see `zzz.R`).
3. The `ext-data/` contains objects which are used in examples within the package.
4. The `data-raw/` directory contains the functions run to generate all data
in the `ext-data/` and `int-data/` directories. 
