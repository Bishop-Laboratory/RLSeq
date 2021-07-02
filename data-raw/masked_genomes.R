library(tidyverse)

# Get the list of genome masks if not already available
if (! file.exists("data-raw/maskLst.rda")) {
  maskLst <- RSeqR:::buildGenomeMasks()
  save(maskLst, file = "misc/maskLst.rda", compress = "xz")
} else {
  load("data-raw/maskLst.rda")
}

# Shrink size of mask genome list
# 1. Only keep masks which have RLFS
maskLst <- maskLst[sapply(names(maskLst), RSeqR:::checkRLFSAnno)]
# 2. Remove large masks
maskLst <- maskLst[sapply(maskLst, function(x) {
  object.size(x) < 500000
})]
names(maskLst) <- paste0(names(maskLst), ".masked")
# 3. Remove ranges of size 1
genomeMasks <- lapply(names(maskLst), function(genome) {
  newGr <- maskLst[[genome]]
  newGr[GenomicRanges::width(newGr) > 1]
})
names(genomeMasks) <- names(maskLst)

# Save the data with xz compression
usethis::use_data(genomeMasks, compress = "xz", overwrite = TRUE)

