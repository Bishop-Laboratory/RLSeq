library(tidyverse)

maskLst <- RSeqR:::buildGenomeMasks()

load("misc/maskLst.rda")
save(maskLst, file = "misc/maskLst2.rda", compress = "xz")

dir.create("misc/mask_genomes", showWarnings = F)
lapply(names(maskLst), function(genome) {
  outFile <- paste0("misc/mask_genomes/", genome, ".masked.bed")
  
  maskLst[[genome]] %>%
    as.data.frame() %>%
    dplyr::select(seqnames, start, end) %>%
    write_tsv(file = outFile,
              col_names = FALSE)
  
  system(paste0("xz ", outFile))
  
  
  
  save(dd, file = "misc/mask_genomes/hg38.masked.rda", compress = "xz")
  
  
  
  ChIPpeakAnno::toGRanges("misc/mask_genomes/hg38.masked.bed.xz")
  
  
  
})


