maskLst <- RSeqR:::buildGenomeMasks()


sesstoken <- system("aws sts get-session-token", intern = TRUE)

lapply(names(maskLst), function(x) {
  aws
})

load("misc/maskLst.rda")
