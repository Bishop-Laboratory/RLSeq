toBed <- function(granges, write = FALSE, filename = NULL) {
  df <- as.data.frame(granges)
  names(df)[1] <- paste0("#", names(df)[1])
  if(write) {
    write.table(df, file = paste0(filename, ".bed"), sep = "\t", col.names = NA)
  }
  return(df)
}
