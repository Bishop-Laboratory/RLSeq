#' To Bed
#'
#' Converts a GRanges object to a .bed formatted DataFrame and writes to file if requested
#'
#' @param granges A GRanges object containing DRIP-Seq peaks
#' @param write A boolean determining if the converted DataFrame will be written to a .bed file
#' @param filename A string containing the desired file name if writing to file
#' @return A DataFrame object containing the GRanges content formatted according to .bed standards
#' @export

toBed <- function(granges, write = FALSE, filename = NULL) {
  df <- as.data.frame(granges)
  names(df)[1] <- paste0("#", names(df)[1])
  if(write) {
    write.table(df, file = paste0(filename, ".bed"), sep = "\t", col.names = NA)
  }
  return(df)
}
