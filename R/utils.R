#' Check if URL exists
#' @param url URL to check
#' @return logical. TRUE if status code 200, FALSE if not
urlExists <- function(url) {
  identical(
    httr::status_code(
      # Checks HEAD only due to size constraints
      httr::HEAD(
        url
      )
    ), 200L  # Checks if response is ok
  )
}


#' Get Chrom Sizes
#' Helper function which downloads chrom sizes from UCSC for a genome.
#' @param genome the UCSC genome for which to download chrom sizes
#' @return A tibble containing chrom sizes
getChromSizes <- function(genome) {
  chrom_sizes <- suppressMessages(readr::read_tsv(paste0(
    BASE_UCSC,
    genome, '/bigZips/', genome, '.chrom.sizes'
  ), col_names = FALSE))
  return(chrom_sizes)
}


#' Check RLFS Anno
#' Helper function that checks whether a genome has RLFS available
#' @param genome the UCSC genome name to check
#' @return A logical, TRUE if available, FALSE if not
checkRLFSAnno <- function(genome) {
  return(
    urlExists(
      paste0(
        RLFS_BED_URL, 
        genome, ".rlfs.bed"
      )
    )
  )
}


#' Get RLFS Anno
#' Helper function that retrieves RLFS
#' @param genome the UCSC genome name to retrieve RLFS for
#' @return A GRanges object with RLFS for that species.
getRLFSAnno <- function(genome) {
  
  # Check if annotations available first
  if (! checkRLFSAnno(genome)) {
    stop("No RLFS annotations available for ", genome)
  }
  
  # Return as a GRanges object
  return(
    regioneR::toGRanges(
      as.data.frame(
        suppressMessages(readr::read_tsv(
          paste0(
            RLFS_BED_URL, 
            genome, ".rlfs.bed"
          ),
          col_names = FALSE))
      )
    )
  )
}


#' Get Chain
#' Helper function that retrieves chain file for liftUtil()
#' @param genomeFrom the UCSC genome name to convert from.
#' @param genomeTo the UCSC genome name to convert to.
getChain <- function(genomeFrom, genomeTo) {
  
  # Get URL
  url <- paste0(BASE_UCSC, genomeFrom, "/liftOver/", 
                genomeFrom, "To",
                paste0(toupper(substring(genomeTo, 1, 1)), 
                       substring(genomeTo, 2)),
                ".over.chain.gz")
  
  # Check if exists
  stopifnot(urlExists(url))
  
  # Check if R.utils available
  if( ! requireNamespace("R.utils", quietly = TRUE)) {
    stop("R.utils is required. Please install it with install.packages('R.utils')")
  }
  
  # Get the chain
  tmp <- tempfile()
  download.file(url, destfile = paste0(tmp, ".gz"))
  R.utils::gunzip(paste0(tmp, ".gz"))
  chain <- rtracklayer::import.chain(tmp)
  
  # Return as a GRanges object
  return(
    chain
  )
}


#' Lift Over Utility
#'
#' Convenience function for converting ranges from between assemblies
#'
#' @param ranges A GRanges object in hg19
#' @param genomeFrom Genome of ranges supplied, in UCSC format (e.g., "hg19")
#' @param genomeTo Genome to convert to (e.g., "hg38")
#' @return A lifted GRanges object
#' @examples
#' 
#' hg38Lift(RSeqR::SRX1025890_peaks_hg19)
#' 
#' @export
liftUtil <- function(ranges, genomeFrom, genomeTo) {
  
  # Get the chain
  chain <- getChain(genomeFrom, genomeTo)
  
  # Lift Over
  lifted <- rtracklayer::liftOver(ranges, chain = chain) %>%
    unlist() 
  
  # Force uniqueness
  nms <- names(lifted) %>% unique()
  lifted <- lifted[nms,]
  
  return(lifted)
}

