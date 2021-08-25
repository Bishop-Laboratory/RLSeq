#' Check if URL exists
#' @param url URL to check
#' @return logical. TRUE if status code 200, FALSE if not
#' @export
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
#' @export
getChromSizes <- function(genome) {
  chrom_sizes <- suppressMessages(readr::read_tsv(paste0(
    'http://hgdownload.soe.ucsc.edu/goldenPath/',
    genome, '/bigZips/', genome, '.chrom.sizes'
  ), col_names = FALSE))
  return(chrom_sizes)
}


#' Check RLFS Anno
#' Helper function that checks whether a genome has RLFS available
#' @param genome the UCSC genome name to check
#' @return A logical, TRUE if available, FALSE if not
#' @export
checkRLFSAnno <- function(genome) {
  return(
    urlExists(
      paste0(
        "https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/", 
        genome, ".rlfs.bed"
      )
    )
  )
}


#' Get RLFS Anno
#' Helper function that retrieves RLFS
#' @param genome the UCSC genome name to retrieve RLFS for
#' @return A GRanges object with RLFS for that species.
#' @export
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
            "https://rmapdb-data.s3.us-east-2.amazonaws.com/rlfs-beds/", 
            genome, ".rlfs.bed"
          ),
          col_names = FALSE))
      )
    )
  )
}


#' Check Gene Annotations
#' Helper function which checks the chrom sizes info from UCSC
#' @param genome the UCSC genome for which to download chrom sizes
#' @return logical. TRUE if status code 200, FALSE if not
#' @export
checkGenes <- function(genome) {
  return(
    urlExists(
      paste0(
        'http://hgdownload.soe.ucsc.edu/goldenPath/',
        genome, '/bigZips/genes/', genome, '.ensGene.gtf.gz'
      )
    )
  )
}


