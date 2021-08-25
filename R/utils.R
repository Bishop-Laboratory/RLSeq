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


#' Check Homer Annotations
#' Helper function which checks whether homer annotations are available
#' @param genome the UCSC genome for which to download chrom sizes
#' @return logical. TRUE if status code 200, FALSE if not
#' @export
checkHomer <- function(genome) {
  return(
    urlExists(
      paste0("http://homer.ucsd.edu/homer/data/genomes/",
             genome, ".v6.4.zip")
    )
  )
}


#' Get Homer Anno
#' Helper function that retrieves Homer Annotations
#' @param genome the UCSC genome name to retrieve annotations for
#' @return A list object with annotations for that species.
#' @export
getHomerAnno <- function(genome) {
  
  # Check if annotations available first
  if (! checkHomer(genome)) {
    stop("No Homer annotations available for ", genome)
  }
  
  # Get the file
  download.file(paste0("http://homer.ucsd.edu/homer/data/genomes/",
                       genome, ".v6.4.zip"), 
                destfile = "homer.zip")
  
  # Return as a GRanges object
  dd <- read_tsv("data/genomes/hg38/hg38.full.annotation", 
                 col_names = c("name", "seqnames", "start", "end", "strand", "type", "unknown"))
  dd2 <- dd %>%
    select(-name, -unknown) %>%
    filter(! grepl(type, pattern = "\\?", perl = TRUE)) %>%
    group_by(type) %>% 
    {setNames(group_split(.), group_keys(.)[[1]])}
    
  
  
  codeKey <- tibble(
    "type" = c("I", "N", "P", "E"),
    "val" = c("Intron", "Intergenic", "Promoter", "Exon")
  )
    
  
}