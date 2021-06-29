#' Build Available Genomes
#' Helper Function which builds the UCSC available genomes dataset
#' @param test if TRUE, will only build info set for first 2 genomes in UCSC.
#' @param ... arguments passed to `getEffectiveGenomeSizes()`
#' @return A data frame of avaialable genomes and corresponding info. 
#' @importFrom dplyr %>%
#' @importFrom rlang .data
buildAvailableGenomes <- function(test=FALSE, ...) {
  
  envName <- buildCondaEnv(packages = "khmer", 
                           envName = "khmerEnv",
                           channel = c("bioconda", "conda-forge"))
  
  # from http://rstudio-pubs-static.s3.amazonaws.com/562103_092b7264f392482e827440cf1521363c.html
  api_genome <- restfulr::RestUri("http://api.genome.ucsc.edu/list")
  response <- restfulr::read(api_genome$ucscGenomes)[['ucscGenomes']]
  
  # Wrangle the data from UCSC Genome API
  available_genomes <- do.call(rbind.data.frame, response) %>%
    tibble::rownames_to_column(var = "UCSC_orgID") %>%
    dplyr::group_by(UCSC_orgID) %>%
    dplyr::mutate(genes_available = checkGenes(UCSC_orgID),
                  homer_anno_available = checkHomer(UCSC_orgID),
                  year = as.numeric(gsub(x = description, 
                                         replacement = "\\1",
                                         pattern = ".+ ([0-9]+) \\(.+")))
  # Get the effective genome sizes
  eff_gen_sizes <- available_genomes %>%
    dplyr::pull(UCSC_orgID) %>%
    lapply(function(genome) {
      getEffectiveGenomeSizes(genome = genome, envPath = envName, ...)
    }) %>%
    dplyr::bind_rows()
  
  # Bind final table and return
  return(
    dplyr::left_join(available_genomes, eff_gen_sizes, by = "UCSC_orgID")
  )
}

# RSeqR:::buildAvailableGenomes()

