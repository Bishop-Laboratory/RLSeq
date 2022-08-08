#' Generates the random genomic windows needed by noiseAnalyze

library(dplyr)

# List of ENCODE black lists for each genome available
bllst_lst <- list(
    "hg19" = "https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz",
    "hg38" = "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz",
    "mm10" = "https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"
)

genomes <- c("mm10", "hg38", "hg19")
randomWindows <- lapply(
    genomes, function(gen) {

        # Get chrom.sizes
        chrompath <- paste0(
            "https://hgdownload.cse.ucsc.edu/goldenpath/",
            gen, "/bigZips/", gen, ".chrom.sizes"
        )
        chromsiz <- valr::read_genome(chrompath) %>% filter(size > 1E7)

        # windows 1000 random genomic windows from the bed file
        windows <- valr::bed_random(
            genome = chromsiz, seed = 42, n = 1000
        )

        # Get the black list
        blst <- valr::read_bed(bllst_lst[[gen]])

        # Remove black list sites
        windows <- valr::bed_subtract(x = windows, blst)

        # Return result
        return(windows)
    }
)

# Set names as genomes
names(randomWindows) <- genomes

# Save
usethis::use_data(randomWindows, overwrite = TRUE)
