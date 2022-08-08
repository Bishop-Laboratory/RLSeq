#' Script for creating the rlbaseNoiseAnalyze.rda object
#' Assumes that you have followed the steps to create "RLBase-data"
#' https://github.com/Bishop-Laboratory/RLBase-data
#' In a separate directory at "../RLBase-data".

## Get the stats for all samples
rlsamps <- RLHub::rlbase_samples()
rlsampsf <- rlsamps %>%
    dplyr::filter(genome %in% c("hg38", "mm10"))
fls <- list.files("../RLBase-data/misc-data/rlranges/", full.names = T)
reslst <- parallel::mclapply(
    seq(nrow(rlsampsf)), function(i) {
        rlr <- fls[grep(fls, pattern = rlsampsf$rlsample[i])]
        if (!length(rlr)) {
            return(NULL)
        }
        rlr <- readRDS(rlr)

        rlres <- rlr@metadata$results
        if (is.null(attr(rlres, "noiseAnalysis"))) {
            attr(rlres, "noiseAnalysis") <- dplyr::tibble()
            rlr@metadata$results <- rlres
            rlr <- noiseAnalyze(rlr)
        }

        rlr <- try(noiseAnalyze(rlr))

        if ("try-error" %in% class(rlr)) {
            return(NULL)
        }

        mean(rlr@metadata$results@noiseAnalysis$value)
    },
    mc.cores = 44
)
names(reslst) <- rlsampsf$rlsample

# Remove nulls
reslstf <- reslst[sapply(reslst, function(x) !is.null(x))]

# Create tbl representation
rlbaseNoiseAnalyze <- reslstf %>%
    lapply(dplyr::as_tibble) %>%
    dplyr::bind_rows(.id = "rlsample")

# Save
usethis::use_data(rlbaseNoiseAnalyze, overwrite = TRUE)
