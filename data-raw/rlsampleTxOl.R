#' Script for creating the rlsampleTxOl.rda object
#' Assumes that you have followed the steps to create "RLBase-data"
#' https://github.com/Bishop-Laboratory/RLBase-data
#' In a separate directory at "../RLBase-data".


## Get the transcript features for hg38
anno <- RLHub::annots_primary_hg38()
txfeat <- anno[sapply(names(anno), grepl, pattern = "Transcript_Features")]
txfeat %>%
    bind_rows(.id = "feature") %>%
    mutate(
        feature = gsub(feature, pattern = ".+__", replacement = "")
    ) -> feats_hg38
## Get the transcript features for mm10
anno <- RLHub::annots_primary_mm10()
txfeat <- anno[sapply(names(anno), grepl, pattern = "Transcript_Features")]
txfeat %>%
    bind_rows(.id = "feature") %>%
    mutate(
        feature = gsub(feature, pattern = ".+__", replacement = "")
    ) -> feats_mm10

## Get the stats for all samples
rlsamps <- RLHub::rlbase_samples()
rlsampsf <- rlsamps %>%
    filter(genome %in% c("hg38", "mm10"))
fls <- list.files("../RLBase-data/rlbase-data/rlpipes-out/peaks/", pattern = "\\.broadPeak$", full.names = T)
reslst <- parallel::mclapply(
    seq(nrow(rlsampsf)), function(i) {
        bed <- fls[grep(fls, pattern = rlsampsf$rlsample[i])]
        if (!length(bed)) {
            return(NULL)
        }
        rlbed <- try(
            bed %>% valr::read_broadpeak()
        )
        if ("try-error" %in% class(rlbed)) {
            return(NULL)
        }

        feats <- if (rlsampsf$genome[i] == "hg38") feats_hg38 else feats_mm10
        olfeats <- valr::bed_intersect(rlbed, feats)
        # Get summarization
        olsum <- olfeats %>%
            group_by(name.x) %>%
            summarise(
                feats = paste0(feature.y, collapse = ", ")
            ) %>%
            mutate(
                feature = case_when(
                    grepl(feats, pattern = "TSS") ~ "TSS",
                    grepl(feats, pattern = "TTS") ~ "TTS",
                    grepl(feats, pattern = "fiveUTR") ~ "fiveUTR",
                    grepl(feats, pattern = "threeUTR") ~ "threeUTR",
                    grepl(feats, pattern = "Exon") ~ "Exon",
                    grepl(feats, pattern = "Intron") ~ "Intron"
                )
            ) %>%
            dplyr::select(name = name.x, feature)

        # Add the intergenic
        olsum_final <- rlbed %>%
            filter(!name %in% olsum$name) %>%
            mutate(
                feature = "Intergenic"
            ) %>%
            dplyr::select(
                name, feature
            ) %>%
            bind_rows(olsum)

        olsum_final %>%
            group_by(feature) %>%
            tally()
    },
    mc.cores = 44
)
names(reslst) <- rlsampsf$rlsample

# Remove nulls
reslstf <- reslst[sapply(reslst, function(x) !is.null(x))]

# Create tbl representation
rlsampleTxOl <- reslstf %>%
    bind_rows(.id = "rlsample") %>%
    group_by(rlsample) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup()

# Save
usethis::use_data(rlsampleTxOl, overwrite = TRUE)
