# ROH calling 
library(tidyverse)
library(data.table)


# # get all ind ids
# ids <- fread("../sheep_imputation/data/hdld_geno_merged.txt", select = 1)
# ids %>% sample_n(1000) %>% mutate(ID2 = ID) %>% fwrite(file = "data/subset_inds.txt", sep = "\t", col.names = FALSE)
# 
# # test ROH on subset
# system("~/programs/plink --bfile data/sheep_imp --keep data/subset_inds.txt --make-bed --out data/geno_sub --sheep")

# calculate ROH according to Joshi et al. (unpruned)
system(paste0("~/programs/plink --bfile ../sheep/data/SNP_chip/sheep_geno_imputed_10062019 --sheep --out output/ROH/roh_nofilt ",
              "--homozyg --homozyg-window-snp 50 --homozyg-snp 50 --homozyg-kb 1000 ",
              "--homozyg-gap 1000 --homozyg-density 50 --homozyg-window-missing 5 ",
              "--homozyg-window-het 1"))

file_path <- "output/ROH/roh_nofilt.hom"
file <- "roh_nofilt"
roh_lengths <- fread(file_path)

froh <- roh_lengths %>%
        dplyr::group_by(FID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/2869490)

hist(froh$FROH, breaks = 100)
# Specifiy length distributon:


# with pruning
system(paste0("~/programs/plink --bfile ../sheep/data/SNP_chip/sheep_geno_imputed_10062019 --sheep --out output/ROH/sheep_geno_imputed_LDpruned ",
              "--indep-pairwise 1000 5 0.7"))
system(paste0("~/programs/plink --bfile ../sheep/data/SNP_chip/sheep_geno_imputed_10062019 --sheep --out output/ROH/roh_nofilt_pruned ",
              "--extract output/ROH/sheep_geno_imputed_LDpruned.prune.in ",
              "--homozyg --homozyg-window-snp 50 --homozyg-snp 50 --homozyg-kb 1000 ",
              "--homozyg-gap 1000 --homozyg-density 250 --homozyg-window-missing 5 ",
              "--homozyg-window-het 1"))

file_path <- "output/ROH/roh_nofilt_pruned.hom"
file <- "roh_nofilt_pruned"
roh_lengths_pruned <- fread(file_path)

froh_pruned <- roh_lengths_pruned %>%
        dplyr::group_by(FID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/2869490)

hist(froh_pruned$FROH, breaks = 100)

plot(froh$FROH, froh_pruned$FROH)
cor(froh$FROH, froh_pruned$FROH)

