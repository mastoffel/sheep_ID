# ROH calling 
library(tidyverse)
library(data.table)
library(snpStats)

# # get all ind ids
# ids <- fread("../sheep_imputation/data/hdld_geno_merged.txt", select = 1)
# ids %>% sample_n(1000) %>% mutate(ID2 = ID) %>% fwrite(file = "data/subset_inds.txt", sep = "\t", col.names = FALSE)
# 
# # test ROH on subset
# system("~/programs/plink --bfile data/sheep_imp --keep data/subset_inds.txt --make-bed --out data/geno_sub --sheep")

# LD pruning ===================================================================
system(paste0("~/programs/plink --bfile ../sheep/data/SNP_chip/ramb_mapping/sheep_geno_imputed_ram_27092019 --sheep --out output/ROH/sheep_geno_imputed_ram_27092019_pruned ",
              "--indep-pairwise 500 50 0.95"))

system(paste0("~/programs/plink --bfile ../sheep/data/SNP_chip/ramb_mapping/sheep_geno_imputed_ram_27092019 --sheep ",
              "--extract output/ROH/sheep_geno_imputed_ram_27092019_pruned.prune.in --make-bed --out output/ROH/sheep_geno_imputed_ram_27092019_pruned"))

# PCA for GWAS =================================================================
system(paste0("~/programs/plink --bfile output/ROH/sheep_geno_imputed_ram_27092019_pruned --sheep ",
              "--pca --out output/sheep_pca")) # sheep_pca

# calculate ROH ================================================================
system(paste0("~/programs/plink --bfile ../sheep/data/SNP_chip/ramb_mapping/sheep_geno_imputed_ram_27092019_pruned --sheep --out output/ROH/roh_nofilt_ram_pruned ",
              "--homozyg --homozyg-window-snp 30 --homozyg-snp 25 --homozyg-kb 600 ",
              "--homozyg-gap 500 --homozyg-density 50 --homozyg-window-missing 2 ",
              "--homozyg-het 1 ",
              "--homozyg-window-het 1"))

# inferred ROH output ==========================================================
# without pruning
# file_path <- "output/ROH/roh_nofilt_ram.hom"
# with pruning
file_path <- "output/ROH/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path)

# distribution
hist(roh_lengths$KB, breaks = 1000, xlim = c(500,5000))

# some transformation
froh <- roh_lengths %>%
        dplyr::group_by(IID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/2869914)

# roh stats
roh_lengths %>% 
        group_by(IID) %>% 
        tally() -> roh_nums
range(roh_nums$n)
roh_lengths[which.max(roh_lengths$KB), ]



#~~~~~~~~~~~~~~ calculate homozygosity in the rest of the genome ~~~~~~~~~~~~~~# 
#~~~~~~~~~~~~~~ not used at the moment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

plink_geno_path <- "data/"
# plink name
sheep_plink_name <- "sheep_geno_imputed_ram_27092019_pruned"
sheep_plink_name <- "sheep_geno_imputed_ram_27092019"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, sheep_plink_name, ".bed")
sheep_bim <- paste0(plink_geno_path, sheep_plink_name, ".bim")
sheep_fam <- paste0(plink_geno_path, sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
summary(full_sample$genotypes)
na_genos <- is.na(full_sample$genotypes)

# which individuals have lots of missing data?
missing_per_ind <- rowSums(na_genos)
inds_with_high_miss <- which(missing_per_ind > 0.05 * ncol(full_sample$genotypes))
df <- data.frame(id = names(inds_with_high_miss), missing = missing_per_ind[inds_with_high_miss])
write_delim(df, "data/ids_more_than_5perc_missing_imputation.txt", delim = " ")

# get df
sheep_geno <- as(full_sample$genotypes, Class = "numeric")

# get roh
file_path <- "output/ROH/roh_nofilt_ram_pruned.hom"
file <- "roh_nofilt_ram_pruned"
roh_lengths <- setDF(fread(file_path))
roh_hom_sum <- fread("output/ROH/roh_nofilt_ram_pruned.hom.indiv")

snp_names <- colnames(sheep_geno)
snp_index <- 1:length(snp_names)

snps_df <- data.frame(SNP1 = snp_names, SNP1_index = snp_index, 
                       stringsAsFactors = FALSE)

#"oar3_OAR2_130321792" %in% snp_names
# find SNPs in ROH for each individual
roh <- roh_lengths %>% 
        left_join(snps_df, by = "SNP1") %>% 
        mutate(SNP2_index = SNP1_index + NSNP - 1) %>% 
        mutate(snp_indices = paste0(SNP1_index, ":", SNP2_index)) %>% 
        as_tibble() %>% 
        mutate(snp_indices_full = purrr::map(snp_indices, function(x) eval(parse(text = x)))) %>% 
        group_by(IID) %>% 
        summarise(snps_in_roh = list(unlist(c(snp_indices_full))))

double_check <- unlist(apply(roh,1, function(x) length(x[2][[1]])))
hist(double_check)

# calculate proportion of homozygous SNPs NOT in ROH
calc_hom <- function(id, sheep_geno, roh) {
        id_index <- which(rownames(sheep_geno) == id)
        roh_snps <- roh[roh$IID == id, ]$snps_in_roh[[1]]
        genos <- sheep_geno[id_index, ]
        genos[roh_snps] <- 1 # give it num for heterozygote so they are included in length of vector calculation 
        sum(genos!=1, na.rm = TRUE)/sum(!is.na(genos))
}

homs <- map(roh$IID,calc_hom, sheep_geno, roh)
df <- data.frame(ID = roh$IID, hom = unlist(homs))

write_delim(df, "output/ROH/roh_nofilt_hom_not_in_roh_ram_pruned.txt", delim = " ")
hist(unlist(homs))

# check geno NAs
ind_missing_snps <- rowSums(is.na(sheep_geno))
hist(ind_missing_snps, breaks = 1000)
which(ind_missing_snps > (0.05*400637))
sum(ind_missing_snps > (0.05*400637)) # 392 individuals
