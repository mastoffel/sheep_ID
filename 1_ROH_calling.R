# (1) filtering the imputed data set for low call rate individuals
# (2) run plink pca on genotypes to get PCs for GWAS
# (3) calling ROH 

# relatedness PCs for GWAS
library(tidyverse)
library(data.table)
library(snpStats)

# imputation based on Oar 3.1 mapping
# retain individuals with survival data
system(paste0("/usr/local/bin/plink --bfile ../sheep/data/SNP_chip/oar31_mapping/sheep_geno_imputed_oar31_17052020 ",
              "--sheep --keep output/ROH/ids_surv.txt ",
              "--make-bed --out data/sheep_geno_imputed_oar_indfilt"))

# check SNP call rates ---------------------------------------------------------
# plink name
sheep_plink_name <- "data/sheep_geno_imputed_oar_indfilt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# snp genotyping summary
snp_sum <- col.summary(full_sample$genotypes) %>% 
                as_tibble(rownames = "snp.name") %>% 
                left_join(full_sample$map, ., by = "snp.name") 
# filter call rate and monomorphic snps
snp_sum %>% 
        as_tibble() %>% 
        filter((Call.rate < 0.95) | (MAF == 0)) %>% 
        .$snp.name %>% 
        write_lines("data/oar_imp_low_call_snps095.txt")

# filter low call rate SNPs
system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar_indfilt --sheep --exclude data/oar_imp_low_call_snps095.txt ",
              "--make-bed --out data/sheep_geno_imputed_oar_filt"))

# check genotyping rates -------------------------------------------------------
# plink name
sheep_plink_name <- "data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
# individual genotyping summary
ind_sum <- row.summary(full_sample$genotypes)
summary(ind_sum)
ind_sum[ind_sum$Call.rate < 0.95, ] # two individuals slightly under 95%, that's alright.
# snp summary
snp_sum <- col.summary(full_sample$genotypes)
summary(snp_sum) 

# PCA for GWAS =================================================================

# for 400k PCA we need to prune data as unpruned data fails for some reason
system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar_filt --sheep ",
              "--indep-pairwise 500 50 0.999 ",
              "--out data/sheep_geno_imputed_oar_filt_pruned"))

system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar_filt --sheep ",
              "--pca --exclude data/sheep_geno_imputed_oar_filt_pruned.prune.out ",
              "--out output/sheep_pca_oar")) # sheep_pca

# calculate ROH 
system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar_filt --sheep --out output/ROH/roh ",
              # "--keep output/ROH/ids_surv.txt ",
              "--homozyg --homozyg-window-snp 50 --homozyg-snp 50 --homozyg-kb 1200 ",
              "--homozyg-gap 300 --homozyg-density 200 --homozyg-window-missing 2 ",
              "--homozyg-het 2 ",
              "--homozyg-window-het 2"))



















# unused content
#~~~~~~~~~~~~~~ calculate homozygosity in the rest of the genome ~~~~~~~~~~~~~~# 
#~~~~~~~~~~~~~~ not used at the moment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

plink_geno_path <- "data/"
# plink name
sheep_plink_name <- "sheep_geno_imputed_ram_398k_filt"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, sheep_plink_name, ".bed")
sheep_bim <- paste0(plink_geno_path, sheep_plink_name, ".bim")
sheep_fam <- paste0(plink_geno_path, sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

summary(full_sample$genotypes)
qc_geno <- row.summary(full_sample$genotypes)
hist(qc_geno$Call.rate)

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


# genotype file for simpleM

# which ids are HD
plink_geno_path <- "../sheep/data/SNP_chip/oar31_mapping/Plates_1-2_HD_QC2"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, ".bed")
sheep_bim <- paste0(plink_geno_path, ".bim")
sheep_fam <- paste0(plink_geno_path, ".fam")
hd_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
write_delim(data.frame(id = hd_sample$fam$member), path = "data/hd_ids.txt", " ")

id <- hd_sample$fam$member
# plink name
#sheep_plink_name <- "sheep_geno_imputed_ram_27092019_pruned"
sheep_plink_name <- "data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

hd_inds <- which(full_sample$fam$member %in% id)

genos_hd <- as(full_sample$genotypes[hd_inds, ], Class = "numeric")
# in the few cases with NA add random genotype so that simpleM works
library(imputeTS)
genos_hd <- na_mean(genos_hd)
#genos_hd[is.na(genos_hd)] <- sample(c(0,1,2), replace=TRUE, sum(is.na(genos_hd)))
fwrite(data.table::transpose(as.data.table(genos_hd)), "data/geno_mat_simpleM_allchr_oar.txt")
