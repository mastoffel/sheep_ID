# # run on server
library(lme4)
library(tidyverse)
library(broom.mixed)
#source("theme_clean.R")
library(snpStats)
library(data.table)
library(furrr)
library(caret)
library(tidyimpute)
library(BGLR)

# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# roh data
file_path <- "data/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path) 

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019_pruned"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

snps_map_sub <- full_sample$map %>% write_delim("data/snp_map.txt")
        #sample_n(2000)

# survival data
annual_survival <- fitness_data %>% 
        dplyr::rename(birth_year = BIRTHYEAR,
                      sheep_year = SheepYear,
                      age = Age,
                      id = ID,
                      twin = TWIN,
                      sex = SEX,
                      mum_id = MOTHER,
                      froh_short = FROH_short,
                      froh_medium = FROH_medium,
                      froh_long = FROH_long,
                      froh_all = FROH_all,
                      froh_not_roh = hom,
                      survival = Survival) %>% 
        # some individuals arent imputed well and should be discarded 
        filter(!(id %in% IDs_lots_missing$id)) %>% 
        #filter(age == 0) %>% 
        filter(!is.na(survival)) %>% 
        filter(!is.na(froh_all)) %>% 
        filter(!(is.na(birth_year) | is.na(sheep_year))) %>%  # no mum_id here
        mutate_at(c("id", "birth_year", "sex", "sheep_year", "survival"), as.factor) %>% 
        mutate(age2 = age^2) %>% 
        mutate(age_std = as.numeric(scale(age)),
               age2_std = as.numeric(scale(age2))) %>% 
        as.data.frame() 

# prepare additive genotypes subset
snps_sub <- full_sample$map %>% 
        # sample_n(2000) %>% 
        .$snp.name
geno_add <- as.data.table(as(full_sample$genotypes[, snps_sub], Class = "numeric"),
                      keep.rownames = "id")

# subset roh 
roh_sub <- roh_lengths %>% filter(KB > 1000)

# define vectorized seq to work with mutate
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

# create indices for all rohs
roh_snps <- roh_sub %>% 
        as_tibble() %>% 
        # sample_frac(0.01) %>% 
        mutate(index1 = as.numeric(match(SNP1, names(geno_add))),
               index2 = as.numeric(index1 + NSNP - 1)) %>% 
        mutate(all_snps = seq2(from = index1, to = index2)) %>% 
        group_by(IID) %>% 
        summarise(all_snps = list(all_snps)) %>% 
        mutate(all_snps = simplify_all(all_snps)) %>% 
        mutate(IID = as.character(IID)) %>% 
        rename(id = IID)

# join roh_snps
roh_snps_reord <- geno_add %>% 
        dplyr::select(id) %>% 
        left_join(roh_snps, by = "id") %>% 
        dplyr::rename(id_roh = id)

# prepare roh yes/no matrix
roh_mat <- matrix(data = 0, nrow = nrow(geno_add), ncol = ncol(geno_add))

# set 1 where SNP is in an roh
roh_list <- pmap(roh_snps_reord, function(id_roh, all_snps) {
        df <- as.matrix(t(as.numeric(c(id_roh, rep(0, ncol(geno_add) - 1)))))
        df[, all_snps] <- 1
        df
}) 

# make a tibble like for genotypes but with 0/1 for whether a SNP is in ROH or not
roh_df <- do.call(rbind, roh_list) %>% 
        as_tibble() %>% 
        mutate(V1 = as.character(V1))
names(roh_df) <- c("id", paste0("roh_", names(geno_add)[-1]))

# make some space
rm(full_sample)
rm(roh_list)
rm(roh_mat)

# join additive and roh data to survival for gwas
annual_survival_gwas <- annual_survival %>% 
        dplyr::select(id, survival, sex, twin, birth_year, sheep_year, mum_id, age_std, age2_std) %>% 
        #left_join(pcs, by = "id") %>% 
        left_join(geno_add, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        as_tibble() %>% 
        filter(!(is.na(twin)) & !is.na(sex) & !is.na(birth_year) & !is.na(sheep_year) & !is.na(age_std))
        
fwrite(annual_survival_gwas, file = "data/annual_survival_gwas_df.txt")
rm(list=setdiff(ls(), "annual_survival_gwas"))






#~~~~~~~~~~~~ Impute missing ROH with means and save as txt ~~~~~~~~~~~~~~~~~~~#
# try to split up the imputation problem
df_names <- names(fread("grep year data/annual_survival_gwas_df.txt",
                                    nrows = 0))
# indices of roh
roh_ind <- which(str_detect(df_names, "roh"))
sub_roh <- split(roh_ind, cut_number(1:length(roh_ind), 30))

impute_roh <- function(roh_indices) {
        rohs <- fread("data/annual_survival_gwas_df.txt", select = roh_indices)
        out <- rohs %>% impute_mean() %>% as.matrix()
}

# get imputed roh matrix
roh_mat_list <- map(sub_roh, impute_roh)
roh_mat <- do.call(cbind, roh_mat_list)

fwrite(x = roh_mat, "data/annual_survival_gwas_roh.txt")

#~~~~~~~~~~~~ Impute missing Genos with means and save as txt ~~~~~~~~~~~~~~~~~~~#
df_names <- names(fread("grep year data/annual_survival_gwas_df.txt",
                        nrows = 0))
not_roh_ind <- which(!str_detect(df_names, "roh"))
# remove roh snps and non genetic variables to get names of additive snps
# find indices of additive snps
add_snps <- which(df_names %in% df_names[not_roh_ind][-c(1:9)])
sub_snps <- split(add_snps, cut_number(1:length(add_snps), 30))

impute_snps <- function(snp_indices) {
        snps <- fread("data/annual_survival_gwas_df.txt", select = snp_indices)
        out <- snps %>% impute_mean() %>% as.matrix()
}

# get imputed roh matrix
snp_mat_list <- map(sub_snps, impute_snps)
snp_mat <- do.call(cbind, snp_mat_list)
fwrite(x = snp_mat, "data/annual_survival_gwas_snps.txt")

#~~~~~~~~~~~~ GWAS non-genetic variables  ~~~~~~~~~~~~~~~~~~~#
df_names <- names(fread("grep year data/annual_survival_gwas_df.txt",
                        nrows = 0))
GWAS_non_gen <- fread("data/annual_survival_gwas_df.txt", select = 1:9)
fwrite(GWAS_non_gen, "data/annual_survival_gwas_vars.txt")

