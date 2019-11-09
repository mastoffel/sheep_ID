# # run on server
library(lme4)
library(tidyverse)
library(broom.mixed)
#source("theme_clean.R")
library(snpStats)
library(data.table)
library(furrr)
library(caret)
# for running on server
chr_inp  <- commandArgs(trailingOnly=TRUE)
if (!(length(chr_inp) == 0)) {
        chr <- as.numeric(chr_inp[[1]])
} else {
        # SNP data
        # which chromosome
        chr <- 20
}


# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# pcas 
pcs <- read_delim("data/ann_surv_pca.txt", " ", col_names = TRUE) %>% 
        mutate(id = as.character(id))

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

snps_map_sub <- full_sample$map %>% 
        filter(chromosome == chr) 

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


# #use only training set (80% of individuals)
# sample_frac_groups = function(tbl, size, replace = FALSE, weight = NULL) {
#         # regroup when done
#         grps = tbl %>% groups %>% lapply(as.character) %>% unlist
#         # check length of groups non-zero
#         keep = tbl %>% summarise() %>% ungroup() %>% sample_frac(size, replace, weight)
#         # keep only selected groups, regroup because joins change count.
#         # regrouping may be unnecessary but joins do something funky to grouping variable
#         tbl %>% right_join(keep, by=grps) %>% group_by(.dots = grps)
# }
# 
# # set.seed(7336)
# early_survival <- early_survival %>%
#         mutate(index = 1:nrow(.)) %>%
#         group_by(id) %>%
#         sample_frac_groups(0.8) %>%
#         ungroup()
#write_lines(early_survival$index, "output/ind_index_test.txt")
#write_lines(unique(early_survival$id), "data/ind_testset_80.txt")

# prepare additive genotypes subset
snps_sub <- full_sample$map %>% 
        filter(chromosome == chr) %>% 
        .$snp.name
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_sub], Class = "numeric"),
                      rownames = "id")

# subset roh 
roh_sub <- roh_lengths %>% filter(CHR == chr) %>% filter(KB > 1000)

# define vectorized seq to work with mutate
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

# create indices for all rohs
roh_snps <- roh_sub %>% 
        as_tibble() %>% 
        # sample_frac(0.01) %>% 
        mutate(index1 = as.numeric(match(SNP1, names(geno_sub))),
               index2 = as.numeric(index1 + NSNP - 1)) %>% 
        mutate(all_snps = seq2(from = index1, to = index2)) %>% 
        group_by(IID) %>% 
        summarise(all_snps = list(all_snps)) %>% 
        mutate(all_snps = simplify_all(all_snps)) %>% 
        mutate(IID = as.character(IID)) %>% 
        rename(id = IID)

# join roh_snps
roh_snps_reord <- geno_sub %>% 
        dplyr::select(id) %>% 
        left_join(roh_snps, by = "id") %>% 
        dplyr::rename(id_roh = id)

# prepare roh yes/no matrix
roh_mat <- matrix(data = 0, nrow = nrow(geno_sub), ncol = ncol(geno_sub))

# set 1 where SNP is in an roh
roh_list <- pmap(roh_snps_reord, function(id_roh, all_snps) {
        df <- as.matrix(t(as.numeric(c(id_roh, rep(0, ncol(geno_sub) - 1)))))
        df[, all_snps] <- 1
        df
}) 

# make a tibble like for genotypes but with 0/1 for whether a SNP is in ROH or not
roh_df <- do.call(rbind, roh_list) %>% 
        as_tibble() %>% 
        mutate(V1 = as.character(V1))
names(roh_df) <- c("id", paste0("roh_", names(geno_sub)[-1]))


# make some space
rm(full_sample)
rm(roh_list)
rm(roh_mat)

annual_survival_gwas_raw <- fread("data/annual_survival_gwas_vars.txt")
# join additive and roh data to survival for gwas
annual_survival_gwas <- annual_survival_gwas_raw %>% 
        sample_frac(1) %>% 
        mutate(id = as.character(id)) %>% 
        dplyr::select(id, survival, sex, twin, birth_year, sheep_year, mum_id, age_std, age2_std) %>% 
        left_join(pcs, by = "id") %>% 
        left_join(geno_sub, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        as_tibble()

roh_snps <- grep("roh", names(annual_survival_gwas))
snps <- annual_survival_gwas[17:3354]

X_roh <- annual_survival_gwas %>% 
                .[roh_snps] %>% 
                impute_mean() %>% 
                as.matrix()

X_add <- snps %>% 
                impute_mean() %>% 
                as.matrix()


svd_roh <- svd(X_roh)
svd_add <- svd(X_add)

roh_US <- svd_roh$u %*% diag(svd_roh$d)
add_US <- svd_add$u %*% diag(svd_add$d)

ETA<-list( fixed = list(~factor(sex)+factor(twin)+age_std+age2_std,
                        data=annual_survival_gwas,model='FIXED'),
           #random = list(~factor(id) + factor(sheep_year) + factor(birth_year), data=annual_survival_gwas, model='BRR'),
           id = list(~factor(id), data=annual_survival_gwas, model='BRR'),
           sheep_year = list(~factor(sheep_year), data=annual_survival_gwas, model='BRR'),
           birth_year = list(~factor(birth_year), data=annual_survival_gwas, model='BRR'),
           roh = list(X_roh=roh_US, model='BayesC'), # , saveEffects=TRUE
           add = list(X_add=add_US, model = 'BayesC') # , saveEffects=TRUE / BayesC
)

y <- annual_survival_gwas$survival

library(BGLR)
#3# Fitting the model / previous 10K / 5k burning
fm <- BGLR(y=y,ETA=ETA, nIter=10000, burnIn=1000, thin = 10, response_type = "ordinal")

# b = Vs
marker_effs <- svd_roh$v %*% fm$ETA$roh$b
plot(marker_effs)




# blup of marker effects
sigma2a <- 0.1
sigma2e <- 0.5
lambda <- sigma2e/sigma2a
y <- as.numeric(as.character(annual_survival_gwas$survival))
W <- X_roh
Wt <- t(W)
WWt <- W%*%Wt
Ivar <- lambda*(diag(nrow(W)))
WWt2 <- WWt + Ivar
WWt2_inv <- solve(WWt2)
blup_a <- (Wt %*% WWt2_inv) %*% y
plot(blup_a ^ 2)

#blup_a <- (Wt %*% solve(W%*%Wt + lambda*(diag(nrow(W))))) %*% y
