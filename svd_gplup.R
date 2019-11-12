# # run on server
library(lme4)
library(tidyverse)
library(broom.mixed)
library(tidyimpute)
library(matrixStats)
#source("theme_clean.R")
library(snpStats)
library(data.table)

chr <- 20
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
        #as_tibble() %>% 
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

#annual_survival_gwas_raw <- fread("data/annual_survival_gwas_vars.txt")
# join additive and roh data to survival for gwas
annual_survival_gwas <- annual_survival %>% 
        #sample_frac(1) %>% 
        mutate(id = as.character(id)) %>% 
        dplyr::select(id, survival, sex, twin, birth_year, sheep_year, mum_id, age_std, age2_std) %>% 
        filter(!is.na(sex)) %>% 
        left_join(geno_sub, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        as_tibble()


roh_snps <- grep("roh", names(annual_survival_gwas), value = TRUE)
snps <- str_remove(roh_snps, "roh_")

library(RSpectra)
X_roh <- annual_survival_gwas %>% 
                .[roh_snps] %>% 
                impute_mean() %>% 
                as.matrix() %>% 
                scale()

X_add <- annual_survival_gwas %>% 
                .[snps] %>% 
                impute_mean() %>% 
                as.matrix() %>% 
                scale()

# svd_roh <- svds(X_roh, k = 1000)
# svd_add <- svds(X_add, k = 1000)

svd_roh <- svd(X_roh)
saveRDS(svd_roh, "output/svd_roh_chr20.rds")
svd_add <- svd(X_add)
saveRDS(svd_add, "output/svd_add_chr20.rds")

roh_US <- svd_roh$u %*% diag(svd_roh$d)
saveRDS(roh_US, "output/svd_roh_chr20_US.rds")
add_US <- svd_add$u %*% diag(svd_add$d)
saveRDS(add_US, "output/svd_add_chr20_US.rds")


ETA<-list( fixed = list(~factor(sex)+factor(twin)+age_std+age2_std,
                        data=annual_survival_gwas,model='FIXED'),
           #random = list(~factor(id) + factor(sheep_year) + factor(birth_year), data=annual_survival_gwas, model='BRR'),
           id = list(~factor(id), data=annual_survival_gwas, model='BRR'),
           sheep_year = list(~factor(sheep_year), data=annual_survival_gwas, model='BRR'),
           birth_year = list(~factor(birth_year), data=annual_survival_gwas, model='BRR'),
           roh = list(X_roh=roh_US, model='BayesC'), # , saveEffects=TRUE # roh_US
           add = list(X_add=add_US, model = 'BayesC') # , saveEffects=TRUE / BayesC # add_US
)

y <- annual_survival_gwas$survival

library(BGLR)
#3# Fitting the model / previous 10K / 5k burning
fm <- BGLR(y=y,ETA=ETA, nIter=100000, burnIn=20000, thin = 80, response_type = "ordinal",
           saveAt = "output/bglr_chr/chr20_markers")

write_rds(fm, path = "output/bglr_chr/chr20_fm.rds")

plot(fm$ETA$roh$b^2)
# b = Vs
marker_effs_roh <- svd_roh$v %*% fm$ETA$roh$b
marker_effs_add <- svd_add$v %*% fm$ETA$add$b
plot(marker_effs_add^2)

marker_df <- tibble(marker_effs_roh, Index = 1:length(marker_effs_roh) )
marker_df_add <- tibble(marker_effs_add, Index = 1:length(marker_effs_roh) )
source("theme_clean.R")
p1 <- ggplot(marker_df, aes(Index, marker_effs_roh^2)) + geom_point(size = 2, alpha = 0.7) +
        theme_clean() +
        ylab("marker effects ^ 2") +
        xlab("") +
        theme(axis.text.x = element_blank())
ggsave("figs/marker_effs_gwas_by_gblup.jpg", p1,  height = 3, width = 5)

p1 <- ggplot(marker_df_add, aes(Index, marker_effs_add^2)) + geom_point(size = 2, alpha = 0.7) +
        theme_clean() +
        ylab("marker effects ^ 2") +
        xlab("") +
        theme(axis.text.x = element_blank())
ggsave("figs/marker_effs_add_gwas_by_gblup.jpg", p1,  height = 3, width = 5)

gwas_roh %>% 
        filter(chromosome == 20) %>% 
        mutate(Index = 1:nrow(.)) %>% 
        ggplot(aes(Index, -log10(p.value))) + geom_point(size = 2, alpha = 0.7) +
        theme_clean() +
        ylab("marker effects pval") +
        xlab("") +
        theme(axis.text.x = element_blank())

# code from another place
gwas_roh <- gwas_full %>% 
        filter(groups == "roh") 

p3 <- gwas_roh %>% 
        filter(chromosome == 20) %>% 
        mutate(Index = 1:nrow(.)) %>% 
        ggplot(aes(Index, estimate^2)) + geom_point(size = 2, alpha = 0.7) +
        theme_clean() +
        ylab("marker effects ^ 2") +
        xlab("") +
        theme(axis.text.x = element_blank())
ggsave("figs/marker_effs_roh_gwas_single.jpg", p3,  height = 3, width = 5)


gwas_add <- gwas_full %>% 
        filter(groups == "add") 

p4 <- gwas_add %>% 
        filter(chromosome == 20) %>% 
        mutate(Index = 1:nrow(.)) %>% 
        ggplot(aes(Index, estimate^2)) + geom_point(size = 2, alpha = 0.7) +
        theme_clean() +
        ylab("marker effects ^ 2") +
        xlab("") +
        theme(axis.text.x = element_blank())
p4
ggsave("figs/marker_effs_add_gwas_single.jpg", p4,  height = 3, width = 5)


out <- read_delim("output/bglr_chr/chr20_markersETA_roh_parBayesC.dat", delim = " ", col_names = F)
plot(out$X2)

fwrite(X_roh, "data/roh_mat_chr20.txt", col.names = FALSE, row.names = FALSE)
fwrite(X_add, "data/add_mat_chr20.txt", col.names = FALSE, row.names = FALSE)


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
