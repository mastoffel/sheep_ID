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
pcs <- read_delim("data/ann_surv_pca_testset80.txt", " ", col_names = TRUE) %>% 
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


#use only training set (80% of individuals)
annual_survival <- annual_survival %>% filter(id %in% pcs$id)


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

# join additive and roh data to survival for gwas
annual_survival_gwas <- annual_survival %>% 
        dplyr::select(id, survival, sex, twin, birth_year, sheep_year, mum_id, age_std, age2_std) %>% 
        left_join(pcs, by = "id") %>% 
        left_join(geno_sub, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        as_tibble()

# time saver function for modeling
nlopt <- function(par, fn, lower, upper, control) {
        .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
                                  opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                              maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
        list(par = res$solution,
             fval = res$objective,
             conv = if (res$status > 0) 0 else res$status,
             message = res$message
        )
}

# GWAS
# run_gwas <- function(snp, data) {
#         formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + ", snp, "+ ", paste0("roh_", snp), "+ (1|birth_year) + (1|sheep_year) + (1|id)"))
#         mod <- glmer(formula = formula_snp,
#                      data = data, family = "binomial",
#                      control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
#         out <- broom.mixed::tidy(mod)
#         out
# }

run_gwas <- function(snp, data) {
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + age_std + age2_std + ", 
                                         "pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + ",
                                         snp, "+ ", paste0("roh_", snp), "+ (1|birth_year) + (1|sheep_year) + (1|id)"))
        mod <- glmer(formula = formula_snp,
                     data = data, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod)
        out
}

safe_run_gwas <- purrr::safely(run_gwas)
# gwas_out <- purrr::map(snps_sub, safe_run_gwas, early_survival_gwas)

# split in pieces of 1000 snps / rohs, each approximately 0.25 Gb)
num_parts <- round(length(seq_along(snps_sub)) / 1000)
snps_pieces <- split(snps_sub, cut(seq_along(snps_sub), num_parts, labels = FALSE))
roh_pieces <- map(snps_pieces, function(x) paste0("roh_", x))

annual_survival_gwas_pieces <- 
        map2(snps_pieces, roh_pieces, function(snps_piece, roh_piece) {
                annual_survival_gwas %>% dplyr::select(id:pc7, one_of(c(snps_piece, roh_piece)))   
        })

# clean up
rm(annual_survival, annual_survival_gwas, fitness_data, geno_sub, roh_lengths, roh_pieces, 
   roh_snps, roh_snps_reord, sheep_ped, snps_map_sub, roh_sub, roh_df)

# set up plan
plan(multiprocess, workers = 8)

# increase maxSize
options(future.globals.maxSize = 3000 * 1024^2)
all_out <- purrr::pmap(list(snps_pieces, annual_survival_gwas_pieces, 1:num_parts),  function(snps, data, num_part) {
        out <- future_map(snps, safe_run_gwas, data)
        # remove one hierarchical level
        all_out_simple <- purrr::flatten(out)
        saveRDS(all_out_simple, file = paste0("output/GWAS_roh_chr", chr, "_", num_part, ".rds"))
})









