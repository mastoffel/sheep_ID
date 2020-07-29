# ROH GWAS. This script runs a modified GWAS which tests an
# effect of ROH status on annual survival at every SNP, controlling
# for a range of other variables.
# Needs to be run on a cluster, every model takes appr. 40 sec to run,
# 417K models overall.

library(lme4)
library(tidyverse)
library(broom.mixed)
library(snpStats)
library(data.table)
library(furrr)

# fitness and pedigree data
load("data/survival_mods_data.RData")
load("data/sheep_ped.RData")

# GRM PCs from plink
pcs <- read_delim("data/ann_surv_pca.txt", " ", col_names = TRUE) %>% 
        mutate(id = as.character(id))

# roh data
file_path <- "output/ROH/roh.hom"
roh_lengths <- fread(file_path)

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map <- full_sample$map


# GWAS results
gwas_res <- read_rds("output/gwas_res_oar_long.rds")

# put into df
gwas_full <- gwas_res %>%
        rename(snp.name = term) %>%
        filter(!str_detect(snp.name, "sd")) %>% 
        mutate(groups = ifelse(str_detect(snp.name, "roh"), "roh", "add")) %>% 
        mutate(snp.name = str_replace(snp.name, "roh_", "")) %>%
        left_join(snps_map) 

# extract roh
gwas_roh <- gwas_full %>% filter(groups == "roh") 
top_snps <- gwas_roh %>% 
        group_by(chromosome) %>% 
        filter(p.value < 0.000005) %>% 
        arrange(p.value)

# index of top snps
snp_indices <- which(colnames(full_sample$genotypes) %in% top_snps$snp.name)

# filter map data
snps_map_sub <- as_tibble(full_sample$map[snp_indices, ])

# additive genotypes
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_map_sub$snp.name], Class = "numeric"),
                      rownames = "id")
# survival data
# survival data preprocessing
annual_survival <- fitness_data %>% 
        # filter na rows
        filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
        mutate(age_cent = age - mean(age, na.rm = TRUE),
               age_cent2 = age_cent^2,
               age_std = as.numeric(scale(age)),
               age_std2 = age_std^2,
               # times 10 to estimate a 10% percent increase
               froh_all10 = froh_all * 10,
               froh_all10_cent = froh_all10 - mean(froh_all10, na.rm = TRUE),
               lamb = ifelse(age == 0, 1, 0),
               lamb_cent = lamb - mean(lamb, na.rm = TRUE),
               lamb = as.factor(lamb)) %>% 
        as.data.frame() 

#roh_lengths <- as.data.table(roh_lengths)
# check whether snp is in ROH for a given individual

setkey(roh_lengths, IID)
roh_id_per_snp <- function(i) {
        position <- as.numeric(snps_map_sub[i, "position"])
        chromosome <- as.numeric(snps_map_sub[i, "chromosome"])
        # varname <- paste0("roh", i)
        #roh <- as.numeric((roh_lengths$POS1 <= position) & (roh_lengths$POS2 >= position) & (roh_lengths$CHR == chromosome))
        #roh_lengths$roh <- roh
        roh_lengths[, roh := as.numeric((CHR == chromosome) & (POS1 <= position) & (POS2 >= position))]
        #roh_lengths[, roh := fifelse((POS1 <= position)&(POS2 >= position)&(CHR == chromosome), 1, 0)]
        roh_id <- roh_lengths[,  .(roh = max(roh)), by = c("IID")]$roh
}

roh_ind <- map(1:nrow(snps_map_sub), roh_id_per_snp)
roh_df <- as.data.frame(do.call(cbind, roh_ind))
names(roh_df) <- paste0("roh_", snps_map_sub$snp.name)
roh_df$id <- as.character(unique(roh_lengths$IID))


# make some space
rm(full_sample)

# which chromosomes do the snps span?
chrs <- unique(snps_map_sub$chromosome)
froh_no_chr <- paste0("froh_no_chr", chrs)
# join additive and roh data to survival for gwas
annual_survival_gwas <- annual_survival %>% 
        #mutate_at(vars(starts_with("froh_no_chr")), scale) %>% 
        dplyr::select(id, survival, sex, twin, lamb, birth_year, sheep_year, mum_id, age_std, age_std2, {{ froh_no_chr }}) %>% 
        left_join(pcs, by = "id") %>% 
        left_join(geno_sub, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        as_tibble()

#test <- annual_survival_gwas[c(1, 22,23, 25, 26)]
#snp_names <- names(test)[c(2,3)]

# make major and minor allele roh
snp_names <- snps_map_sub$snp.name
#snp_names <- top_snps %>% group_by(chromosome) %>% top_n(-2, p.value) %>% .$snp.name
for (i in snp_names) {
        annual_survival_gwas[[paste0("roh_", i, "_0")]] <- as.numeric((annual_survival_gwas[[i]] == 0) & (annual_survival_gwas[[paste0("roh_", i)]] == 1))
        annual_survival_gwas[[paste0("roh_", i, "_2")]] <- as.numeric((annual_survival_gwas[[i]] == 2) & (annual_survival_gwas[[paste0("roh_", i)]] == 1))
        snp_fac <- paste0("roh_fac_", i)
        # make a factor
        annual_survival_gwas[snp_fac] <- case_when(
                (annual_survival_gwas[[paste0("roh_", i)]] == 0) ~ 0, # (annual_survival_gwas[[i]] == 1) & 
                (annual_survival_gwas[[i]] == 0) & (annual_survival_gwas[[paste0("roh_", i)]] == 1) ~ 1,
                (annual_survival_gwas[[i]] == 2) & (annual_survival_gwas[[paste0("roh_", i)]] == 1) ~ 2
        ) %>% 
                as.factor()
}


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

# focal SNP, chromosome of focal snp, data
run_gwas <- function(snp, data) {
        # for mean froh without focal chr
        chr <- as.numeric(snps_map_sub[snps_map_sub$snp.name == snp, "chromosome"])
        froh_no_chr <- paste0("froh_no_chr", chr)
        
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + age_std + age_std2 + ", 
                                         froh_no_chr, " + ",
                                         "pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + ",
                                         #"pc1 + pc2 + pc3 + pc4 +",
                                         snp, "+ ", paste0("roh_fac_", snp), "+ (1|birth_year) + (1|sheep_year) + (1|id)"))
                                        # snp, "+ ", paste0("roh_", snp, "_0"), "+ ", paste0("roh_", snp, "_2"), "+ (1|birth_year) + (1|sheep_year) + (1|id)"))
        #snp, "+ ", paste0("roh_", snp), " + (1|sheep_year) + (1|id)"))
        mod <- glmer(formula = formula_snp,
                     data = data, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod)
        out
}


out <- map(snp_names, run_gwas, annual_survival_gwas)
out %>% 
        bind_rows(.id = "snp") %>% 
        filter(str_detect(term, "^roh")) %>% 
        arrange(p.value)


#1 oar3_OAR10_85579446   -1.08  0.0000000355 roh            10 105.   85579446 G        A       
#2 s23340.1              -0.822 0.000000176  roh            23  NA    36164028 G        A       
#3 DU289160_204.1        -0.739 0.000000771  roh 


safe_run_gwas <- purrr::safely(run_gwas)

# 

snps_sub <- snps_map_sub$snp.name
# split into pieces of 50 SNPs 
num_parts <- round(length(seq_along(snps_sub )) / 50)
snps_pieces <- split(snps_sub, cut(seq_along(snps_sub), num_parts, labels = FALSE))
roh_pieces <- map(snps_pieces, function(x) paste0("roh_", x))

annual_survival_gwas_pieces <- 
        map2(snps_pieces, roh_pieces, function(snps_piece, roh_piece) {
                annual_survival_gwas %>% dplyr::select(id:pc7, one_of(c(snps_piece, roh_piece)))   
        })

# clean up
rm(annual_survival, annual_survival_gwas, fitness_data, geno_sub, 
   roh_lengths, roh_pieces, sheep_ped, roh_df)

# set up plan
plan(multiprocess, workers = 3)

# increase maxSize
options(future.globals.maxSize = 3000 * 1024^2)

all_out <- future_map2(snps_pieces, annual_survival_gwas_pieces, function(snps, data) {
        out <- purrr::map(snps, safe_run_gwas, data)
})

all_out_simple <- purrr::flatten(all_out)
saveRDS(all_out_simple, file = paste0("output/GWAS_roh_part", "_", part, ".rds"))








