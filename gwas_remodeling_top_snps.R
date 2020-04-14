# check models for top SNPs

# which were the top snps before

library(tidyverse)
library(snpStats)
source("theme_simple.R")
library(viridis)
library(magrittr)
library(patchwork)
library(furrr)
library(data.table)

# first, find top snps from old gwas -------------------------------------------
chr_info <- read_delim("../sheep/data/sheep_genome/chromosome_info_ram.txt", "\t") %>% 
        .[-1, ] %>% 
        rename(chromosome = Part) %>% 
        mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
        mutate(chromosome = as.integer(chromosome)) %>% 
        filter(!is.na(chromosome))

gwas_files <- list.files("output/gwas_full_pca_with_f", pattern = "*.rds", full.names = TRUE)

all_gwas <- purrr::map(gwas_files, readRDS) %>% 
        purrr::flatten() %>% 
        .[seq(1,length(.),by=2)] %>% 
        purrr::compact()
# remove snps that didnt work

# set up plan
plan(multiprocess, workers = 8)
gwas_res <- future_map_dfr(all_gwas, function(x) x %>% .[c(15), ] %>% 
                           dplyr::select(term, estimate, p.value))
gwas_res <- gwas_res[-which(str_detect(gwas_res$term, "sd")), ]
top_snps <- gwas_res %>% filter(-log10(p.value) > 5.6) %>% 
                mutate(term = str_replace(term, "roh_", ""))

# now, model the same SNPs in new gwas -----------------------------------------

# data
load("data/survival_mods_data.RData")
load("data/sheep_ped.RData")
#IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# pcs 
pcs <- read_delim("data/ann_surv_pca.txt", " ", col_names = TRUE) %>% 
        mutate(id = as.character(id))

# roh data
file_path <- "data/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path)

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_pruned"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# filter map data
snps_map_sub <- full_sample$map %>% 
                        filter(snp.name %in% top_snps$term) %>% 
                        as_tibble()
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
        roh_lengths[, roh := as.numeric((CHR == chromosome) & (POS1 <= position) & (POS2 >= position))]
        roh_id <- roh_lengths[,  .(roh = max(roh)), by = c("IID")]$roh
}

roh_ind <- map(1:nrow(snps_map_sub), roh_id_per_snp)
roh_df <- as.data.frame(do.call(cbind, roh_ind))
names(roh_df) <- paste0("roh_", snps_map_sub$snp.name)
roh_df$id <- as.character(unique(roh_lengths$IID))

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
        
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + age_std + age_std2 + lamb + ", 
                                         froh_no_chr, " + ",
                                         "pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + ",
                                         #"pc1 + pc2 + pc3 + pc4 +",
                                         snp, "+ ", paste0("roh_", snp), "+ (1|birth_year) + (1|sheep_year) + (1|id)"))
        #snp, "+ ", paste0("roh_", snp), " + (1|sheep_year) + (1|id)"))
        mod <- glmer(formula = formula_snp,
                     data = data, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod)
        out
}

plan(multiprocess, workers = 8)
gwas1 <- future_map(snps_map_sub$snp.name, run_gwas, annual_survival_gwas)
gwas1 %>% bind_rows() %>% filter(str_detect(term, "^roh"))


# focal SNP, chromosome of focal snp, data
run_gwas2 <- function(snp, data) {
        # for mean froh without focal chr
        chr <- as.numeric(snps_map_sub[snps_map_sub$snp.name == snp, "chromosome"])
        froh_no_chr <- paste0("froh_no_chr", chr)
        
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + age_std + age_std2 + ", 
                                         froh_no_chr, " + ",
                                         "pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + ",
                                         #"pc1 + pc2 + pc3 + pc4 +",
                                         snp, "+ ", paste0("roh_", snp), "+ (1|birth_year) + (1|sheep_year) + (1|id)"))
        #snp, "+ ", paste0("roh_", snp), " + (1|sheep_year) + (1|id)"))
        mod <- glmer(formula = formula_snp,
                     data = data, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod)
        out
}

plan(multiprocess, workers = 8)
gwas2 <- future_map(snps_map_sub$snp.name, run_gwas2, annual_survival_gwas)
gwas2 %>% bind_rows() %>% filter(str_detect(term, "^roh")) %>% mutate(log10 = -log10(p.value))

# focal SNP, chromosome of focal snp, data
run_gwas3 <- function(snp, data) {
        # for mean froh without focal chr
        chr <- as.numeric(snps_map_sub[snps_map_sub$snp.name == snp, "chromosome"])
        froh_no_chr <- paste0("froh_no_chr", chr)
        
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + age_std + age_std2 + ", 
                                         froh_no_chr, " + ",
                                         "pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + ",
                                         #"pc1 + pc2 + pc3 + pc4 +",
                                         snp, "+ ", paste0("roh_", snp), "+ (1|birth_year) + (1|sheep_year) + (1|id) "))
        
        formula_snp <- as.formula(paste0("survival ~ 1 + ", 
                                         froh_no_chr, " + ",
                                         "pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + ",
                                         #"pc1 + pc2 + pc3 + pc4 +",
                                         snp, "+ ", paste0("roh_", snp), "+ (1|birth_year) + (1|sheep_year) + (1|id) + sex + twin + age_std + age_std2 "))
        
        #snp, "+ ", paste0("roh_", snp), " + (1|sheep_year) + (1|id)"))
        mod <- glmer(formula = formula_snp,
                     data = data, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod)
        out
}

plan(multiprocess, workers = 8)
gwas3 <- future_map(snps_map_sub$snp.name, run_gwas3, annual_survival_gwas)
gwas3 %>% bind_rows() %>% filter(str_detect(term, "^roh")) %>% mutate(log10 = -log10(p.value))
gwas31 %>% bind_rows() %>% filter(str_detect(term, "^roh")) %>% mutate(log10 = -log10(p.value))

snps <- gwas31 %>% map_chr(function(x) x[14, ]$term)
full_sample$map[full_sample$map$snp.name %in% snps, ]

library(ggplot2)
ggplot(annual_survival_gwas, aes(fill = as.factor(roh_oar3_OAR23_11362245), x = lamb)) + 
        geom_bar(position = "fill") +
        labs(y = "Proportion")
       