# This script is to plot genetic diversity and model estimates
# in regions around GWAS peaks. 
# Contains some exploratory analyses which are not in the manuscript.

library(tidyverse)
library(windowscanr)
library(data.table)
source("theme_simple.R")
library(snpStats)
library(viridis)
library(gt)

# prepare data -----------------------------------------------------------------
chr_info <- read_delim("data/chromosome_info_oar31.txt", "\t") %>% 
        .[-1, ] %>% 
        rename(chromosome = Part) %>% 
        mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
        mutate(chromosome = as.integer(chromosome)) %>% 
        filter(!is.na(chromosome))

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map <- full_sample$map 
table(full_sample$map$chromosome, useNA = "always")

# fitness data
load("data/survival_mods_data.RData")

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

# detailed look at significant regions -----------------------------------------
gwas_res <- read_rds("output/gwas_res_oar_roh_sep.rds")
#add map positions
gwas_out <- gwas_res %>% 
        rename(snp.name = snp) %>% 
        left_join(snps_map) %>% 
        rename(snp = snp.name)

# plot 1: Regional estimates, log10p and ROH counts ---------------------------- 

# filter top snp from every peak
top_snps <- gwas_out %>% 
        filter(state != "add") %>% 
        filter(-log10(p.value) > -log10(0.05/(2*39184))) %>%  # -log10(0.05/39149)
        mutate(peak = ifelse((chromosome == 3) & (position > 13000000 & position < 14000000), "3a", chromosome)) %>% 
        mutate(peak = ifelse((chromosome == 3) & (position > 170000000 & position < 180000000), "3b", chromosome)) %>% 
        group_by(peak) %>% 
        # mutate(pos_diff = position - lag(position, n = 1, default = position[1]))
        top_n(-1, wt = p.value) %>% 
        # two snps have same pvalue
        filter(!(snp %in% c("oar3_OAR10_85722487", "oar3_OAR3_177455437")))

# check MAF of top snps
maf <- col.summary(full_sample$genotypes[, top_snps$snp])$MAF
#write_delim(bind_cols(top_snps, maf) %>% rename(maf = `...11`), path = "data/top_snps_maf.txt")

# rerun models for positive effects
file_path <- "output/ROH/roh.hom"
roh_lengths <- fread(file_path)
snp_names <- top_snps$snp[c(5,7)]
pcs <- read_delim("data/ann_surv_pca.txt", " ", col_names = TRUE) %>% 
        mutate(id = as.character(id))
geno_sub <- as_tibble(as(full_sample$genotypes[, snp_names], 
                         Class = "numeric"),
                      rownames = "id")
snps_map_sub <- as_tibble(full_sample$map %>% filter(snp.name %in% snp_names))

ind_roh <- function(chromosome, position) {
        out <- (roh_lengths$CHR == chromosome) & (roh_lengths$POS1 <= position) & (roh_lengths$POS2 >= position)
}
rohs <- top_snps %>% 
        filter(estimate > 0) %>% 
        ungroup() %>% 
        select(chromosome, position) %>% 
        pmap(ind_roh) %>% 
        setNames(paste0("roh_", c(top_snps$snp[5], top_snps$snp[7]))) %>% 
        bind_cols()

surv_mod <- roh_lengths %>% 
        cbind(rohs) %>% 
        group_by(IID) %>% 
        summarise(roh_oar3_OAR19_36967641 = sum(roh_oar3_OAR19_36967641),
                  roh_oar3_OAR3_13845652 = sum(roh_oar3_OAR3_13845652)) %>% 
        rename(id =IID) %>% 
        mutate(id = as.factor(id)) %>% 
        right_join(annual_survival) %>% 
        right_join(geno_sub) %>% 
        left_join(pcs)

for (i in snp_names) {
        # dummy coding
        surv_mod[[paste0("roh_0_", i)]] <- as.numeric((surv_mod[[i]] == 0) & ( surv_mod[[paste0("roh_", i)]] == 1))
        surv_mod[[paste0("roh_2_", i)]] <- as.numeric(( surv_mod[[i]] == 2) & ( surv_mod[[paste0("roh_", i)]] == 1))
        surv_mod[[paste0("roh_", i)]] <- NULL
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

library(lme4)
# first positive snp
mod <- glmer(survival ~ 1 + sex + twin + age_std + age_std2 + froh_no_chr19 +
                     pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + roh_0_oar3_OAR19_36967641 +
                     roh_2_oar3_OAR19_36967641 + oar3_OAR19_36967641 + (1|birth_year) + (1|sheep_year) + (1|id),
             data = surv_mod, family = "binomial",
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
tidy(mod, conf.int = TRUE)

sum(surv_mod$roh_0_oar3_OAR19_36967641, na.rm = TRUE)

surv_mod %>% 
        filter(roh_0_oar3_OAR19_36967641 == 1) %>% 
        ggplot(aes(birth_year)) +
        geom_histogram(stat="count")

surv_mod_sub <- surv_mod %>% 
        filter(sex == "F")

mod2 <- glmer(survival ~ 1 + twin + sex + age_std + age_std2  +
                      pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + 
                      roh_0_oar3_OAR19_36967641 + roh_2_oar3_OAR19_36967641 + oar3_OAR19_36967641 + froh_no_chr19 +
                      #roh_0_oar3_OAR3_13845652 + roh_2_oar3_OAR3_13845652 + oar3_OAR3_13845652 + froh_no_chr3 +
                      (1|birth_year) + (1|sheep_year) + (1|id),
              data = surv_mod, family = "binomial",
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
out <- tidy(mod2, conf.int = TRUE)
print(out, n = 100)

ggplot(surv_mod, aes(as.factor(roh_0_oar3_OAR19_36967641), offspr_surv, col = sex)) +
        geom_jitter() +
        scale_y_log10()

surv_mod %>% group_by(id) %>% 
        summarise(lrs = sum(offspr_surv, na.rm = TRUE), 
                  roh_0_oar3_OAR19_36967641=mean(roh_0_oar3_OAR19_36967641, na.rm=TRUE),
                  sex) %>% 
        group_by(sex, roh_0_oar3_OAR19_36967641) %>% 
        summarise(mean(lrs, na.rm = TRUE))

# test LRS
surv_mod <- surv_mod %>% 
        mutate(olre = 1:nrow(surv_mod)) %>% 
        mutate(weight = scale(weight)) %>% 
        mutate(across(starts_with("pc"), ~(.x - mean(.x, na.rm=TRUE))/sd(.x, na.rm=TRUE)))

surv_mod_m <- surv_mod %>% 
        filter(sex == "M", age > 1) %>% 
        mutate(horn = as.factor(horn))

surv_mod_f <- surv_mod %>% 
        filter(sex == "F") %>% 
        mutate(horn = as.factor(horn))

# surv_mod %>% 
#         group_by(sex, horn) %>% 
#         tally()
# weight?
mod_lrs1 <- glmer(offspr_born ~ 1 + age_std + age_std2 + weight + horn + froh_no_chr3 +  roh_0_oar3_OAR3_13845652  + oar3_OAR3_13845652 +
                          (1|olre) + (1|birth_year) + (1|sheep_year) + (1|id),
                  data = surv_mod_m, family = "poisson",
                  control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod_lrs1)

mod_lrs2 <- glmer(offspr_born ~ 1 + age_std + age_std2 + froh_no_chr19 + horn + weight + oar3_OAR19_36967641 + 
                          roh_0_oar3_OAR19_36967641  + (1|olre) + (1|birth_year) + (1|sheep_year) + (1|id),
                  data = surv_mod_f, family = "poisson",
                  control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod_lrs2)



