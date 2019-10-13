library(tidyverse)
library(snpStats)
source("theme_clean.R")
gwas_files <- list.files("output/gwas_full_roh/", pattern = "*.rds", full.names = TRUE)
all_gwas <- map(gwas_files, readRDS) %>% 
            flatten() %>% 
            map(function(x) x$result)
# check how many snps didn't work
not_working <- all_gwas %>% map(function(x) ifelse(is.null(x), 1, 0)) %>% 
                        unlist() %>% 
                        sum()

# filter empty elements
all_gwas <- all_gwas %>% purrr::compact()

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map <- full_sample$map 
table(full_sample$map$chromosome, useNA = "always")

# get roh pval
gwas_res <- map_df(all_gwas, function(x) x %>% .[c(4,5), ] %>% 
                           dplyr::select(term,estimate, p.value))
# gwas_res_roh <- map_df(all_gwas, function(x) x %>% .[c(5), ] %>% dplyr::select(term,estimate, p.value))
# gwas_res_snp <- map_df(all_gwas, function(x) x %>% .[c(4), ] %>% dplyr::select(term,estimate, p.value))

# check how often roh var couldnt be estimated
# gwas_roh <- gwas_res %>% filter(str_detect(term, "roh"))
# gwas_roh[which(!(gwas_snp$term %in% full_sample$map$snp.name)), ]

all_gwas[[11842]]
# put into df
gwas_full <- gwas_res %>%
        rename(snp.name = term) %>%
        # roh variable couldnt be estimated in some models, so the next
        # row was extracted from the tidy output
        filter(!str_detect(snp.name, "sd")) %>% 
        mutate(groups = ifelse(str_detect(snp.name, "roh"), "roh", "add")) %>% 
        mutate(snp.name = str_replace(snp.name, "roh_", "")) %>%
        left_join(snps_map) 
# 
qqman::qq(gwas_full$p.value)


# manhattan
## computing new x axis
gwas_roh <- gwas_full %>% 
                        filter(chromosome %in% c(4:26)) %>% 
                        group_by(groups) %>% 
                        arrange(chromosome, position) %>% 
                        dplyr::mutate(tmp = 1, cumsum.tmp = cumsum(tmp))
## calculating x axis location for chromosome label
med.dat <- gwas_roh %>% dplyr::group_by(groups, chromosome) %>% 
                dplyr::summarise(median.x = median(cumsum.tmp))


ggplot(data = gwas_roh) + 
        geom_point(aes(x = cumsum.tmp, y = -log10(p.value), color = chromosome %%2 == 0),
                   size = 0.01) + ## create alternate coloring
        geom_hline(yintercept = -log10(5e-06)) + ## add horizontal line
        scale_x_continuous(breaks = med.dat$median.x, labels = med.dat$chromosome) + ## add new x labels 
        guides(colour=FALSE) +  ## remove legend
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
        scale_color_manual(values = c(gray(0.5), gray(0))) +## instead of colors, go for gray
        theme_clean() +
        facet_wrap(groups~., nrow = 2)






# prepare genotypes for simpleM
snps_geno <- full_sample$map %>% 
                filter(chromosome == 2)
sheep_geno <- as(full_sample$genotypes[, snps_geno$snp.name], Class = "numeric")
missings <- rowSums(is.na(sheep_geno))
sheep_geno <- sheep_geno[missings < (0.01*ncol(sheep_geno)), ] # remove inds with more than 1% missing
sheep_geno[is.na(sheep_geno)] <- sample(c(0,1,2), 1)
sheep_geno_t <- t(sheep_geno)

library(data.table)
fwrite(sheep_geno_t, file = "data/geno_mat_simpleM_chr2.txt", col.names = FALSE, row.names = FALSE)

# library(synbreed)
# ?codeGeno
# impute_matrix <- function(sheep_geno) {
#         col_means <- colMeans(sheep_geno, na.rm = TRUE)
# }

gwas_roh %>% arrange(p.value)











library(data.table)
library(snpStats)
# check model
chr <- 15
# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# roh data
file_path <- "data/roh_nofilt_ram.hom"
roh_lengths <- fread(file_path)

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

snps_map_sub <- full_sample$map %>% 
        filter(chromosome == chr)  

# survival data
early_survival <- fitness_data %>% 
        dplyr::rename(birth_year = BIRTHYEAR,
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
        filter(age == 0) %>% 
        filter(!is.na(froh_all)) %>% 
        filter(!(is.na(birth_year) | is.na(mum_id))) %>% 
        mutate_at(c("id", "birth_year", "mum_id", "sex"), as.factor) %>% 
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
        as_tibble() %>% 
        mutate(MB = KB/1000) %>% 
        # sample_frac(0.01) %>% 
        mutate(index1 = as.numeric(match(SNP1, names(geno_sub))),
               index2 = as.numeric(index1 + NSNP - 1)) %>% 
        mutate(roh_snp_index = seq2(from = index1, to = index2)) %>% 
        group_by(IID) %>% 
        summarise(roh_snp_index = list(roh_snp_index), roh_lengths = list(MB)) %>% 
       # mutate(roh_snp_index = simplify_all(roh_snp_index)) %>% 
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
roh_list <- pmap(roh_snps_reord, function(id_roh, roh_snp_index, roh_lengths) {
        df <- as.matrix(t(as.numeric(c(id_roh, rep(0, ncol(geno_sub) - 1)))))
        tibble_iter <- tibble(roh_snp_index, roh_lengths)
        
        out <- pmap(tibble_iter, function(roh_snp_index, roh_lengths){
                df[, roh_snp_index] <<- roh_lengths
        })
        
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
early_survival_gwas <- early_survival %>% 
        dplyr::select(id, survival, sex, twin, birth_year, mum_id) %>% 
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

gwas_roh %>% arrange(p.value)

early_survival_gwas_top <- early_survival_gwas %>% select(id:mum_id, oar3_OAR15_63762123, roh_oar3_OAR15_63762123) 
early_survival_gwas_top_imp <- early_survival_gwas_top %>% mutate(survival = ifelse(is.na(survival), 0, 1))
early_survival_gwas_mod <- early_survival_gwas_top %>% 
                               # mutate(roh_roh_oar3_OAR15_63762123 = scale(roh_oar3_OAR15_63762123))
                                mutate(roh_oar3_OAR15_63762123 = ifelse(roh_oar3_OAR15_63762123 > 0, 1, 0))
library(lme4)
#gwas_f <- early_survival_gwas_top %>% filter(sex == "F")
mod <- glmer(survival ~ 1 + sex + twin + roh_oar3_OAR15_63762123 + oar3_OAR15_63762123 + (1|birth_year) + (1|mum_id),
             data = early_survival_gwas_mod, family = "binomial",
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod)
mod2 <- glmer(survival ~ 1 + sex + twin + oar3_OAR15_63762123 + (1|birth_year) + (1|mum_id),
             data = early_survival_gwas_top, family = "binomial")

early_survival_gwas_top %>% group_by(sex, roh_oar3_OAR15_63762123) %>% summarise(survival = mean(survival, na.rm = TRUE))

ggplot(early_survival_gwas_top, aes(roh_oar3_OAR15_63762123, mean(survival))) + geom_col()

early_survival_gwas_top %>% 
        filter(!is.na(sex) & !is.na(twin) & !is.na(birth_year) & !is.na(mum_id) & !is.na(roh_oar3_OAR15_63762123)) %>% 
        filter(!is.na(survival)) %>% 
        filter(!is.na(roh_oar3_OAR15_63762123))

summary(mod)
tidy(mod)
anova(mod, mod2)

ggplot(early_survival_gwas_top, aes(roh_oar3_OAR15_63762123, survival)) + geom_point()

# GWAS
run_gwas <- function(snp, data) {
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + ", snp, "+ ", paste0("roh_", snp), "+ (1|birth_year) + (1|mum_id)"))
        mod <- glmer(formula = formula_snp,
                     data = data, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod)
        out
}













# all survival data



library(data.table)
library(snpStats)
# check model
chr <- 15
# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# roh data
file_path <- "data/roh_nofilt_ram.hom"
roh_lengths <- fread(file_path)

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

snps_map_sub <- full_sample$map %>% 
        filter(chromosome == chr)  

# survival data
early_survival <- fitness_data %>% 
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
        filter(!(is.na(birth_year) | is.na(mum_id) | is.na(sheep_year))) %>% 
        mutate_at(c("id", "birth_year", "mum_id", "sex", "sheep_year"), as.factor) %>% 
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
        as_tibble() %>% 
        mutate(MB = KB/1000) %>% 
        # sample_frac(0.01) %>% 
        mutate(index1 = as.numeric(match(SNP1, names(geno_sub))),
               index2 = as.numeric(index1 + NSNP - 1)) %>% 
        mutate(roh_snp_index = seq2(from = index1, to = index2)) %>% 
        group_by(IID) %>% 
        summarise(roh_snp_index = list(roh_snp_index), roh_lengths = list(MB)) %>% 
        # mutate(roh_snp_index = simplify_all(roh_snp_index)) %>% 
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
roh_list <- pmap(roh_snps_reord, function(id_roh, roh_snp_index, roh_lengths) {
        df <- as.matrix(t(as.numeric(c(id_roh, rep(0, ncol(geno_sub) - 1)))))
        tibble_iter <- tibble(roh_snp_index, roh_lengths)
        
        out <- pmap(tibble_iter, function(roh_snp_index, roh_lengths){
                df[, roh_snp_index] <<- roh_lengths
        })
        
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
early_survival_gwas <- early_survival %>% 
        dplyr::select(id, survival, sex, twin, birth_year, mum_id, sheep_year, froh_all) %>% 
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

gwas_roh %>% arrange(p.value)

early_survival_gwas_top <- early_survival_gwas %>% select(id:mum_id,sheep_year, oar3_OAR15_63762123, roh_oar3_OAR15_63762123) 
early_survival_gwas_top_imp <- early_survival_gwas_top %>% mutate(survival = ifelse(is.na(survival), 0, 1))
early_survival_gwas_mod <- early_survival_gwas_top %>% 
        # mutate(roh_roh_oar3_OAR15_63762123 = scale(roh_oar3_OAR15_63762123))
        mutate(roh_oar3_OAR15_63762123 = ifelse(roh_oar3_OAR15_63762123 > 0, 1, 0))
library(lme4)
#gwas_f <- early_survival_gwas_top %>% filter(sex == "F") (1|mum_id) +
mod <- glmer(survival ~ 1 + sex + twin + oar3_OAR15_63762123 + roh_oar3_OAR15_63762123 + (1|mum_id) + (1|birth_year) + (1|sheep_year) + (1|id),
             data = early_survival_gwas_mod, family = "binomial",
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

summary(mod)
mod2 <- glmer(survival ~ 1 + sex + twin + oar3_OAR15_63762123 + (1|birth_year) + (1|mum_id),
              data = early_survival_gwas_top, family = "binomial")

early_survival_gwas_top %>% group_by(sex, roh_oar3_OAR15_63762123) %>% summarise(survival = mean(survival, na.rm = TRUE))

ggplot(early_survival_gwas_top, aes(roh_oar3_OAR15_63762123, mean(survival))) + geom_col()

early_survival_gwas_top %>% 
        filter(!is.na(sex) & !is.na(twin) & !is.na(birth_year) & !is.na(mum_id) & !is.na(roh_oar3_OAR15_63762123)) %>% 
        filter(!is.na(survival)) %>% 
        filter(!is.na(roh_oar3_OAR15_63762123))

summary(mod)
tidy(mod)
anova(mod, mod2)

ggplot(early_survival_gwas_top, aes(roh_oar3_OAR15_63762123, survival)) + geom_point()

# GWAS
run_gwas <- function(snp, data) {
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + ", snp, "+ ", paste0("roh_", snp), "+ (1|birth_year) + (1|mum_id)"))
        mod <- glmer(formula = formula_snp,
                     data = data, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod)
        out
}

