library(tidyverse)
library(data.table)
library(AnimalINLA)
library(rlang)
library(furrr)
library(lme4)
library(broom.mixed)
library(partR2)
library(magrittr)
# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# roh data
file_path <- "data/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path) 

chr_info <- read_delim("../sheep/data/sheep_genome/chromosome_info_ram.txt", "\t") %>% 
        .[-1, ] %>% 
        rename(chromosome = Part) %>% 
        mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
        mutate(chromosome = as.integer(chromosome)) %>% 
        filter(!is.na(chromosome)) %>% 
        rename(chr = chromosome)

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
        filter(!(is.na(survival) | is.na(froh_all) | is.na(birth_year))) %>% 
        filter(!(is.na(sheep_year))) %>% 
        mutate_at(c("id", "birth_year", "sex", "sheep_year"), as.factor) %>% 
        mutate(age2 = age^2) %>% 
        mutate(age_std = as.numeric(scale(age)),
               age2_std = as.numeric(scale(age2))) %>% 
        as.data.frame() 


froh_per_chr <- roh_lengths %>% 
        rename(chr = CHR) %>% 
        group_by(chr, IID) %>% 
        summarise(KB_sum = sum(KB)) %>% 
        left_join(chr_info, by = "chr") %>% 
        mutate(length_chr_KB = Length / 1000) %>% 
        mutate(froh_chr = KB_sum / length_chr_KB) %>% 
        rename(id = IID) %>% 
        mutate(id = as.factor(id)) %>% 
        select(chr, id, froh_chr) %>% 
        pivot_wider(id_cols = id, names_from = chr, names_prefix = "froh_chr_",
                    values_from = froh_chr)

annual_survival %<>% 
        left_join(froh_per_chr, by = "id")
        
# lme4 
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

# run model for each chromosome
chrom_part <- function(chr, annual_survival) {
        
        formula_surv <- reformulate(c(paste0("froh_chr_",chr), "sex", "age_std", 
                                      "age2_std", "twin", paste0("FROH_no_chr", chr),
                                      "(1|birth_year)", "(1|sheep_year)", "(1|id)"),
                                    response = "survival")
        mod <- glmer(formula_surv, data = annual_survival, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
       
}

mod_per_chr <- map(1:26, chrom_part, annual_survival)

# get tidy output
mod_per_chr %>% 
        map(tidy, conf.int = TRUE) %>% 
        bind_rows(.id = "chr") %>% 
        filter(str_detect(term, "froh_chr")) %>% 
        #filter(grepl('FROH_no_chr', term)) %>% 
        mutate(chr = as.numeric(chr)) -> froh_b_per_chr
        #ggplot(aes(chr, estimate)) + geom_point()
        
# get average ROH length per chromosome
roh_lengths %>% 
        #filter(KB > 5000) %>% 
        group_by(CHR, IID) %>% 
        summarise(sum_roh_MB = sum(KB) / 1000) %>% 
        ungroup() %>% 
        group_by(CHR) %>% 
        summarise(mean_roh_MB = mean(sum_roh_MB),
                  sd_roh_MB = sd(sum_roh_MB)) %>% 
        rename(chr = CHR) %>% 
        right_join(froh_b_per_chr, by = "chr") %>% 
        #right_join(froh_r2_per_chr, by = "chr") %>% 
        left_join(chr_info, by = "chr") %>% 
        mutate(chr_length = Length / 1000000,
               roh_length_ratio = mean_roh_MB/chr_length) -> chrom_part_df

ggplot(chrom_part_df, aes(roh_length_ratio, estimate, color = as.factor(chr))) + 
        geom_point() 
        #geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) 
       

















# inla prep
sheep_ped_inla <- sheep_ped %>% 
        as_tibble() %>% 
        rename(id = ID,
               mother = MOTHER,
               father = FATHER) %>% 
        mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "F", "888")) %>% 
        mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "M", "999")) %>% 
        #filter(!is.na(id)) %>% 
        mutate(father = ifelse(is.na(father), 0, father)) %>% 
        mutate(mother = ifelse(is.na(mother), 0, mother)) %>% 
        mutate_if(is.character, list(as.numeric)) %>% 
        as.data.frame() 

comp_inv <- AnimalINLA::compute.Ainverse(sheep_ped_inla)
ainv <- comp_inv$Ainverse
ainv_map <- comp_inv$map
Cmatrix <- sparseMatrix(i=ainv[,1],j=ainv[,2],x=ainv[,3])

add_index_inla <- function(dat) {
        Ndata <- dim(dat)[1]
        dat$IndexA <- rep(0, times = Ndata)
        for(i in 1:Ndata) dat$IndexA[i] <- which(ainv_map[,1]==dat$id[i])
        dat
}

# add some more stuff to data 
annual_survival <- add_index_inla(annual_survival)
annual_survival <- annual_survival %>% 
        mutate(IndexA2 = IndexA) %>% 
        mutate(froh_all_std = scale(froh_all),
               froh_short_std = scale(froh_short),
               froh_medium_std = scale(froh_medium),
               froh_long_std = scale(froh_long))

# modeling
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))

chrom_part <- function(chr, annual_survival, prec_prior, Cmatrix) {
        # add froh for a chromosome
        dat_mod <- annual_survival %>% 
                mutate(froh_chr = froh_all - .data[[paste0("FROH_no_chr", chr)]])
        formula_surv <- as.formula(paste('survival ~ 1 + froh_chr + sex + age_std + age2_std + twin ',
                                         paste0("FROH_no_chr", chr), 
                                         'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                         'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                         'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                         # 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))',
                                         'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
        mod_inla <- inla(formula = formula_surv, family="binomial",
                         data = dat_mod,
                         control.compute = list(dic = TRUE))
}

#plan(multiprocess, workers = 4)
all_chr <- map(1, chrom_part, annual_survival, prec_prior, Cmatrix)


