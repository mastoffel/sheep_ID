# Run survival models 

# # run on server
library(lme4)
library(tidyverse)
library(broom.mixed)
source("theme_clean.R")
library(INLA)
library(snpStats)
library(data.table)
library(furrr)
library(brinla)

# data
#load("data/fitness_roh_df.RData") # formerly
load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")

# roh data
file_path <- "data/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Annual survival~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# survival data preprocessing
annual_survival <- fitness_data %>% 
        # filter na rows
        filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
        mutate(age2 = age^2,
               age_std = as.numeric(scale(age)),
               age2_std = as.numeric(scale(age2))) %>% 
        as.data.frame() 

# prepare pedigree for animal inla
sheep_ped_inla <- sheep_ped %>% 
        as_tibble() %>% 
        dplyr::rename(id = ID,
               mother = MOTHER,
               father = FATHER) %>% 
        mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "F", "888")) %>% 
        mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "M", "999")) %>% 
        mutate(father = ifelse(is.na(father), 0, father),
               mother = ifelse(is.na(mother), 0, mother)) %>% 
        mutate_if(is.character, list(as.numeric)) %>% 
        as.data.frame() 

# compute Ainverse and map
comp_inv <- AnimalINLA::compute.Ainverse(sheep_ped_inla)

# make Cmatrix and map
ainv <- comp_inv$Ainverse
ainv_map <- comp_inv$map
Cmatrix <- sparseMatrix(i=ainv[,1],j=ainv[,2],x=ainv[,3])

# link id to cmatrix
annual_survival$IndexA <- match(annual_survival$id, ainv_map[, 1])


# ~~~~~~~~~~~~~~~~~~
# add standardised and centered variables
annual_survival <- annual_survival %>% 
                        mutate(IndexA2 = IndexA) %>% 
                        mutate(froh_all_std = scale(froh_all),
                               froh_short_std = scale(froh_short),
                               froh_medium_std = scale(froh_medium),
                               froh_long_std = scale(froh_long),
                               froh_all_cent = froh_all - mean(froh_all),
                               froh_short_cent = froh_short - mean(froh_short),
                               froh_medium_cent = froh_medium - mean(froh_medium),
                               froh_long_cent = froh_long - mean(froh_long),
                               age_cent = age - mean(age),
                               age2_cent = age_cent ^ 2)

# model 1, FROH 
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
formula_surv <- as.formula(paste('survival ~ froh_all + sex + age_cent + age2_cent + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 #'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla <- inla(formula=formula_surv, family="binomial",
                 data=annual_survival, 
                 control.compute = list(dic = TRUE),
                 control.inla = list(correct = TRUE))

summary(mod_inla)
# bri.hyperpar.summary(mod_inla)
# hrtbl <- bri.hyperpar.summary(mod_inla)[4, 1] / sum(sum(bri.hyperpar.summary(mod_inla)[,1], pi^2/3))

# model 2, FROH long, medium and short
formula_surv2 <- as.formula(paste('survival ~ froh_long + froh_medium + froh_short + sex + age_cent + age2_cent + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 # 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla2 <- inla(formula=formula_surv2, family="binomial",
                 data=annual_survival , 
                 control.compute = list(dic = TRUE))

summary(mod_inla2)
saveRDS(list(mod_inla, mod_inla2), file = "output/inla_survival_models.rds")



# ~~~~~~~~~~~~~~~~~~~~~ Survival across different life-stages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

survival_inla <- function(ageclass) {
        dat_survival <- annual_survival %>% filter(age %in% c(ageclass)) %>% 
                        add_index_inla() %>% 
                        mutate(IndexA2 = IndexA)
        prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
        formula_surv <- as.formula(paste(#'survival ~ froh_long + froh_medium + froh_short + sex + twin + 1', 
                                         'survival ~ froh_all + sex + twin + 1', 
                                         'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                         'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                         #'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                         #'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                         'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
        mod_inla <- inla(formula=formula_surv, family="binomial",
                         data= dat_survival, 
                         control.compute = list(dic = TRUE))
}

all_inla_mods <- map(1:9, survival_inla)
saveRDS(all_inla_mods, file = "output/inla_survival_models_diff_ages_froh_all_new.rds")



# # extract all fixed effects
# ran_effs <- all_inla_mods %>% 
#         compact() %>% 
#         map(bri.hyperpar.summary) %>% 
#         map(as_tibble, rownames = "var") %>% 
#         bind_rows(.id = "age") %>% 
#         mutate(age = as.numeric(age) - 1)# %>% 
# #filter(str_detect(var, "froh"))
# names(ran_effs) <- c("age", "var", "mean", "sd", "lower", "median", "upper", "mode")
# ggplot(ran_effs, aes(age, mean, color = var)) + 
#         geom_point() +
#         geom_smooth(method="lm", se = FALSE) +
#         geom_errorbar(aes(ymin = lower, ymax = upper)) +
#         facet_wrap(~var, scales = "free_y")





# inbreeding depression changes across life?? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# lme4 models
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
annual_survival %>% filter(age > 0) -> annual_survival_1_and_over
mod <- glmer(survival ~ froh_long * age_std + froh_medium * age_std + froh_short * age_std + sex + twin + (1|mum_id) + (1|birth_year) + (1|sheep_year) + (1|id), 
             data = dat, family = "binomial",
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod)


####### inla ##############
annual_survival %>% filter(age > 0) -> annual_survival_1_and_over

prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
formula_surv <- as.formula(paste('survival ~ froh_all_cent * age_cent + age2_cent + sex + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 #'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))

mod_inla_std_inc_age0 <- inla(formula=formula_surv, family="binomial",
                 data=annual_survival, 
                 control.compute = list(dic = TRUE))

summary(mod_inla_std)
bri.hyperpar.summary(mod_inla)
hrtbl <- bri.hyperpar.summary(mod_inla)[4, 1] / sum(sum(bri.hyperpar.summary(mod_inla)[,1], pi^2/3))

# mod2
formula_surv2 <- as.formula(paste('survival ~ froh_long_cent * age_cent + froh_medium_cent * age_cent + froh_short_cent * age_cent + age2_cent + sex + twin + 1', 
                                  'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                  'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                  'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 # 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                  'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla2_std_inc_age0 <- inla(formula=formula_surv2, family="binomial",
                  data=annual_survival, 
                  control.compute = list(dic = TRUE))

summary(mod_inla2_std)

saveRDS(list(mod_inla_std_inc_age0, mod_inla2_std_inc_age0), file = "output/inla_survival_models_interaction_inc_age0.rds")


# repeat and exclude age = 0
annual_survival %>% filter(age > 0) -> annual_survival_1_and_over
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
formula_surv <- as.formula(paste('survival ~ froh_all_cent * age_cent + age2_cent + sex + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 #'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla_std_excl_age0 <- inla(formula=formula_surv, family="binomial",
                              data=annual_survival_1_and_over, 
                              control.compute = list(dic = TRUE))

# mod2
formula_surv2 <- as.formula(paste('survival ~ froh_long_cent * age_cent + froh_medium_cent * age_cent + froh_short_cent * age_cent + age2_cent + sex + twin + 1', 
                                  'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                  'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                  'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                  # 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                  'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla2_std_excl_age0 <- inla(formula=formula_surv2, family="binomial",
                               data=annual_survival_1_and_over, 
                               control.compute = list(dic = TRUE))

summary(mod_inla2_std)

saveRDS(list(mod_inla_std_excl_age0, mod_inla2_std_excl_age0), file = "output/inla_survival_models_interaction_excl_age0.rds")














# plot inla models
inla_mods_std <- list(mod_inla_std, mod_inla2_std)
# extract all fixed effects
fix_effs <- inla_mods_std %>% 
        map("summary.fixed") %>% 
        map(rownames_to_column, var = "var") %>% 
        bind_rows() %>% 
        filter(str_detect(var, c("froh"))) 

#filter(str_detect(var, "froh"))
names(fix_effs) <- c("var", "mean", "sd", "lower", "median", "upper", "mode", "kld")

library(boot)
library(wesanderson)

# only FROH
ggplot(fix_effs[c(1,2), ], aes(mean, as.factor(var))) + 
        geom_vline(aes(xintercept = 0), color = "lightgrey") +
        geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, color = "#3B4252") +
        geom_point(size = 2.5, shape = 21, fill = "darkgrey") +
        theme_clean() +
        scale_y_discrete(labels = c('FROH', "FROH * age")) +
        ylab("standardized predictors\n") + 
        xlab("\nFROH model estimate\n(logit scale)")
        # scale_color_manual(values = wes_palette("Darjeeling2")) +
        theme() -> p_roh_froh



# repeat the exercise with lme4

calc_id_est <- function(ageclass) {
        dat <- annual_survival %>% filter(age %in% c(ageclass)) 
        mod <- glmer(survival ~ froh_long + froh_medium + froh_short + sex + twin + (1|birth_year) + (1|sheep_year), 
                     data = dat, family = binomial(link = "probit"))
       # control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
}

all_mods <- map(0:7, calc_id_est)

all_mods %>% 
        map(tidy, conf.int = TRUE) %>% 
        bind_rows(.id = "age") %>% 
        filter(str_detect(term, "froh")) %>% 
        ggplot(aes(age, estimate)) +
        geom_point() +
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
        facet_wrap(~term, scales = "free_y")


