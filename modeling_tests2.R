# Modelling inbreeding depression in survival
# Using binomial mixed (animal) models with logit link, annual survival as response
# main modeling package is INLA

library(lme4)
library(tidyverse)
library(broom.mixed)
source("theme_simple.R")
library(INLA) # Downloaded from http://www.r-inla.org/download
#library(AnimalINLA) # Downloaded from http://www.r-inla.org/related-projects/animalinla
library(MCMCglmm)
library(sjPlot)
library(performance)
library(Hmisc)
#library(brinla)
# data
load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")
ped <- sheep_ped

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Annual survival~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# survival data preprocessing
annual_survival <- fitness_data %>% 
        # filter na rows
        filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
        mutate(age_cent = age - mean(age, na.rm = TRUE),
               age_cent2 = age_cent^2,
               age_std = as.numeric(scale(age)),
               age_std2 = age_std^2,
               froh_all_cent = froh_all - mean(froh_all, na.rm = TRUE),
               # times 10 to estimate a 10% percent increase
               froh_all10 = froh_all * 10,
               froh_all10_cent = froh_all10 - mean(froh_all10, na.rm = TRUE),
               lamb = ifelse(age == 0, 1, 0),
               lamb_cent = lamb - mean(lamb, na.rm = TRUE),
               lamb = as.factor(lamb)) %>% 
        as.data.frame() 


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

#levels(annual_survival$sex) <- c("M", "F")
annual_survival <- annual_survival %>% 
        mutate(lamb = as.numeric(factor(lamb, levels = c(1, 0)))) %>% 
        mutate(life_stage = case_when(
                age == 0 ~ "lamb",
                age > 0 & age < 2 ~ "yearling",
                age >= 2 ~ "adult"
        )) %>% 
        mutate(life_stage = factor(life_stage, levels = c( "yearling","lamb", "adult"))) %>% 
        mutate(life_stage2 = case_when(
                age == 0 ~ "lamb",
                age > 0 & age < 3 ~ "juvenile",
                age >= 3 ~ "adult"
        )) %>% 
        mutate(life_stage2 = factor(life_stage2, levels = c("juvenile", "lamb", "adult")))

fit1 <- glmer(survival ~ froh_all10_cent * age_cent + froh_all10 * lamb + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

fit2b <- glmer(survival ~ froh_all10_cent * age_cent + age_cent2 + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

fit3 <- glmer(survival ~ froh_all10_cent * life_stage + sex + twin + (1|birth_year) +  (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

fit4 <- glmer(survival ~ froh_all10_cent * life_stage2 + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

fit5 <- glmer(survival ~ froh_all10_cent * lamb + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

perf <- compare_performance(fit1, fit2, fit2b, fit3, fit4, fit5)
plot(perf)
check_model(fit3, check = c("vif", "qq", "normality", "ncv", "homogeneity", "outliers"))

summary(fit3)

tab_model(fit1)

mod_lme4_2class <- glmer(survival ~ froh_all10 * life_stage + sex + twin +  (1|id) + (1|sheep_year),
                         family = binomial, data = annual_survival,
                         control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

# mod_lme42 <- glmer(survival ~ froh_all10 * age_cent + lamb + sex + twin + (1|sheep_year) + (1|id),
#                   family = binomial, data = annual_survival,
#                   control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
# mod_lme43 <- glmer(survival ~ froh_all10 * age_cent + froh_all10 * lamb + sex + twin + (1|sheep_year) + (1|id),
#                    family = binomial, data = annual_survival,
#                    control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
#mod_lme44 <- 
compare_performance(mod_lme4, mod_lme4_2class)
plot_model(fit4, type = "pred", terms = c("froh_all10_cent[all]","life_stage2"))


surv_per_F <- function(iter) {
        out <- annual_survival %>% 
                dplyr::mutate(binned_froh = cut2(froh_all10, g = 4)) %>% 
                dplyr::group_by(life_stage, binned_froh) %>% 
                dplyr::sample_n(size = nrow(annual_survival), replace = TRUE) %>% 
                dplyr::select(froh_all10, binned_froh, survival) %>% 
                dplyr::summarise(binned_survival = mean(survival)) 
        out
}

all_boot <- map_df(1:10, surv_per_F)

ggplot(all_boot , aes(binned_froh, binned_survival, color = life_stage)) +
        geom_point()
