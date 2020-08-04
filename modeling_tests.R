# Modelling inbreeding depression in survival
# Using binomial mixed (animal) models with logit link, annual survival as response
# main modeling package is INLA

library(lme4)
library(tidyverse)
library(broom.mixed)
source("theme_simple.R")
library(INLA) # Downloaded from http://www.r-inla.org/download
library(AnimalINLA) # Downloaded from http://www.r-inla.org/related-projects/animalinla
library(MCMCglmm)
library(sjPlot)
library(brinla)
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


mod_lme4 <- glmer(survival ~ froh_all10_cent * age_cent + froh_all10_cent * age_cent2 + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
                  family = binomial, data = annual_survival,
                  control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

tidy(mod_lme4)

library(sjPlot)
plot_model(mod_lme4, type = "pred", terms = c("froh_all10_cent", "age_cent2"))

annual_survival <- annual_survival %>% 
        mutate(life_stage = case_when(
                age < 2 ~ "lamb",
                age >= 2 & age < 4 ~ "adult",
                age >=4 ~ "old"
        )) %>% 
        mutate(life_stage = factor(life_stage, levels = c("lamb", "adult", "old")))

mod_lme4 <- glmer(survival ~ froh_all10_cent * life_stage + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
                  family = binomial, data = annual_survival,
                  control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
tidy(mod_lme4)
plot_model(mod_lme4, terms = c("froh_all10_cent", "life_stage"))
