# try to write a function for marginal effects

# Run survival models 
library(lme4)
library(tidyverse)
library(broom.mixed)
source("theme_simple.R")
library(INLA)
library(snpStats)
library(data.table)
library(furrr)
library(brinla)
library(performance)
library(sjPlot)
library(ggeffects)
library(patchwork)
library(effects)
library(MCMCglmm)
library(MasterBayes)
# data
load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")
ped <- orderPed(sheep_ped)
# roh data
file_path <- "data/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Annual survival~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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


# lme4 -------------------------------------------------------------------------
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

# age interaction and lamb
mod <- glmer(survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
mod2 <- glmer(survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
             family = binomial(link = "probit"), data = annual_survival,
             control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

df1 <- get_model_data(mod, type = "eff", terms = c("froh_all10_cent [all]", "age_cent[-1.4, 2.6, 5.6]", "lamb")) 
df2 <- ggeffect(mod, type = "fe", terms = c("froh_all10_cent [0, 1]", "age_cent[-1.4, 2.6]"))

p1 <- plot_model(mod, type = "eff", terms = c("froh_all10_cent [all]", "age_cent[-1.4, 1.6, 2.6, 3.6, 7.6]"))
p2 <- plot_model(mod2, type = "eff", terms = c("froh_all10_cent [all]", "age_cent[-1.4, 1.6, 2.6, 3.6, 7.6]"))

# get predicted values
inv_logit <- function(x) exp(x) / (1+exp(x))
preds <- data.frame(df2)
preds

tidy(mod, effects = "fixed")
?new_data
df_new <- new_data(mod, terms = c("froh_all10_cent [-0.56, 0, 2.5]", "age_cent[-1.4, 2.6, 4.6]")) %>% 
          mutate(survival = NA)

newdata <- annual_survival
newdata$froh_all10_cent <- seq(min(annual_survival$froh_all10_cent), max(annual_survival$froh_all10_cent), length = nrow(annual_survival))
newdata$animal <- newdata$id
predict(mod, type = "response",  df_new, re.form = NA)

inv_logit(3.42 + 1.4*0.187)



mod_AS <- readRDS("output/AS_mod_MCMCglmm_cat_slice")
plot(mod_AS)
p1 <- predict(mod_AS, newdata = newdata, interval = "confidence", approx = "taylor2")
