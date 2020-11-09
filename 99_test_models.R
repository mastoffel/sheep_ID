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
library(performance)
library(ggeffects)
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
               lamb = as.factor(lamb),
               log_age = log(age + 0.00001)) %>% 
        as.data.frame() 


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

mod1 <- glmer(survival ~ froh_all10_cent + twin + sex + (1|sheep_year) + (1|birth_year) + (1|id),
              family = binomial(link = 'logit'), data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

summary(mod1)
tidy(mod1, conf.int = TRUE)
ggpredict(mod1, terms = c("froh_all10_cent [all]", "age [0, 1, 2, 3, 4, 8,9,10]"), type="re", ci.lvl = 0.95) %>% plot()
plot_model(mod1)


mod2 <- glmer(survival ~ froh_all10_cent * age + froh_all10_cent * lamb + twin + sex + (1|sheep_year) + (1|birth_year) + (1|id),
              family = binomial(link = 'logit'), data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

tidy(mod2, conf.int = TRUE)
sims <- simulate(mod2, nsim = 1000, re.form = NULL)
emp_ratio <- sum(annual_survival$survival) / nrow(annual_survival)
hist(colSums(sims) / nrow(sims), breaks = 100) + abline(v = emp_ratio)
mean(colSums(sims) / nrow(sims))
emp_ratio
plot_model(mod2)
binned_residuals(mod2, n_bins = 100)
ggsave("figs/surv_mod_binned_resid.jpg", width = 5, height = 3)
binned_residuals(mod2, term = "log_age")
check_model(mod2)
plot(mod2)

mod3 <- glmer(survival ~ froh_all10_cent * as.factor(age) + twin + sex + (1|sheep_year) + (1|birth_year) + (1|id),
              family = binomial(link = 'logit'), data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
binned_residuals(mod3)

simulationOutput <- simulateResiduals(fittedModel = mod2)
plot(simulationOutput, asFactor = T)
plotResiduals(mod2)

annual_survival_mod <- annual_survival %>% mutate(age_cent = as.factor(age_cent))
mod3 <- glmer(survival ~ froh_all10_cent * age_cent + twin + sex + (1|sheep_year) + (1|birth_year) + (1|id),
              family = binomial(link = 'logit'), data = annual_survival_mod,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

tidy(mod3, conf.int = TRUE)
sims <- simulate(mod3, nsim = 1000, re.form = NULL)
emp_ratio <- sum(annual_survival$survival) / nrow(annual_survival)
hist(colSums(sims) / nrow(sims), breaks = 100) + abline(v = emp_ratio)
mean(colSums(sims) / nrow(sims))
emp_ratio
