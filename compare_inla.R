# this script is to test and compare inla

library(tidyverse)
library(INLA)
library(brinla)
library(lme4)
library(broom.mixed)


load("data/survival_mods_data.RData") 

# survival data preprocessing
annual_survival <- fitness_data %>% 
        # filter na rows
        filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
        mutate(age2 = age^2,
               age_std = as.numeric(scale(age)),
               age2_std = as.numeric(scale(age2))) %>% 
        as.data.frame() 

# INLA
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
formula_surv <- as.formula(paste('survival ~ froh_all + sex + age_std + age2_std + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(id, model = "iid", hyper = list(prec = prec_prior))', sep = " + "))

mod_inla <- inla(formula=formula_surv, family="binomial",
                 data=annual_survival, 
                 control.compute = list(dic = TRUE))
# with copula-based correction to laplace approximation
mod_inla2 <- inla(formula=formula_surv, family="binomial",
                 data=annual_survival, 
                 control.compute = list(dic = TRUE),
                 control.inla = list(correct = TRUE))

# lme4
mod_lme4 <- glmer(survival ~ 1 + froh_all + sex + age_std + age2_std + twin + (1|birth_year) + (1|sheep_year) + (1|id),
                  data = annual_survival, family = "binomial",
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

tidy_lme4 <- tidy(mod_lme4, conf.int = TRUE, conf.level = 0.95)
summary(mod_lme4)

# fixed effect estimates v similar
mod_inla$summary.fixed %>% 
        rownames_to_column(var = "term") %>% 
        rename(estimate = mean,
               conf.low = `0.025quant`,
               conf.high = `0.975quant`) %>% 
        bind_rows(tidy_lme4, .id = "mods") %>% 
        filter(term != "sd__(Intercept)") %>% 
        ggplot(aes(estimate, term, color = mods)) +
                geom_point(size = 2) + 
                geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                               height = 0.5) 


# roh split models
formula_surv2 <- as.formula(paste('survival ~ froh_long + froh_medium + froh_short + sex + age_std + age2_std + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(id, model = "iid", hyper = list(prec = prec_prior))', sep = " + "))

# with copula-based correction to laplace approximation
mod_inla2<- inla(formula=formula_surv2, family="binomial",
                  data=annual_survival, 
                  control.compute = list(dic = TRUE),
                  control.inla = list(correct = TRUE))

# lme4
mod_lme42 <- glmer(survival ~ froh_long + froh_medium + froh_short + sex + age_std + age2_std + twin + (1|birth_year) + (1|sheep_year) + (1|id),
                  data = annual_survival, family = "binomial",
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

tidy_lme42 <- tidy(mod_lme42, conf.int = TRUE, conf.level = 0.95)
summary(mod_lme4)

# fixed effect estimates v similar
mod_inla2$summary.fixed %>% 
        rownames_to_column(var = "term") %>% 
        rename(estimate = mean,
               conf.low = `0.025quant`,
               conf.high = `0.975quant`) %>% 
        bind_rows(tidy_lme42, .id = "mods") %>% 
        filter(term != "sd__(Intercept)") %>% 
        ggplot(aes(estimate, term, color = mods)) +
        geom_point(size = 2) + 
        geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                       height = 0.5) 

# with interaction


