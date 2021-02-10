# Modelling inbreeding depression in survival
# Using binomial mixed (animal) models with logit link, annual survival as response
# main modeling package is INLA

library(lme4)
library(tidyverse)
library(broom.mixed)
source("theme_simple.R")
library(MCMCglmm)
library(sjPlot)
library(performance)
library(Hmisc)
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
        mutate(life_stage = case_when(
                age == 0 ~ "lamb",
                age > 0 & age <= 2 ~ "early_life",
                age > 2 & age <= 4 ~ "mid_life",
                age > 4 ~ "late_life",
        )) %>% 
        mutate(life_stage = factor(life_stage, levels = c( "lamb","early_life", "mid_life", "late_life"))) %>% 
        mutate(life_stage1 = case_when(
                age == 0 ~ "lamb",
                age > 0 & age <= 2 ~ "early_life",
                age >= 3 ~ "late_life"
        )) %>% 
        mutate(life_stage1 = factor(life_stage1, levels = c("early_life", "lamb", "late_life"))) %>% 
        mutate(life_stage2 = case_when(
                age == 0 ~ "lamb",
                age > 0 & age <= 3 ~ "early_life",
                age >= 4 ~ "late_life"
        )) %>% 
        mutate(life_stage2 = factor(life_stage2, levels = c("early_life", "lamb", "late_life"))) 


as1 <- annual_survival %>% 
        filter(age > 2)

fit1 <- glmer(survival ~ froh_all10_cent * age_cent + sex + twin + (1|birth_year) + (1|sheep_year),
              family = binomial, data = as1,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
tidy(fit1, conf.int = TRUE)



fit1 <- glmer(survival ~ froh_all10_cent * (age_cent + lamb) + sex + twin + (1|birth_year) + (1|sheep_year),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

fit2 <- glmer(survival ~ froh_all10_cent * age_cent + lamb + sex + twin + (1|birth_year) + (1|sheep_year) ,
               family = binomial, data = annual_survival,
               control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

fit3 <- glmer(survival ~ froh_all10_cent * life_stage1 + sex + twin + (1|birth_year) + (1|sheep_year),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

fit4 <- glmer(survival ~ froh_all10_cent * life_stage2 + sex + twin + (1|birth_year) + (1|sheep_year),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

fit5 <- glmer(survival ~ froh_all10_cent * life_stage + sex + twin + (1|birth_year) + (1|sheep_year),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

perf <- compare_performance(fit1, fit2, fit3, fit4, fit5,
                            metrics = c("AIC", "BIC", "RMSE", "SIGMA", "LOGLOSS", "SCORE"))
perf %>% 
        arrange(AIC, BIC)

tidy(fit5, conf.int = TRUE)
check_model(fit5, check = c("vif", "qq", "normality", "ncv", "homogeneity", "outliers"))


plot_model(fit5, type = "pred", terms = c("froh_all10_cent[all]", "life_stage"),
           ci.lvl = 0.95)
plot_model(fit3, type = "pred", terms = c("froh_all10_cent[all]", "life_stage1"),
           ci.lvl = 0.95)

plot_model(fit1, type = "pred", terms = c("froh_all10_cent[all]", "age_cent"),
           ci.lvl = 0.80)

# plot raw data
surv_per_F <- function(iter) {
        out <- annual_survival %>% 
                #filter(sex == "F") %>% 
                #filter(twin == 0) %>% 
                dplyr::mutate(binned_froh = cut2(froh_all10, cuts = c(1.8, 2, 2.2, 2.4, 2.6, 2.8, 5.1))) %>% 
                #dplyr::mutate(binned_froh = cut2(froh_all10, cuts = c(1.8, 2.1, 2.4, 2.7, 3.0, 5.1))) %>% 
               # dplyr::mutate(binned_froh = cut2(froh_all10, m = 700)) %>% 
                #dplyr::mutate(binned_froh = cut2(froh_all10, g = 6)) %>% 
                #dplyr::mutate(binned_froh = cut(froh_all10, breaks = c(1.8, 2.4, 3, 5.1))) %>% 
                dplyr::group_by(life_stage, binned_froh) %>% 
                #filter(length(unique(id)) >= 20) %>% 
                # filter by unique individuals
                filter(n() >= 20) %>% 
                #dplyr::sample_frac(0.9) %>% 
                dplyr::sample_n(size = nrow(annual_survival), replace = TRUE) %>% 
                dplyr::select(froh_all10, binned_froh, survival) %>% 
                dplyr::summarise(binned_survival = mean(survival)) 
        out
}

all_boot <- map_df(1:10, surv_per_F)
all_boot <- all_boot %>% 
        mutate(life_stage1 = fct_relevel(life_stage1, "lamb", after = 0))

p <- ggplot(all_boot , aes(binned_froh, binned_survival, fill = life_stage)) +
        geom_jitter(width = 0.05, size = 2, alpha = 0.7, shape = 21, stroke=0.2) +
        scale_fill_viridis_d("Life Stage (age)", 
                             labels = c("Lamb (0)", "Early life (1,2)",
                                        "Mid life (2,3)",
                                        "Late life (4+)"), option = "D") +
        scale_y_continuous(breaks = seq(0.3, 1, 0.1), labels = seq(30, 100, 10)) +
        scale_x_discrete(labels = c("[0.18,0.20)", "[0.20,0.22)", "[0.22,0.24)", 
                                    "[0.24,0.26)", "[0.26,0.28)", "[0.28,0.50)"),
                         guide = guide_axis(n.dodge = 2)) +
        theme_simple() +
        xlab(expression(F[ROH]~class)) +
        ylab("% annual survivors")
        #geom_line(mapping = aes(group =age), size = 0.2, alpha = 1) 
p
ggsave("figs/annual_survivors.jpg", width = 6, height =3.5)
