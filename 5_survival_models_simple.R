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

# random_slopes
mod1 <- glmer(survival ~ froh_all10_cent + sex + twin + (1 + froh_all10_cent|age) + (1|birth_year) + (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
# plot 
ggpredict(mod1, type = "re", terms = c("froh_all10_cent", "age [0,1,2,3]")) %>% plot()
performance(mod1)
plot(-0.9 - ranef(mod6)$age[,2])

library(boot)
# age interaction and lamb
mod2 <- glmer(survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
saveRDS(mod2, file = "output/survival_mod_full_lme4")
mod2 <- readRDS("output/survival_mod_full_lme4")
plot_model(mod2)
summary(mod2)
inv.logit(-1.197)
odds_to_prob <- function(x) exp(x) / (1+exp(x))
get_model_data(mod2, type = "est")
check_model(mod2)
performance(mod2)
plot_model(mod2, type = "eff", terms = c("froh_all10_cent [all]", "age_cent[-1.4, 1.6, 2.6, 3.6, 7.6]"))
plot_model(mod2, type = "eff", terms = c("froh_all10_cent [all]",  "age_cent[-2.4, -1.4]", "lamb"))
df1 <- get_model_data(mod2, type = "eff", 
                      terms = c("froh_all10_cent [all]", "age_cent[-1.4, 2.6, 5.6]", "lamb")) %>% 
        as_tibble() %>% 
        filter(facet == 0)
df2 <- get_model_data(mod2, type = "eff", 
                      terms = c("froh_all10_cent [all]",  "age_cent[-2.4, -1.4]", "lamb")) %>% 
        as_tibble() %>% 
        filter(facet == 1) %>% 
        filter(group == -2.4)
df_full <- bind_rows(df1, df2) %>% 
        mutate(group = as.character(as.numeric(group) + 2.4))

ggplot() +
        geom_line(data = df_full, aes(x = x, predicted, color = group), size = 1) +
        geom_ribbon(data= df_full, aes(x=x, ymin=conf.low, ymax=conf.high, fill = group, color = group), alpha= 0.2,
                    linetype = 2, size = 0.1) +
        scale_color_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        scale_fill_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) 




# MCMCglmm ---------------------------------------------------------------------
prior1 <-list(R=list(V=1, fix=1), 
             G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), 
                    G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                    G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                    G4=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)
                    ))
annual_survival$animal <- annual_survival$id

# ~ 8 minutes for 5000/10/1000
mod_AS <- MCMCglmm(survival~froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin, 
                   random=~birth_year + sheep_year + id + animal,
                   data=annual_survival,
                   family="threshold", 
                   #family="categorical",
                   trunc=TRUE,
                   prior=prior1, 
                   pr=TRUE,
                   nitt=200000,thin=80,burnin=20000,
                   #nitt=5000,thin=10,burnin=1000,
                   pedigree=ped)

saveRDS(mod_AS, file = "output/AS_mod_MCMCglmm")

mod_AS2 <- MCMCglmm(survival~froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin, 
                   random=~birth_year + sheep_year + id + animal,
                   data=annual_survival,
                   family="threshold", 
                   family="categorical",
                   #trunc=TRUE,
                   prior=prior1, 
                   pr=TRUE,
                   nitt=200000,thin=80,burnin=20000,
                   #nitt=5000,thin=10,burnin=1000,
                   pedigree=ped)

saveRDS(mod_AS2, file = "output/AS_mod_MCMCglmm_cat")


# load mcmc
mod_AS <- readRDS("output/AS_mod_MCMCglmm_logit_08m_nosclice")


plot(mod_AS)
tidy(mod_AS, conf.int = TRUE)
tidy(mod2, conf.int = TRUE)
summary(mod_AS)

test <- predict(mod_AS)
tidy(mod_AS, conf.int = TRUE)
tidy(mod_AS, conf.int = TRUE)
predict(mod_AS, type = "response")[1:5, ]

out <- predict(mod_AS2, type = "response", newdata = annual_survival,
               interval = "prediction")[1:100, ]

ggemmeans(mod_AS, terms = "froh_all10_cent")

df1 <- ggemmeans(mod_AS, terms = c("froh_all10_cent [all]", "age_cent[-1.4, 2.6, 5.6]", "lamb")) %>% 
        as_tibble() %>% 
        filter(facet == 0)
df2 <- ggemmeans(mod_AS, terms = c("froh_all10_cent [all]",  "age_cent[-2.4, -1.4]", "lamb")) %>% 
        as_tibble() %>% 
        filter(facet == 1) %>% 
        filter(group == -2.4)
df_full <- bind_rows(df1, df2) %>% 
        mutate(group = as.character(as.numeric(group) + 2.4))

ggplot() +
        #geom_jitter(data = annual_survival, aes(froh_all10_cent, y = -0.1), height = 0.05,
        #            alpha = 0.2, size = 3) +
        #geom_point(data = df_full, aes(x = froh_all10_cent, fit, color = age_cent)) +
        geom_line(data = df_full, aes(x = x, predicted, color = group), size = 1) +
        geom_ribbon(data= df_full, aes(x=x, ymin=conf.low, ymax=conf.high, fill = group, color = group), alpha= 0.2,
                    linetype = 2, size = 0.1) +
        scale_color_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        scale_fill_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) 




# INLA
# prepare pedigree for animal inla
sheep_ped_inla <- ped %>% 
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
annual_survival <- annual_survival %>% 
        mutate(IndexA2 = IndexA) 

#annual_survival_test <- annual_survival %>% sample_frac(0.2)
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
# model 2, with lamb 
formula_surv <- as.formula(paste('survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 #'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
control.family1 = list(control.link=list(model="logit"))
mod_inla <- inla(formula=formula_surv, family="binomial",
                 data=annual_survival, 
                 control.family = control.family1,
                 control.compute = list(dic = TRUE),
                 control.inla = list(correct = TRUE))
summary(mod_inla)
tidy(mod_AS)

mod_inla$summary.fixed
tidy(mod_AS, conf.int = TRUE)
# 
# inv.phylo <- MCMCglmm::inverseA(ped)
# A <- solve(inv.phylo$Ainv)
# rownames(A) <- rownames(inv.phylo$Ainv)
# 
# library(brms)
# 
# model_simple <- brm(
#         survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin +(1|animal) + (1|id) + (1|birth_year) + (1|sheep_year),
#         data = annual_survival, 
#         family = bernoulli(), cov_ranef = list(animal = A),
#         chains = 1, cores = 1, iter = 100

