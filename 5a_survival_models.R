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


# INLA -------------------------------------------------------------------------
# format pedigree
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
# set prior
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
# model
formula_surv <- as.formula(paste('survival ~ froh_all10_cent * age + froh_all10_cent * lamb + sex + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                # 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))

control.family1 = list(control.link=list(model="logit"))
mod_inla <- inla(formula=formula_surv, family="binomial",
                 data=annual_survival, 
                 control.family = control.family1,
                 control.compute = list(dic = TRUE, cpo=TRUE, waic = TRUE, po=TRUE, config=TRUE),
                 control.inla = list(correct = TRUE))

saveRDS(mod_inla, file = "output/AS_mod_oar2.rds")
mod_inla <- readRDS("output/AS_mod_oar2.rds")

mod_inla$summary.fixed
summary(mod_inla)
summary(mod_inla$cpo)
pit_adj <- mod_inla$cpo$pit - 0.5 * mod_inla$cpo$cpo
hist(pit_adj)
## posterior predictive check
#bri.fixed.plot(mod_inla)
sim <- inla.posterior.sample(100, mod_inla)
sim[[1]]$latent
library(INLAutils)
ggplot_inla_residuals(mod_inla, annual_survival$survival)

# POISSON model for inbreeding load --------------------------------------------
# to calculate lethal equivalents according to Nietlisbach et al 2019
# this is a simplified model for better cross-study comparability of the inbreeding load
formula_surv <- as.formula(paste('survival ~ froh_all10_cent + age_std + age_std2 + sex + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 #'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla_pois <- inla(formula=formula_surv, family="poisson",
                       data=annual_survival, 
                       control.compute = list(dic = TRUE, config=TRUE)
)

saveRDS(mod_inla_pois, file = "output/AS_mod_INLA_oar_poisson_full.rds")
summary(mod_inla_pois)

# inbreeding load (lethal equivalents) 2B
#  froh_all10_cent  4.57 [2.61, 6.55]  
mod_inla_pois$summary.fixed %>% 
        rownames_to_column(var = "pred") %>% 
        as_tibble() %>% 
        dplyr::select(c(1,2,4,6)) %>% 
        setNames(c("pred", "est", "lower_ci", "upper_ci")) %>% 
        filter(pred == "froh_all10_cent") %>% 
        mutate(across(is.numeric, function(x) abs(2 * x/0.10)))


# INLA models across age classes
# ~~~~~~~~~~~~~~~~~~~~~ Survival across different life-stages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# here, we fit one model per age class until the age of 10 (sample size becomes too small afterwards)
survival_inla <- function(ageclass) {
        dat_survival <- annual_survival %>% filter(age %in% c(ageclass)) 
        prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
        formula_surv <- as.formula(paste(#'survival ~ froh_long + froh_medium + froh_short + sex + twin + 1', 
                'survival ~ froh_all10_cent + sex + twin + 1', 
                'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                #'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                #'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
        mod_inla <- inla(formula=formula_surv, family="binomial",
                         data= dat_survival, 
                         control.inla = list(correct = TRUE),
                         control.compute = list(dic = TRUE))
}

all_inla_mods <- map(0:10, survival_inla)
saveRDS(all_inla_mods, file = "output/inla_survival_models_diff_ages_froh_all_oar.rds")
all_inla_mods <- readRDS(file = "output/inla_survival_models_diff_ages_froh_all_oar.rds")

# plot for supplementary
inla_plot <- map_df(all_inla_mods, function(x) as_tibble(x$summary.fixed[2, ]), .id = "age") %>% 
        mutate(age = as.numeric(age) - 1) %>% 
        filter(age < 10) %>% 
        mutate(age = fct_inorder(as.factor(age))) 
   
sample_size <- annual_survival %>% group_by(age) %>% tally() %>% mutate(age = as.character(age))
inla_plot <- inla_plot %>% left_join(sample_size, by = "age")

ggplot(inla_plot, aes(age, mean)) +
        geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`),  width = 0.5) +
        geom_point(size = 3, shape = 21, col = "#4c566a", fill = "#eceff4", # "grey69"
                   alpha = 1, stroke = 0.7) +
        geom_text(aes(y = `0.975quant` + 0.3, label = n), size = 3) + 
        xlab("Age class") + 
        ylab(expression(beta[F[ROH]])) +
        geom_hline(yintercept = 0, linetype='dashed', colour =  "#4c566a", size = 0.3) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) -> p_froh_across_life
p_froh_across_life
ggsave("figs/sup_froh_across_life.jpg", width = 5, height = 3)










# Modeling with lme4 and MCMCglmm ----------------------------------------------
# here, we compared INLA results to those from lme4 and MCMCglmm
# code might be a bit messy
# main conclusion: fixed and random effect estimates highly similar

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


levels(annual_survival$sex) <- c("M", "F")
mod_lme4 <- glmer(survival ~ froh_all * age_cent + froh_all * lamb + sex + twin + (1|sheep_year) + (1|id),
                  family = binomial, data = annual_survival,
                  control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))


summary(mod_lme4)
df <- as_tibble(emmeans::ref_grid(mod_lme4))
predict(emmeans::ref_grid(mod_lme4))
#plot_model(mod_lme4, type = "pred", terms = c("froh_all10_cent", "age_cent", "lamb"))

df1 <- get_model_data(mod_lme4, type = "pred", transform = "plogis",
                      terms = c("froh_all [all]", "age_cent[-1.4, 2.6, 5.6]", "lamb")) %>% 
        as_tibble() %>% 
        filter(facet == 0)
df2 <- get_model_data(mod_lme4, type = "pred", transform = "plogis",
                      terms = c("froh_all [all]",  "age_cent[-2.4, -1.4]", "lamb")) %>% 
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


# across age classes
library(broom.mixed)
ageclass_inbreeding <- function(ageclass, sex1) {
        df <- annual_survival %>% mutate(sex = as.character(sex)) %>% filter(age == ageclass)# %>% 
                #filter(sex == sex1)
        mod <- glmer(survival ~ froh_all10_cent + twin + sex + (1|sheep_year) + (1|birth_year) + (1|mum_id) , # + (1|birth_year) + (1|mum_id)
                          family = binomial, data = df,
                          control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        tidy(mod, conf.int = TRUE)
}


grid_df <- expand_grid(age = 0:10, sex = c("M", "F"))
all_mods <- pmap(grid_df, ageclass_inbreeding)

b_age <- map_df(all_mods, function(x) x[x$term == "froh_all10_cent", c("term", "estimate", "conf.low", "conf.high")]) %>% 
                cbind(grid_df)

ggplot(b_age, aes(age, estimate)) +
        geom_point() +
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) #+
       # facet_wrap(~sex)
        





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

# load mcmc
mod_AS <- readRDS("output/AS_mod_MCMCglmm_logit_08m_slice")


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


# model checking

mod_inla <- readRDS("output/AS_mod_oar.rds")
mod_inla$summary.fixed
library(brinla)

ppp <- vector(mode = "numeric", length = nrow(annual_survival))
for (i in (1:nrow(annual_survival))) {
        ppp[i] <- inla.pmarginal(q = annual_survival$survival[i], marginal = mod_inla$marginals.fitted.values[[i]])
}

hist(ppp)
