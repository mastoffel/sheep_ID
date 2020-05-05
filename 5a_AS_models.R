# Run survival models 
library(lme4)
library(tidyverse)
library(broom.mixed)
source("theme_simple.R")
library(INLA)
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
                 control.compute = list(dic = TRUE, config=TRUE),
                 control.inla = list(correct = TRUE)
                 )

saveRDS(mod_inla, file = "output/AS_mod_INLA_398k.rds")
mod_inla <- readRDS("output/AS_mod_INLA_398k.rds")
mod_inla$summary.fixed

# POISSON model for inbreeding load --------------------------------------------
# to calculate lethal equivalents according to Nietlisbach et al 2019

formula_surv <- as.formula(paste('survival ~ froh_all10_cent + age_std + age_std2 + sex + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 #'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla_pois3 <- inla(formula=formula_surv, family="poisson",
                       data=annual_survival, 
                       control.compute = list(dic = TRUE, config=TRUE)
)

saveRDS(mod_inla_pois3, file = "output/AS_mod_INLA_398k_poisson_full.rds")
summary(mod_inla_pois3)

# inbreeding load (lethal equivalents) 2B
#  froh_all10_cent  4.78 [2.76, 6.81]  
mod_inla_pois3$summary.fixed %>% 
        rownames_to_column(var = "pred") %>% 
        as_tibble() %>% 
        select(c(1,2,4,6)) %>% 
        setNames(c("pred", "est", "lower_ci", "upper_ci")) %>% 
        filter(pred == "froh_all10_cent") %>% 
        mutate(across(is.numeric, function(x) abs(2 * x/0.10)))

# plot INLA marginal effects ---------------------------------------------------
mod_inla <- readRDS("output/AS_mod_INLA_398k.rds")
fun <- function(...) {
        one <-  invlogit(Intercept + 
                                 df1$x1 * froh_all10_cent + 
                                 df1$x2 * age_cent + 
                                 df1$x3 * lamb1 + 
                                 df1$x4 * twin1 + 
                                 df1$x5 * sexM + 
                                 df1$x6 * `froh_all10_cent:lamb1` + 
                                 df1$x7 * `froh_all10_cent:age_cent`)
        
        return (list(one))
}

invlogit <- function(x) exp(x)/(1+exp(x))

froh <- seq(from = min(annual_survival$froh_all10_cent), to = (max(annual_survival$froh_all10_cent)), by = 0.1)
age <- c(-2.4, -1.4, 1.6, 4.6)
combined_df <- expand_grid(froh, age) %>% 
        mutate(lamb = ifelse(age == -2.4, 1, 0),
               twin = 0.15,
               sex = 0.4,
               frohxlamb = froh*lamb,
               frohxage = froh*age) 
names(combined_df) <- paste0("x", 1:7)

xx <- inla.posterior.sample(1000, mod_inla)
marg_means <- purrr::map(1:nrow(combined_df), function(x) {
        df1 <<- combined_df[x, ]
        out <- inla.posterior.sample.eval(fun, xx)
}) 


d <- marg_means %>% 
        map(as_tibble) %>% 
        bind_rows() %>% 
        pmap_df(function(...) {
                samp <- as.numeric(unlist(list(...)))
                c(mean = mean(samp), quantile(samp, probs = c(0.025, 0.975)))
        }) %>% 
        bind_cols(combined_df) %>% 
        .[, 1:5] %>% 
        setNames(c("prediction", "ci_lower", "ci_upper", "froh", "age")) %>% 
        mutate(age = as.factor(round(age + mean(annual_survival$age), 0)),
               froh = (froh + mean(annual_survival$froh_all10))/10)
        

#saveRDS(d, file = "output/AS_mod_INLA_predictions_for_plot.rds")
ggplot(d, aes(froh, prediction)) +
        geom_line(aes(color = age), size = 1.5) +
        geom_ribbon(aes(x=froh, ymin = ci_lower, ymax = ci_upper, fill = age, color = age),
                    alpha = 0.2, linetype = 2, size = 0.1)+
        scale_color_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        scale_fill_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) 




# INLA models across age classes

# ~~~~~~~~~~~~~~~~~~~~~ Survival across different life-stages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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
saveRDS(all_inla_mods, file = "output/inla_survival_models_diff_ages_froh_all_398k.rds")

# plot for supplementary
inla_plot <- map_df(all_inla_mods, function(x) as_tibble(x$summary.fixed[2, ]), .id = "age") %>% 
        mutate(age = as.numeric(age) - 1) %>% 
        filter(age < 10) %>% 
        mutate(age = fct_inorder(as.factor(age))) 
   

ggplot(inla_plot, aes(age, mean)) +
        geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`),  width = 0.5) +
        geom_point(size = 3, shape = 21, col = "#4c566a", fill = "#eceff4", # "grey69"
                   alpha = 1, stroke = 0.7) +
        xlab("Age class") + 
        ylab(expression(beta[F[ROH]])) +
        geom_hline(yintercept = 0, linetype='dashed', colour =  "#4c566a", size = 0.3) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) -> p_froh_across_life

ggsave("figs/sup_froh_across_life.jpg", width = 5, height = 3)




# Modeling with lme4 and MCMCglmm ----------------------------------------------


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

start_time <- Sys.time()
mod_lme4 <- glmer(survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
                  family = binomial, data = annual_survival,
                  control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
end_time <- Sys.time()
end_time - start_time
summary(mod_lme4)

df1 <- get_model_data(mod_lme4, type = "eff", 
                      terms = c("froh_all10_cent [all]", "age_cent[-1.4, 2.6, 5.6]", "lamb")) %>% 
        as_tibble() %>% 
        filter(facet == 0)
df2 <- get_model_data(mod_lme4, type = "eff", 
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




