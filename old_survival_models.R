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
library(performance)
library(sjPlot)
library(ggeffects)
library(patchwork)
library(effects)
library(MCMCglmm)
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


# a few lme4 trials
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
plot_model(mod2)
summary(mod2)
inv.logit(-1.197)
odds_to_prob <- function(x) exp(x) / (1+exp(x))
get_model_data(mod2, type = "est")
summary(mod2)
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
  #geom_jitter(data = annual_survival, aes(froh_all10_cent, y = -0.1), height = 0.05,
  #            alpha = 0.2, size = 3) +
  #geom_point(data = df_full, aes(x = froh_all10_cent, fit, color = age_cent)) +
  geom_line(data = df_full, aes(x = x, predicted, color = group), size = 1) +
  geom_ribbon(data= df_full, aes(x=x, ymin=conf.low, ymax=conf.high, fill = group, color = group), alpha= 0.2,
              linetype = 2, size = 0.1) +
  scale_color_viridis_d("Age", labels = c(0, 1, 4, 7)) +
  scale_fill_viridis_d("Age", labels = c(0, 1, 4, 7)) +
  theme_simple(grid_lines = FALSE, axis_lines = TRUE) 



# calc ID per ageclass
#  sex +
est_id_per_age <- function(ageclass) {
  annual_survival_sub <- annual_survival[(annual_survival$age == ageclass), ]
  mod_sub <- glmer(survival ~ froh_all10_cent + sex + twin + (1|birth_year) + (1|sheep_year),
                family = binomial, data = annual_survival_sub,
                control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
}

id_per_age <- map(0:9, est_id_per_age)
id_per_age %>% 
  map(tidy, conf.int = TRUE) %>% 
  bind_rows(.id = "age") %>% 
  filter(term == "froh_all10_cent") %>% 
  mutate(age = as.factor(as.numeric(age) - 1)) %>% 
        # sex = rep(rep(c("M", "F") ,each =9))) %>% 
  ggplot(aes(age, estimate)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.5) +
  geom_point(size = 3, shape = 21, col = "#4c566a", fill = "#eceff4", # "grey69"
             alpha = 1, stroke = 0.7) +
  xlab("Age class") + 
  ylab(expression(beta[F[ROH]])) +
  geom_hline(yintercept = 0, linetype='dashed', colour =  "#4c566a", size = 0.3) +
  theme_simple(grid_lines = FALSE, axis_lines = TRUE) -> p_froh_across_life

ggsave("figs/sup_froh_across_life.jpg", width = 5, height = 3)






mod8 <- glmer(survival ~ froh_all10_cent * age_cent + lamb + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
performance(mod8)
check_collinearity(mod8)
tidy(mod8, conf.int = TRUE)

plot_model(mod8, type = "eff", terms = c("froh_all10_cent [all]", "lamb"))

roh_eff_lamb <- ggeffect(mod8, type = "eff", terms = c("froh_all10_cent [all]", "lamb [1]")) %>% 
          as_tibble() #%>% 
         #filter(group == 1)
summary(roh_eff_lamb)
roh_eff_adult <- ggeffect(mod8, type = "eff", 
                          terms = c("froh_all10_cent [all]", "age_cent[-1.4, 3.6, 7.6]")) %>%  # -1.4, 1.6, 3.6, 5.6, 7.6
                as_tibble() 
summary(roh_eff_adult)
p_lamb <- ggplot(roh_eff_lamb, aes(x, predicted)) +
          geom_line(size = 1) +
          #geom_point(data = annual_survival %>% filter(age == 0), 
          #           aes(x = froh_all10_cent, y = survival),
          #           size = 4, alpha = 0.1) +
          geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                      alpha = 0.1, linetype = 2, color = "grey") +
          theme_simple() +
          scale_y_continuous(labels = paste0(seq(from = 0, to = 100, by = 25), "%"), limits = c(0,1)) +
          scale_x_continuous(labels = c("0.24", "0.34", "0.44"), breaks = c(0, 1, 2)) +
          ylab("Survival probability") +
          #ggtitle("Lambs") +
          xlab(expression(F[ROH])) 
p_lamb
p_adult <-  ggplot(roh_eff_adult, aes(x, predicted)) +
                geom_line(aes(color = group), size = 1) +
                geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group, color = group),
                            alpha = 0.05, size = 0.2, linetype = 2) +
                theme_simple() +
                scale_color_viridis_d("age", labels = c(1, 4, 8)) +
                scale_fill_viridis_d("age", labels = c(1, 4, 8)) +
                scale_y_continuous(labels = paste0(seq(from = 0, to = 100, by = 25), "%"), limits = c(0,1)) +
                scale_x_continuous(labels = c("0.24", "0.34", "0.44"), breaks = c(0, 1, 2)) +
                ylab("Survival probability") +
        xlab(expression(F[ROH])) 
p_adult  

p_lamb + p_adult + plot_annotation(tag_levels = "A") 

p <- plot_model(mod6)
p +
  theme_simple() +
  scale_x_continuous(limits = c(0.05, 5), breaks = c(0.05, 1, 5))
plot_model(mod7, type = "eff",  terms = c("age_cent [-1.4, 1.6, 3.6, 5.6]", "froh_all10_cent [-1, 0, 1, 2]"))
plot_model(mod8, type = "eff",  terms = c("froh_all10_cent [-1, 0, 1, 2]", "age_cent [-1.4, 2.6, 5.6]"))

plot_model(mod8, type = "eff", terms = c("age_cent [-2.4, -0.4, 1.6, 3.6]", "froh_all10_cent [-0.5, 0, 0.5, 1]"))
plot_model(mod8)
out <- ggeffect(mod8, type = "eff", terms = c("age_cent [-2.4, -0.4, 1.6, 3.6]", "froh_all10_cent [-0.5, 0, 0.5, 1]"))

plot_model(mod8, type = "eff", terms = c("froh_all10_cent", "age_cent", "lamb"))
plot_model(mod8, type = "pred", terms = c("froh_all10_cent", "lamb"))
plot_model

p2 <- plot_model(mod6, type = "eff", terms = c("age_cent [-2.4, -0.4, 1.6, 3.6, 5.6]", "froh_all10_cent  [-1, 0, 1, 2]"))
plot_model(mod6, type = "pred")

library(patchwork)
p1+p2
plot_model(mod6, type = "pred", terms = c( "froh_all10_cent [all]", "age_cent [-2.4, -1.4, -0.4, 1.4, 2.4]", "lamb_cent"))

library(effects)
effects_roh <- Effect(focal.predictors = c("froh_all10_cent", "age_cent"), mod = mod6)
summary(effects_roh)
x_roh <- as.data.frame(effects_roh) %>%
         mutate(age_cent = as.factor(age_cent))
         #mutate(froh = froh_all10_cent + mean(annual_survival$froh_all),
         #       age = age_cent + mean(annual_survival$age))
ggplot() +
        geom_point(data = annual_survival, aes(froh_all10_cent, survival)) +
        geom_point(data = x_roh, aes(x = froh_all10_cent, fit, color = age_cent)) +
        geom_line(data = x_roh, aes(x = froh_all10_cent, fit, color = age_cent)) +
        geom_ribbon(data= x_roh, aes(x=froh_all10_cent, ymin=lower, ymax=upper, fill = age_cent), alpha= 0.3) +
        theme_simple() 
        

# depression across life

run_across_life <- function(age, data) {
        dat <- data %>% filter(age == {{ age }})
        formula_age <- as.formula(paste0("survival ~ 1 + froh_all10 + sex + twin ", 
                                   "+ (1|birth_year) + (1|sheep_year) + (1|id)"))
        mod <- glmer(formula = formula_age,
                     data = dat, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod, conf.int = TRUE)
        out
}

all_mods <- purrr::map(0:8, run_across_life, annual_survival)
all_mods %>% bind_rows(.id = "age") %>% filter(term == "froh_all10") %>% 
        ggplot(aes(age, estimate)) +
        geom_point() +
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high))



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
annual_survival <- annual_survival %>% 
        mutate(IndexA2 = IndexA) 

#annual_survival_test <- annual_survival %>% sample_frac(0.2)

# model 2, with lamb 
formula_surv <- as.formula(paste('survival ~ froh_all10_cent * (lamb + age_cent) + sex + twin + 1', 
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
str(mod_inla, max.level = 1)
saveRDS(mod_inla, file = "output/inla_survival_model_full.rds")

mod_inla <- readRDS("output/inla_survival_model_full.rds")
plot(mod_inla, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=TRUE, plot.q=FALSE, plot.cpo=FALSE,
     single=FALSE)

mod_inla$summary.fixed
library(INLAutils)
autoplot(mod_inla)

trans_link_to_dat <- function(pred, mod_inla) {
        inla.rmarginal(n = 10000, marginal = mod_inla$marginals.fixed[[pred]]) %>% 
                exp() %>% 
                quantile(probs = c(0.025,0.5, 0.975))  %>% 
                t() %>% 
                as.data.frame() %>% 
                bind_cols() %>% 
                add_column(pred = pred, .before = 1)
}

fix_eff <- map_df(names(mod_inla$marginals.fixed), trans_link_to_dat, mod_inla) %>% 
           .[c(2,4,5,6,8), ] %>% 
           filter(!(pred == "(Intercept)")) %>% 
           mutate(pred = fct_inorder(as.factor(pred)))
names(fix_eff) <- c("Predictor", "lower_CI", "median", "upper_CI")

p_surv_mod <- ggplot(fix_eff, aes(median, Predictor, xmax = upper_CI, xmin = lower_CI)) +
        geom_vline(xintercept = 1, linetype='dashed', colour =  "#4c566a", size = 0.7) +
        geom_errorbarh(alpha = 1, height = 0.5,
                       size = 0.5) +
        geom_point(size = 3, shape = 21, col = "#4c566a", fill = "#eceff4", # "grey69"
                   alpha = 1, stroke = 0.7) + 
        scale_y_discrete(limits = rev(levels(fix_eff$Predictor)),
                         #labels = rev(c("FROH", "Age", "Sex(Male)", "Twin", "FROH:Age"))) +
                         labels = rev(c(expression(F[ROH]), "Age", expression(Sex[Male]), "Twin", expression(F[ROH]:Age)))) + 
                         #labels = rev(c("Froh", "Lamb", "Age", "Sex(Male)", "Twin", "Froh:lamb", "Froh:age"))) +
        scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1, 1.25)) +
        theme_simple(axis_lines = TRUE) +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.line.x = element_line(colour = "black",size = 0.8),
                axis.ticks.x = element_line(size = 0.8),
                axis.text = element_text(),
                axis.title.x = element_text(margin=margin(t=8))
        ) +
        xlab(expression(beta~and~95*"%"*~CI~(odds~of~survival)))
p_surv_mod
ggsave("figs/surv_mod_fixs.jpg", p_surv_mod, width = 3.5, height = 3)















# model 3, with age * roh interaction 
formula_surv <- as.formula(paste('survival ~ froh_all10_cent * age_cent + sex + lamb + age_cent2 + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 #'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla3 <- inla(formula=formula_surv, family="binomial",
                 data=annual_survival, 
                 control.compute = list(dic = TRUE),
                 control.inla = list(correct = TRUE))

summary(mod_inla3)

# model 4, without age_cent2
formula_surv <- as.formula(paste('survival ~ froh_all10_cent * age_cent + sex + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 #'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla4 <- inla(formula=formula_surv, family="binomial",
                  data=annual_survival, 
                  control.compute = list(dic = TRUE),
                  control.inla = list(correct = TRUE))

summary(mod_inla4)



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


