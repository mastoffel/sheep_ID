
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
        as.data.frame()# %>% 
      #  sample_frac(0.05)
        
emp_ratio <- sum(annual_survival$survival)/nrow(annual_survival)

# INLA -------------------------------------------------------------------------
# posterior predictive check
mod_inla <- readRDS("output/AS_mod_oar.rds")

set.seed(144)
xx <- inla.posterior.sample(1000, mod_inla)

mod_inla$.args$control.family$control.link$model
contents <- mod_inla$misc$configs$contents

sample1 <- xx[[1]]$latent
inds <- str_detect(rownames(sample1), "sheep_year")
sum(inds)
1/var(sample1[inds])
xx[[1]]$hyperpar


1/var(xx[[1]]$latent)

get_ratio <- function(xx_sample) {
  pred_ind <- str_detect(rownames(xx_sample$latent), pattern = "Predictor")
  pred <- xx_sample$latent[pred]
  #sum(resp)/length(resp)
}

sim_ratios <- map_df(xx, get_ratio)

ggplot(as_tibble(sim_ratios), aes(value)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = emp_ratio)

# which names available
mod_inla$misc$configs$contents 
fun <- function(df1) {
        resp <- plogis(Intercept + 
                df1$froh_all10_cent * froh_all10_cent + 
                df1$age_cent * age_cent + 
                model.matrix(~lamb, data = df1)[, 2] * lamb1 + 
                model.matrix(~twin, data = df1)[, 2] * twin1 + 
                model.matrix(~sex, data = df1)[, 2] * sexM + 
                df1$froh_all10_cent * model.matrix(~lamb, data = df1)[, 2] * `froh_all10_cent:lamb1` + 
                df1$froh_all10_cent * df1$age_cent * `froh_all10_cent:age_cent`)
        return(resp)
}

df1 <- annual_survival
mod_inla$misc$configs$contents
out <- inla.posterior.sample.eval(fun, xx, df1 = df1)
out_bin <- round(out)
emp_ratio <- sum(annual_survival$survival)/nrow(annual_survival)
sim_ratio <- as_tibble(colSums(out)/nrow(out))

ggplot(sim_ratio, aes(value)) +
        geom_histogram(bins = 30) +
        geom_vline(xintercept = emp_ratio)


marg_means <- purrr::map(1:nrow(combined_df), function(x) {
        df1 <<- combined_df[x, ]
        out <- inla.posterior.sample.eval(fun, xx)
}) 
marg_means <- purrr::map(1:nrow(combined_df), function(x) {
        df1 <<- combined_df[x, ]
        out <- inla.posterior.sample.eval(fun, xx)
}) 



str(sims[[1]], max.level = 1)
test <- sims[[1]]$logdens
str(as_tibble(test))
xpost <- generate(mod_inla, annual_survival, survival ~ ., n.samples = 5, n = 10)
?predict.inla

preds <- predict(mod_inla, formula = ~froh_all10_cent)


# try rstanarm
library(rstanarm)
library(bayesplot)
mod_stan <- stan_glmer(survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb_cent + sex + twin + 
                              (1|birth_year) + (1|sheep_year) + (1|id), data = annual_survival,
                        family = binomial)
summary(mod_stan)
posterior <- as.matrix(mod_stan)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars = c("froh_all10_cent"),
           prob = 0.8) + plot_title

library("dplyr")
color_scheme_set("brightblue")
mod_stan %>%
        posterior_predict(draws = 500) %>%
        ppc_stat(y = annual_survival$survival,
                #         group = mtcars$carb,
                         stat = "median")


# posterior predictive distr. with lme4

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

mod1 <- glmer(survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin + + (1|sheep_year) + (1|birth_year) + (1|id),
              family = binomial(link = 'logit'), data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

emp_ratio <- sum(annual_survival$survival) / nrow(annual_survival)
new_resps <- simulate(mod1, nsim = 1000, re.form = NULL)
sim_ratio <- as_tibble(colSums(new_resps)/nrow(new_resps))

# posterior predictive check
p1 <- ggplot(sim_ratio, aes(value)) +
  geom_histogram(bins = 40, fill = "grey") +
  geom_vline(xintercept = emp_ratio, colour = "cornflowerblue") +
  xlab("Proportion of survival") +
  ggtitle("Distribution of simulated survival (grey) and empirical
survival (blue line)") +
  theme_minimal()

# prediction
predict_accuracy <- function(rep) {
  test <- sample(1:nrow(annual_survival), round(0.2 * nrow(annual_survival)))
  train <- c(1:nrow(annual_survival))[-test]
  df_train <- annual_survival[train, ]
  df_test <- annual_survival[test, ]
  mod <- glmer(survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin + (1|sheep_year) + (1|birth_year) + (1|id) + (froh_all10_cent|age),
               family = binomial(link = 'logit'), data = df_train,
               control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
  
  preds <- predict(mod, newdata = df_test, allow.new.levels = TRUE, re.form = NULL,
                   type = "response")
  preds <- ifelse(preds > 0.5, 1, 0)
  cm <- as.matrix(table(actual = annual_survival[test, "survival"], predicted = preds))
  accuracy <- sum(diag(cm)) /  sum(cm)
  accuracy
  
}

library(future)
plan(multiprocess, workers = 3)
all_acc <- map_dbl(1:100, predict_accuracy)
all_acc
p1 <- ggplot(as_tibble(all_acc), aes(value*100)) +
  geom_histogram(bins = 14, fill = "grey") +
  theme_minimal() +
  xlab("Survival prediction accuracy %")
p1
ggsave("figs/pred_accuracy.jpg", p1, height = 2, width = 3.4)
