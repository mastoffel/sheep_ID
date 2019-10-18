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
# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# roh data
file_path <- "data/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path) 

# survival data
early_survival <- fitness_data %>% 
        dplyr::rename(birth_year = BIRTHYEAR,
                      sheep_year = SheepYear,
                      age = Age,
                      id = ID,
                      twin = TWIN,
                      sex = SEX,
                      mum_id = MOTHER,
                      froh_short = FROH_short,
                      froh_medium = FROH_medium,
                      froh_long = FROH_long,
                      froh_all = FROH_all,
                      froh_not_roh = hom,
                      survival = Survival) %>% 
        # some individuals arent imputed well and should be discarded 
        filter(!(id %in% IDs_lots_missing$id)) %>% 
        filter(!(is.na(survival) | is.na(froh_all) | is.na(birth_year))) %>% 
        filter(!(is.na(sheep_year))) %>% 
        mutate_at(c("id", "birth_year", "sex", "sheep_year"), as.factor) %>% 
        mutate(age2 = age^2) %>% 
        mutate(age_std = as.numeric(scale(age)),
               age2_std = as.numeric(scale(age2))) %>% 
        as.data.frame() 



# inla prep
sheep_ped_inla <- sheep_ped %>% 
        as_tibble() %>% 
        rename(id = ID,
               mother = MOTHER,
               father = FATHER) %>% 
        mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "F", "888")) %>% 
        mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "M", "999")) %>% 
        #filter(!is.na(id)) %>% 
        mutate(father = ifelse(is.na(father), 0, father)) %>% 
        mutate(mother = ifelse(is.na(mother), 0, mother)) %>% 
        mutate_if(is.character, list(as.numeric)) %>% 
        as.data.frame() 

comp_inv <- AnimalINLA::compute.Ainverse(sheep_ped_inla)
ainv <- comp_inv$Ainverse
ainv_map <- comp_inv$map
Cmatrix <- sparseMatrix(i=ainv[,1],j=ainv[,2],x=ainv[,3])

add_index_inla <- function(dat) {
        Ndata <- dim(dat)[1]
        dat$IndexA <- rep(0, times = Ndata)
        for(i in 1:Ndata) dat$IndexA[i] <- which(ainv_map[,1]==dat$id[i])
        dat
}
early_survival <- add_index_inla(early_survival)
early_survival <- early_survival %>% mutate(IndexA2 = IndexA)

prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))
formula_surv <- as.formula(paste('survival ~ froh_all + sex + age_std + age2_std + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 # 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla <- inla(formula=formula_surv, family="binomial",
                 data=early_survival , 
                 control.compute = list(dic = TRUE))
summary(mod_inla)
bri.hyperpar.summary(mod_inla)
hrtbl <- bri.hyperpar.summary(mod_inla)[4, 1] / sum(sum(bri.hyperpar.summary(mod_inla)[,1], pi^2/3))

# mod2
formula_surv2 <- as.formula(paste('survival ~ froh_long + froh_medium + froh_short + sex + age_std + age2_std + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                 # 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
mod_inla2 <- inla(formula=formula_surv2, family="binomial",
                 data=early_survival , 
                 control.compute = list(dic = TRUE))

summary(mod_inla2)

saveRDS(list(mod_inla, mod_inla2), file = "output/inla_survival_models.rds")

out <- readRDS("output/inla_survival_models.rds")
inla_mod1 <- out[[1]]
inla_mod2 <- out[[2]]

summary(inla_mod1)
summary(inla_mod2)

sampvars <- 1/inla.hyperpar.sample(1000, inla_mod1)
sampicc <- sampvars[, 3]/ (rowSums(sampvars))
quantile(sampicc, c(0.025, 0.5, 0.975))

# inv logit exp(x)/(1+exp(x))
inv_logit <- function(x) exp(x)/(1+exp(x))

# plot
surv_mod_df <- inla_mod2$summary.fixed %>% 
                rbind(inla_mod1$summary.fixed[2, ]) %>% 
                as_tibble(rownames = "predictor") %>% 
                rename(lower_ci = `0.025quant`,
                       upper_ci = `0.975quant`) %>% 
                filter(predictor != "(Intercept)") %>% 
                mutate(predictor = factor(predictor, 
                                          levels = rev(c("froh_all", "froh_long", "froh_medium", "froh_short",
                                                     "sexM", "twin", "age_std", "age2_std"))))

formatC(inv_logit(seq(0, -50, by = -10)), format = "e", digits = 1)

p_fix_eff <- ggplot(surv_mod_df, aes(mean, predictor)) +
        geom_point(size = 2) +
        geom_errorbarh(aes(xmax = upper_ci, xmin = lower_ci), height = 0.5, size = 0.1) +
        geom_vline(xintercept = 0) +
        xlab("estimate") +
        theme_clean() +
        scale_x_continuous(labels = )
ggsave(filename = "figs/surv_mod_fixeff.jpg", p_fix_eff, width = 4.5, height = 3)        


# rpt and herit
sampvars <- 1/inla.hyperpar.sample(1000, inla_mod1)
sampicc <- sampvars[, 4]/ (rowSums(sampvars) + (pi^2)/3)
quantile(sampicc, c(0.025, 0.5, 0.975))




library(lme4)
mod_lme4 <- glmer(survival ~ 1 + froh_all + sex + twin + age_std + age2_std + (1|birth_year) + (1|sheep_year) + (1|id),
                  data = early_survival, family = "binomial")
summary(mod_lme4)

mod_lme42 <- glmer(survival ~ 1 + froh_long + froh_medium + froh_short + sex + twin + age_std + age2_std + (1|birth_year) + (1|sheep_year) + (1|id),
                  data = early_survival, family = "binomial")
summary(mod_lme42)
