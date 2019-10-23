library(tidyverse)
library(lme4)
library(performance)

# get df
survival_top_snps <- read_delim("output/early_survival_top_snps_pca.txt", delim = " ") %>% 
                                mutate(age_std = as.numeric(scale(age)),
                                       age2_std = as.numeric(scale(age2))) %>% 
                                mutate(sex = as.factor(sex), twin = as.factor(twin),
                                       birth_year = as.factor(birth_year),
                                       sheep_year = as.factor(sheep_year),
                                       id = as.factor(id))

# how much variation do all ROH snps explain?
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

# get roh snps
roh_snps <- names(early_survival_top_snps)[str_detect(names(early_survival_top_snps), "roh")][-1]

# check
test <- survival_top_snps %>% 
        filter(!is.na(survival_top_snps$survival)) %>% 
        summarise_at(vars(contains("roh")), function(x) length(table(x)))

# with roh_snps
glme_form <- reformulate(c("sex", "twin", "age_std", "age2_std", roh_snps, "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
mod1 <- glmer(glme_form, data = survival_top_snps, family = "binomial",
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod1)
r2_nakagawa(mod1)

glme_form2 <- reformulate(c("sex", "twin", "age_std", "age2_std",  "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
mod2 <- glmer(glme_form2, data = survival_top_snps, family = "binomial",
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod2)
library(performance)
r2_nakagawa(mod2)


# get all inla mods
top_snps_gwas <- read_delim("output/top_snps_gwas_pca.txt", delim = " ")
mods <- list.files("output/inla_mods_gwas_top/", full.names = TRUE)
all_inla <- map(mods, read_rds)
all_eff <- map_df(all_inla, function(mod) {
        mod <- as_tibble(mod$summary.fixed, rownames = "pred")
        mod[c(7), ]
        }) %>% 
        mutate(pred = str_replace(pred, "roh_", ""))

top_snps_gwas %>% 
        select(snp.name, estimate) %>% 
        rename(pred=snp.name) %>% 
        inner_join(all_eff, by = "pred") %>% 
        ggplot(aes(estimate, mean)) + geom_point() + geom_smooth(method = "lm")









# split 
library(caret)
library(lme4)
roh_snps <- paste0("roh_", top_snps$snp.name)
trainIndex <- createDataPartition(early_survival_top_snps$survival, p = .7, 
                                  list = FALSE, 
                                  times = 1)

survival_train <- early_survival_top_snps[trainIndex, ]
survival_test <- early_survival_top_snps[-trainIndex, ]
survival_test %>% filter((birth_year%in%survival_train$birth_year)) %>% 
        filter((sheep_year %in% survival_train$sheep_year))


lme_form <- reformulate(c(roh_snps, "sex", "twin", "age_std", "age2_std",  "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
mod1 <- glmer(lme_form, data = survival_train, family = "binomial")
summary(mod1)
?predict

newdat <- survival_test %>% rename(survival_org = survival) %>% 
        mutate(survival = predict(mod1, newdat, 
                                  allow.new.levels = TRUE, type = "response")) %>% 
        mutate(surv = ifelse(survival > 0.5, 1, 0)) %>% 
        mutate(surv = as.factor(surv),
               survival_org = as.factor(survival_org))

confusionMatrix(data = newdat$surv, reference = newdat$survival_org)

library(merTools)
predictInterval(mod1, newdata = survival_test,
                level = 0.95, n.sims = 5)

mm <- model.matrix(terms(mod1),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod1),mm))










library(INLA)
# inla
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

add_index_inla <- function(traits) {
        Ndata <- dim(traits)[1]
        traits$IndexA <- rep(0, times = Ndata)
        for(i in 1:Ndata) {
                traits$IndexA[i] <- which(ainv_map[,1]==traits$id[i])
        }
        traits
}

early_survival_gwas <- add_index_inla(early_survival_gwas)
# for all other random effects too
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))

formula_surv <- as.formula(paste('survival ~ oar3_OAR15_63743364 + roh_oar3_OAR15_63743364 + sex + age_std + age2_std + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(id, model = "iid", hyper = list(prec = prec_prior))',
                                 # 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))

mod_inla <- inla(formula=formula_surv, family="binomial",
                 data=early_survival_gwas, 
                 control.compute = list(dic = TRUE))
summary(mod_inla)

formula_surv_norel <- as.formula(paste('survival ~ oar3_OAR15_63743364 + roh_oar3_OAR15_63743364 + sex + age_std + age2_std + twin + 1', 
                                       'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                       'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                       'f(id, model = "iid", hyper = list(prec = prec_prior))', sep = " + "))
# 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
#'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))

mod_inla2 <- inla(formula=formula_surv_norel, family="binomial",
                  data=early_survival_gwas, 
                  control.compute = list(dic = TRUE))
summary(mod_inla2)
