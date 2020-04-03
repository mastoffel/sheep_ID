library(MCMCglmm)
library(MasterBayes)
library(tidyverse)
# data
load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")
ped <- orderPed(sheep_ped)

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

# MCMCglmm ---------------------------------------------------------------------
prior1 <-list(R=list(V=1, fix=1), 
              G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), 
                     G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                     G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000),
                     G4=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)
              ))
annual_survival$animal <- annual_survival$id

# ~ 8 minutes for 5000/10/1000
# mod_AS <- MCMCglmm(survival~froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin, 
#                    random=~birth_year + sheep_year + id + animal,
#                    data=annual_survival,
#                    family="threshold", 
#                    #family="categorical",
#                    trunc=TRUE,
#                    prior=prior1, 
#                    pr=TRUE,
#                    nitt=200000,thin=80,burnin=20000,
#                    #nitt=5000,thin=10,burnin=1000,
#                    pedigree=ped)
# 
# saveRDS(mod_AS, file = "output/AS_mod_MCMCglmm")

mod_AS2 <- MCMCglmm(survival~froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin, 
                    random=~birth_year + sheep_year + id + animal,
                    data=annual_survival,
                    #family="threshold", 
                    family="categorical",
                    #trunc=TRUE,
                    slice = TRUE,
                    prior=prior1, 
                    pr=TRUE,
                    nitt=200000,thin=80,burnin=20000,
                    #nitt=5000,thin=10,burnin=1000,
                    pedigree=ped)

saveRDS(mod_AS2, file = "output/AS_mod_MCMCglmm_cat_slice")