# run on server
library(brms)
library(dplyr)
library(MCMCglmm)
library(tidyverse)
library(asreml)
library(broom.mixed)
load("model_in/lrt_roh_df.RData")
load("model_in/sheep_ped.RData")

ggplot(lrt_roh_df, aes(FROH, n_offspring)) +
        geom_smooth(method = "lm") +
        geom_point(size = 0.1, alpha = 0.3) +
        theme_classic() +
        facet_wrap(Sex~., scales = "free_y") +
        scale_y_sqrt()

# prep
lrt_roh_df <- lrt_roh_df %>% 
       # rename(animal = ID) %>% 
        mutate_at(c("DeathYear", "BirthYear", "MumID"), as.factor)

inv.phylo <- MCMCglmm::inverseA(sheep_ped)
sheep_ped_inv <- solve(inv.phylo$Ainv)
rownames(sheep_ped_inv) <- rownames(inv.phylo$Ainv)

lrt_roh_df_mod1 <- lrt_roh_df %>%
        # males
        filter(Sex == 1) %>%  
        filter(!is.na(MumID)) %>% 
        filter(!is.na(BirthYear)) %>% 
        filter(!is.na(DeathYear)) %>% 
        mutate(MumID = as.factor(MumID),
               BirthYear = as.factor(BirthYear),
               DeathYear = as.factor(DeathYear)) %>% 
        #sample_frac(0.4) %>% 
        as.data.frame() %>% 
        mutate(animal = as.factor(animal)) %>% 
        mutate(olre = as.factor(1:nrow(.)))

# model0 MCMCglmm
sheep_ped_MCMCglmm <- sheep_ped %>% 
                        rename(animal = ID)

#phenvar <- var(lrt_roh_df_mod1$n_offspring)
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002),
                                                G2 = list(V=1, nu=0.002),
                                                G3 = list(V=1, nu=0.002),
                                                G4 = list(V=1, nu=0.002)))

nitt <- 55000
burnin <- 5000
thin <- 100
mod1 <- MCMCglmm(n_offspring ~ 1 + FROH, random = ~animal + MumID + BirthYear + DeathYear, 
                            family = "poisson",
                  prior = prior, pedigree = sheep_ped_MCMCglmm, 
                  data = lrt_roh_df_mod1, 
                  nitt = nitt, burnin = burnin, thin = thin)

plot(mod1$Sol)
plot(mod1$VCV)
summary(mod1)
tidy(mod1)

library(broom.mixed)
tidy(model1_MCMCglmm, conf.int = TRUE)

# asreml
library(asreml)

# inverse relationship mat
sheep_ainv <- asreml::ainverse(sheep_ped)

# 1 + FROH, random = ~vm(sheep_ped) + idv(MumID) + idv(BirthYear) + idv(DeathYear)
mod <- asreml(fixed = n_offspring ~ 1 + FROH, random = ~vm(animal, sheep_ainv) + 
                      idv(MumID) + idv(BirthYear)+ idv(DeathYear) + idv(olre),  #  + idv(DeathYear)
              data = lrt_roh_df_mod1,
              family = asreml::asr_negative.binomial(),
              #family = list(asr_poisson()),
              na.action = na.method(x=c("omit")))

mod_sum <- summary(mod, coef = TRUE)
mod_sum$varcomp
mod_sum$coef.fixed


