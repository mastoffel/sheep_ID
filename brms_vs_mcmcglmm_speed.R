# run on server
library(brms)
library(dplyr)
library(MCMCglmm)
library(tidyverse)
library(asreml)
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

#save(lrt_roh_df, file = "model_in/lrt_roh_df.RData")
#save(sheep_ped_inv, file = "model_in/sheep_ped_inv.RData")

lrt_roh_df_mod1 <- lrt_roh_df %>% 
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

start_time <- Sys.time()
model0_brms <- brm(
        n_offspring ~ 1 + FROH + (1|animal) + (1|MumID) + (1|BirthYear), data = lrt_roh_df_mod1 , 
        family = gaussian(), cov_ranef = list(animal = sheep_ped_inv),
        chains = 1, cores = 1, iter = 5000, warmup = 1000, thin = 10
)
end_time <- Sys.time()
model0_brms_time <- end_time - start_time


# model0 MCMCglmm
sheep_ped_MCMCglmm <- sheep_ped %>% 
                        rename(animal = ID)

phenvar <- var(lrt_roh_df_mod1$n_offspring)
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=phenvar/5, nu=0.002),
                                                G2 = list(V=phenvar/5, nu=0.002),
                                                G3 = list(V=phenvar/5, nu=0.002),
                                                G4 = list(V=phenvar/5, nu=0.002)))

start_time <- Sys.time()
model0_MCMCglmm <- MCMCglmm(n_offspring ~ 1 + FROH, random = ~animal + MumID + BirthYear + DeathYear, family = "gaussian",
                  prior = prior, pedigree = sheep_ped_MCMCglmm, 
                  data = lrt_roh_df_mod1, 
                  nitt = 110000, burnin = 10000, thin = 100)
end_time <- Sys.time()
model0_MCMCglmm_time <- end_time - start_time
save(model0_MCMCglmm, file = "model_out/lrt_mcmcglmm_gaussian.RData")


summary(model0_MCMCglmm)
summary(model0_brms)


start_time <- Sys.time()
model1_MCMCglmm <- MCMCglmm(n_offspring ~ 1 + FROH, random = ~animal + MumID + BirthYear + DeathYear, 
                            prior = prior, pedigree = sheep_ped_MCMCglmm, 
                            data = lrt_roh_df_mod1, family = "poisson",
                            nitt = 1100000, burnin = 100000, thin = 100)
end_time <- Sys.time()
model1_MCMCglmm_time <- end_time - start_time
save(model1_MCMCglmm, file = "model_out/lrt_mcmcglmm_poisson.RData")
summary(model1_MCMCglmm)

plot(model1_MCMCglmm$Sol)
plot(model1_MCMCglmm$VCV)

library(broom.mixed)
tidy(model1_MCMCglmm, conf.int = TRUE)

# asreml

library(asreml)
# inverse relationship mat
sheep_ainv <- asreml::ainverse(sheep_ped)

# 1 + FROH, random = ~vm(sheep_ped) + idv(MumID) + idv(BirthYear) + idv(DeathYear)
mod <- asreml(fixed = n_offspring ~ 1 + FROH, random = ~vm(animal, sheep_ainv) + idv(MumID) + idv(BirthYear), 
              data = lrt_roh_df_mod1,
             # family = asr_negative.binomial(link = "identity", phi = 1),
              family = list(asr_gaussian()),
              na.action = na.method(x=c("omit")))

mod_sum <- summary(mod, coef = TRUE)
mod_sum$varcomp
mod_sum$coef.fixed
