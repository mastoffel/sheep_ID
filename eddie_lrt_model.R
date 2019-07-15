# run on server
library(brms)
library(dplyr)
library(MCMCglmm)
library(tidyverse)
library(INLA)
library(AnimalINLA)

load("model_in/lrt_roh_df.RData")
load("model_in/sheep_ped.RData")

inv.phylo <- MCMCglmm::inverseA(sheep_ped)
sheep_ped_inv <- solve(inv.phylo$Ainv)
rownames(sheep_ped_inv) <- rownames(inv.phylo$Ainv)

# males
lrt_roh_df_mod1 <- lrt_roh_df %>% 
        filter(Sex == 1) %>% 
        filter(!is.na(MumID)) %>% 
        filter(!is.na(BirthYear)) %>% 
        filter(!is.na(DeathYear)) 


ggplot(lrt_roh_df, aes(FROH, n_offspring)) +
        geom_smooth(method = "lm") +
        geom_point(size = 0.1, alpha = 0.3) +
        theme_classic() +
        facet_wrap(Sex~., scales = "free_y") +
        scale_y_sqrt()


# prep
lrt_roh_df <- lrt_roh_df %>% 
        rename(animal = ID)

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
        sample_frac(0.05) %>% 
        as.data.frame()

start_time <- Sys.time()
model0_brms <- brm(
        n_offspring ~ 1 + FROH + (1|animal), data = lrt_roh_df_mod1 , 
        family = negbinomial(), cov_ranef = list(animal = sheep_ped_inv),
        chains = 1, cores = 1, iter = 2000, warmup = 1000, thin = 10
)
end_time <- Sys.time()
model0_brms_time <- end_time - start_time

summary(model0)
plot(model0)

# model0 MCMCglmm
sheep_ped_MCMCglmm <- sheep_ped %>% 
                        rename(animal = ID)

phenvar <- var(lrt_roh_df_mod1$n_offspring)
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=phenvar/2, nu=0.002)))

start_time <- Sys.time()
model0_MCMCglmm <- MCMCglmm(n_offspring ~ 1 + FROH, random = ~animal, family = "gaussian",
                  prior = prior, pedigree = sheep_ped_MCMCglmm, 
                  data = lrt_roh_df_mod1, 
                  nitt = 2000, burnin = 1000, thin = 10)
end_time <- Sys.time()
model0_MCMCglmm_time <- end_time - start_time







plot(model0_MCMCglmm$Sol)
plot(model0_MCMCglmm$VCV)
summary(model0_MCMCglmm)


model1 <- brm(
        n_offspring ~ 1 + FROH + (1|MumID) + (1|BirthYear) + (1|DeathYear) + (1|animal), data = lrt_roh_df_mod1 , 
        family = negbinomial(), cov_ranef = list(animal = sheep_ped_inv),
        chains = 2, cores = 2
        #control = list(adapt_delta = 0.9)
)

model2 <- brm(
        n_offspring ~ 1 + FROH + (1|MumID) + (1|BirthYear) + (1|DeathYear) + (1|animal), data = lrt_roh_df_mod1 , 
        family = poisson(), cov_ranef = list(animal = A),
        chains = 2, cores = 2
        #control = list(adapt_delta = 0.9)
)

model3 <- brm(
        n_offspring ~ 1 + FROH + (1|MumID) + (1|BirthYear) + (1|DeathYear) + (1|animal), data = lrt_roh_df_mod1 , 
        family = hurdle_negbinomial(), cov_ranef = list(animal = A),
        chains = 2, cores = 2
        #control = list(adapt_delta = 0.9)
)






summary(model2)
plot(model2)
coef(model2)

pp_check(model2)


# inverse relationship mat
sheep_ainv <- asreml::ainverse(sheep_ped)

#attr(sheep_ainv, "rowNames")
# 1 + FROH, random = ~vm(sheep_ped) + idv(MumID) + idv(BirthYear) + idv(DeathYear)
mod <- asreml(fixed = n_offspring ~ 1 + FROH, random = ~vm(sheep_ainv) + idv(MumID), 
              data = lrt_roh_df_mod,
              #family = asr_negative.binomial(link = "log", dispersion = 1, phi = NA),
              #   family = list(asr_poisson()),
              na.action = na.method(x=c("omit")))

mod_sum <- summary(mod, coef = TRUE)
mod_sum$varcomp
mod_sum$coef.fixed
