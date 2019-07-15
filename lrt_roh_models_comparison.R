# # run on server
library(MCMCglmm)
library(tidyverse)
library(asreml)
library(nadiv)
load("model_in/lrt_roh_df.RData")
load("model_in/sheep_ped.RData")

# ped
# inv.phylo <- MCMCglmm::inverseA(sheep_ped)
# sheep_ped_inv <- solve(inv.phylo$Ainv)
# rownames(sheep_ped_inv) <- rownames(inv.phylo$Ainv)

lrt_roh_df <- lrt_roh_df %>% 
        mutate_at(c("animal", "Sex", "MumID", "BirthYear", "DeathYear"), as.factor) %>% 
        #sample_frac(0.4) %>% 
        as.data.frame() 

# males
lrt_roh_df_males <- lrt_roh_df %>% 
        filter(Sex == 1) %>%
        filter(!is.na(MumID)) %>%
        filter(!is.na(BirthYear)) %>%
        filter(!is.na(DeathYear))


################################ asreml ########################################

# inverse for asreml
sheep_ainv <- asreml::ainverse(sheep_ped)

# olre
lrt_roh_df_males_asreml <- 
        lrt_roh_df_males %>% 
        mutate(olre = factor(1:nrow(.))) %>% 
        as.data.frame()

#attr(sheep_ainv, "rowNames")
# 1 + FROH, random = ~vm(sheep_ped) + idv(MumID) + idv(BirthYear) + idv(DeathYear)
asreml.options(maxit = 30)
mod1 <- asreml(fixed = n_offspring ~ 1 + FROH, random = ~vm(animal, sheep_ainv) + idv(MumID) + idv(BirthYear) + idv(DeathYear), 
              data = lrt_roh_df_males_asreml,
              #family = asr_negative.binomial(link = "log", dispersion = 1, phi = NA),
              family = asr_poisson(),
              na.action = na.method(x=c("omit")))

mod_sum <- summary(mod1, coef = TRUE)
mod_sum$varcomp
mod_sum$coef.fixed

mod2 <- asreml(fixed = n_offspring ~ 1 + FROH, random = ~vm(animal, sheep_ainv) + idv(MumID) + idv(BirthYear) + idv(DeathYear) + idv(olre), 
               data = lrt_roh_df_males_asreml,
               #family = asr_negative.binomial(link = "log", dispersion = 1, phi = NA),
               family = asr_poisson(),
               na.action = na.method(x=c("omit")))

mod_sum <- summary(mod2, coef = TRUE)
mod_sum$varcomp
mod_sum$coef.fixed

lrt_roh_df_males_asreml2 <- lrt_roh_df_males %>% 
        mutate(log(n_offspring + 0.00001))
mod3 <- asreml(fixed = n_offspring ~ 1 + FROH, random = ~vm(animal, sheep_ainv) + idv(MumID) + idv(BirthYear) + idv(DeathYear), 
               data = lrt_roh_df_males_asreml2 ,
               #family = asr_negative.binomial(link = "log", dispersion = 1, phi = NA),
               family = asr_gaussian(),
               na.action = na.method(x=c("omit")))

mod_sum <- summary(mod3, coef = TRUE)
mod_sum$varcomp
mod_sum$coef.fixed





############################# MCMCglmm #########################################
lrt_roh_df_males
# model0 MCMCglmm
sheep_ped_MCMCglmm <- sheep_ped %>% 
        rename(animal = ID)

#phenvar <- var(lrt_roh_df_mod1$n_offspring)
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002),
                                                G2 = list(V=1, nu=0.002),
                                                G3 = list(V=1, nu=0.002),
                                                G4 = list(V=1, nu=0.002)))

prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002)))
nitt <- 10000
burnin <- 2000
thin <- 10
mod1 <- MCMCglmm(n_offspring ~ 1 + FROH, random = ~animal, #+ MumID + BirthYear + DeathYear, 
                 family = "poisson",
                 prior = prior, pedigree = sheep_ped_MCMCglmm, 
                 data = lrt_roh_df_males, 
                 nitt = nitt, burnin = burnin, thin = thin)

plot(mod1$Sol)
plot(mod1$VCV)
summary(mod1)
tidy(mod1)





