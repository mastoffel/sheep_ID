# # run on server
library(MCMCglmm)
library(tidyverse)
load("model_in/lrt_roh_df.RData")
load("model_in/sheep_ped.RData")

# prep
lrt_roh_df <- lrt_roh_df %>% 
        # rename(animal = ID) %>% 
        mutate_at(c("DeathYear", "BirthYear", "MumID"), as.factor)

inv.phylo <- MCMCglmm::inverseA(sheep_ped)
sheep_ped_inv <- solve(inv.phylo$Ainv)
rownames(sheep_ped_inv) <- rownames(inv.phylo$Ainv)

lrt_roh_df_mod1 <- lrt_roh_df %>% 
        filter(Sex == 1) %>% 
        filter(!is.na(MumID)) %>% 
        filter(!is.na(BirthYear)) %>% 
        filter(!is.na(DeathYear)) %>% 
        #sample_frac(0.4) %>% 
        as.data.frame() %>% 
        mutate(animal = as.factor(animal)) 

# model0 MCMCglmm
sheep_ped_MCMCglmm <- sheep_ped %>% 
        rename(animal = ID)

phenvar <- var(lrt_roh_df_mod1$n_offspring)
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=phenvar/5, nu=0.002),
                                                G2 = list(V=phenvar/5, nu=0.002),
                                                G3 = list(V=phenvar/5, nu=0.002),
                                                G4 = list(V=phenvar/5, nu=0.002)))

start_time <- Sys.time()

nitt   <- 1100000
burnin <- 100000
thin <- 1000

mod1_chain1 <- MCMCglmm(n_offspring ~ 1 + FROH, random = ~animal + MumID + BirthYear + DeathYear, 
                            prior = prior, pedigree = sheep_ped_MCMCglmm, 
                            data = lrt_roh_df_mod1, family = "poisson",
                            nitt = nitt, burnin = burnin, thin = thin)

# mod1_chain2 <- MCMCglmm(n_offspring ~ 1 + FROH, random = ~animal + MumID + BirthYear + DeathYear, 
#                        prior = prior, pedigree = sheep_ped_MCMCglmm, 
#                        data = lrt_roh_df_mod1, family = "poisson",
#                        nitt = nitt, burnin = burnin, thin = thin)

end_time <- Sys.time()
mod1_time <- end_time - start_time

mod1 <- list(mod1_chain1, mod1_time)
save(mod1, file = "model_out/lrt_mcmcglmm_poisson.RData")

