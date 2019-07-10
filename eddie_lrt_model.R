# run on server
library(brms)
library(dplyr)
library(MCMCglmm)
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

model1 <- brm(
        n_offspring ~ 1 + FROH + (1|MumID) + (1|BirthYear) + (1|DeathYear) + (1|animal), data = lrt_roh_df_mod1 , 
        family = negbinomial(), cov_ranef = list(animal = sheep_ped_inv),
        chains = 2, cores = 2, warmup = 100, iter = 200, thin = 1
        #control = list(adapt_delta = 0.9)
)
save(model1, file = "model_out/lrt_negbinomial")
# model2 <- brm(
#         n_offspring ~ 1 + FROH + (1|MumID) + (1|BirthYear) + (1|DeathYear) + (1|animal), data = lrt_roh_df_mod1 , 
#         family = poisson(), cov_ranef = list(animal = A),
#         chains = 2, cores = 2
#         #control = list(adapt_delta = 0.9)
# )
# 
# model3 <- brm(
#         n_offspring ~ 1 + FROH + (1|MumID) + (1|BirthYear) + (1|DeathYear) + (1|animal), data = lrt_roh_df_mod1 , 
#         family = hurdle_negbinomial(), cov_ranef = list(animal = A),
#         chains = 2, cores = 2
#         #control = list(adapt_delta = 0.9)
# )