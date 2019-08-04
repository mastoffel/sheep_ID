library(broom.mixed)
library(MCMCglmm)
load("model_out/lrt_mcmcglmm_poisson.RData")

mod <- mod1[[1]]
plot(mod$Sol)
plot(mod$VCV)

tidy(mod, conf.int = TRUE)
autocorr(mod$VCV)

summary(mod)
HPDinterval(mod$VCV)
