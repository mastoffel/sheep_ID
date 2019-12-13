library(tidyverse)
dat <- read_delim("output/bglr/BRR/marker_effects2.txt", " ")
dat_old <- read_delim("output/bglr/BRR/marker_effects_old2.txt", " ")

dat2 <- read_delim("output/bglr/BayesC_vars_sel/marker_effects2.txt", " ")
plot(dat2$effs^2)

plot(dat$effs, dat2$effs)
chain <- read_delim("output/bglr/BRR/BRR_runETA_roh_varB.dat", " ",
                    col_names = FALSE)
ggplot(chain, aes(1:nrow(chain), X1)) + geom_line()

mod <- readRDS("output/bglr/BRR/BRR_run_mod.rds")
str(mod)

cor(dat$effs, dat2$effs)

plot(abs(sample(dat2$effs, size = 1000, replace = F)))
plot(dat_old$effs^2)

par(mfrow=c(2,1))
plot(dat$effs^2)
plot(dat2$effs^2)
