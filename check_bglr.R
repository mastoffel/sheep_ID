library(BGLR)
library(tidyverse)
library(lobstr)
library(coda)

first_run_mod[[2]]$birth_year

chain_fix<-read_delim('output/bglr/first_run_ETA_fixed_b.dat', " ") %>% 
                mutate(iter = 1:nrow(.))
ggplot(chain_fix, aes(iter, `factor(sex)M`)) + geom_line()

chain_mu<-read_delim('output/bglr/first_run_mu.dat', " ", col_names = FALSE) %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_mu, aes(iter, X1)) + geom_line()

chain_thresh<-read_delim('output/bglr/first_run_thresholds.dat', " ", col_names = FALSE) %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_mu, aes(iter, X1)) + geom_line()

chain_roh<-read_delim('output/bglr/first_run_ETA_roh_parBayesC.dat', " ", col_names = FALSE) %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_roh, aes(iter, X2)) + geom_line()

chain_add<-read_delim('output/bglr/first_run_ETA_add_parBayesC.dat', " ", col_names = FALSE) %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_add, aes(iter, X1)) + geom_line()

chain_random<-read_delim('output/bglr/first_run_ETA_random_varB.dat', " ", col_names = FALSE) %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_roh, aes(iter, X1)) + geom_line()

chain_random<-read_delim('output/bglr/first_run_ETA_birth_year_varB.dat', " ", col_names = FALSE) %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_random, aes(iter, X1)) + geom_line()


roh <- read_delim("output/bglr/first_run_ETA_roh_parBayesC.dat", " ", col_names = FALSE)
plot(roh[, 2])

plot(scan("output/bglr/first_run_varE.dat"), type='o')

Sys.setenv('R_MAX_VSIZE'=32000000000)

mod <- readRDS("output/bglr/first_run_mod.rds")

install.packages("benchmarkme")
benchmarkme::get_ram()

str(mod)

marker_effs <- read_delim("output/bglr/marker_effects_third_run_mod.txt", " ") %>% 
                mutate(num_snp = 1:nrow(.))
ggplot(marker_effs, aes(num_snp, b_roh^2)) + geom_point()

bHat <- mod$ETA$roh$b
plot(bHat^2, ylab='Estimated Squared-Marker Effect',
        type='o',cex=.5,col=4,main='Marker Effects')
