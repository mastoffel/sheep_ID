library(BGLR)
library(tidyverse)
library(lobstr)
library(coda)

first_run_mod[[2]]$birth_year

chain_roh<-read_delim('output/bglr/first_run_svd_ETA_roh_parBayesC.dat', " ", col_names = FALSE)  %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_roh, aes(iter, X2)) + geom_line()
chain_roh_old <- read_delim('output/bglr/first/first_run_svd_ETA_roh_parBayesC.dat', " ", col_names = FALSE)  %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_roh_old, aes(iter, X2)) + geom_line()

markers <-read_delim('output/bglr/var_sel/marker_effects_bglr_2nd.txt', " ", col_names = TRUE)  %>% 
        mutate(iter = 1:nrow(.))
ggplot(markers, aes(iter, effs^2)) + geom_point()

chain <- read_delim('output/bglr/var_sel/first_run_svd_ETA_roh_parBayesC.dat', " ", col_names = FALSE)  %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_roh_old, aes(iter, X1)) + geom_line()

mod <- readRDS("output/bglr/first_run_svd_mod.rds")

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
