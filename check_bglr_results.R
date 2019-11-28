library(BGLR)
library(tidyverse)
library(lobstr)
library(coda)

chain_roh<-read_delim('output/bglr/first_run_svd_ETA_roh_parBayesC.dat', " ", col_names = FALSE)  %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_roh, aes(iter, X2)) + geom_line()

chain_roh<-read_delim('output/bglr/var_sel/var_sel_svd_ETA_add_parBayesC.dat', " ", col_names = FALSE)  %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_roh, aes(iter, X2)) + geom_line()

hist(chain_roh$X2)
chain_birth_year<-read_delim('output/bglr/first_run_svd_ETA_birth_year_varB.dat', " ", col_names = FALSE)  %>% 
        mutate(iter = 1:nrow(.))
ggplot(chain_birth_year, aes(iter, X1)) + geom_line()
hist(chain_birth_year$X1)

# get beta
beta <- readBinMat("output/bglr/first_run_svd_ETA_roh_b.bin") %>% 
        colMeans()
b <- read_delim('output/bglr/estimates_third_run_svd_mod.txt', " ") %>% .$b_roh 

# bglr chr wise
effs <- list.files("output/bglr_chr",pattern = "*.txt", full.names = TRUE)
all_effs <- map_dfr(effs, read_delim, " ", .id = "chr") %>% 
                mutate(chr = as.numeric(chr)) %>% 
                mutate(chr2 = case_when(
                        chr == 1 ~ 10,
                        chr == 2 ~ 11,
                        TRUE ~ chr + 1
                )) %>% 
                select(-chr) %>% 
                arrange(chr2) %>% 
                filter(chr2 == 10)
ggplot(all_effs, aes(1:nrow(all_effs), b_roh^2)) + geom_point()

plot(beta, b)

markers <- read_delim('output/bglr/var_sel/marker_effects.txt', " ")  %>% 
        mutate(iter = 1:nrow(.))
ggplot(markers, aes(iter, effs^2)) + geom_point()


markers <- read_delim('output/bglr/marker_effects_bglr_third.txt', " ")  %>% 
        mutate(iter = 1:nrow(.))
ggplot(markers, aes(iter, effs)) + geom_point()


markers <-read_delim('output/bglr/var_sel/marker_effects_var_sel.txt', " ", col_names = TRUE)  %>% 
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
