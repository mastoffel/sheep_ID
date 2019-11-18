# calculate marker effects from svd
library(tidyverse)
library(data.table)

estimates <- read_delim("/exports/eddie/scratch/mstoffel/bglr/estimates_first_run_svd_mod.txt", " ") %>% 
                .$b_roh %>% 
                as.numeric()

V <- fread("data/gen_mats/roh_V.txt") %>% as.matrix()

effs <- V%*%estimates

write_delim(tibble(effs = effs), path = "data/marker_effects_bglr.txt", delim = " ")