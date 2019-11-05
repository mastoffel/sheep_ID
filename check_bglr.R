library(BGLR)
library(tidyverse)
library(lobstr)
roh <- read_delim("output/bglr/first_run_ETA_roh_parBayesC.dat", " ", col_names = FALSE)
plot(roh[, 1])

plot(scan("output/bglr/first_run_varE.dat"), type='o')

Sys.setenv('R_MAX_VSIZE'=32000000000)

mod <- readRDS("output/bglr/first_run_mod.rds")

install.packages("benchmarkme")
benchmarkme::get_ram()

str(mod)

bHat <- mod$ETA$roh$b
plot(bHat^2, ylab='Estimated Squared-Marker Effect',
        type='o',cex=.5,col=4,main='Marker Effects')
