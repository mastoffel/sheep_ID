library(tidyverse)
library(BGLR)
# check chains
roh_chain <- read_delim("output/bglr_chr/chr_20ETA_roh_parBayesC.dat", " ", col_names = FALSE)
plot(roh_chain$X2, type = "line")

fixed <- BGLR::readBinMat("output/bglr_chr/chr_1ETA_id_b.bin")

load_effs <- function(chr) {
        dat <- read_delim(paste0("output/bglr_chr/marker_effects_chr_", chr, ".txt"), " ")
        dat %>% mutate(chr = chr)
}

save_load <- safely(load_effs)

all_effs <- map(20, save_load) %>% 
                map("result") %>% 
                bind_rows() %>% 
                group_by(chr) %>% 
                #filter(chr != 4) %>% 
                mutate(b_roh_stand = b_roh/mean(b_roh))

ggplot(all_effs, aes(1:nrow(all_effs), b_roh^2, color = as.factor(chr))) + geom_point()

map(1:26, save_load) %>% 
        map("result") %>% 
        bind_rows() %>% 
        group_by(chr) %>% 
        summarise(mean_eff = median(b_roh), n_eff = n()) %>% 
        ggplot(aes(mean_eff, n_eff)) + geom_point()

