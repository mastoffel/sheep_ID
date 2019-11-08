library(tidyverse)

# check chains
roh_chain <- read_delim("output/bglr_chr/chr_26ETA_roh_parBayesC.dat", " ", col_names = FALSE)
plot(roh_chain$X1, type = "line")

load_effs <- function(chr) {
        dat <- read_delim(paste0("output/bglr_chr/marker_effects_chr_", chr, ".txt"), " ")
        dat %>% mutate(chr = chr)
}

save_load <- safely(load_effs)

all_effs <- map(2, save_load) %>% 
                map("result") %>% 
                bind_rows() %>% 
                group_by(chr) %>% 
                mutate(b_roh_stand = b_roh * n(),
                       n = n())

ggplot(all_effs, aes(1:nrow(all_effs), b_roh^2)) + geom_point()


