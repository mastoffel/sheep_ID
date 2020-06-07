# create example files for reviewers to run code
library(tidyverse)
# choose subset of individuals
load("data/survival_mods_data.RData") 
# sample 100 random ids for example dataset
annual_survival <- fitness_data %>% 
        # filter na rows
        filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
        .$id %>% 
        unique() %>% 
        sample(100) %>% 
        as.data.frame() %>% 
        setNames("id") %>% 
        write_delim("example_data/example_ids.txt")

# example_ids_plink
read_delim("example_data/example_ids.txt", " ") %>% 
        cbind(1, .) %>% 
        write_delim("example_data/example_ids_plink.txt", col_names = FALSE)

example_ids <- read_delim("example_data/example_ids.txt", " ")

# make example ROH
#~~ ROH for survival data subset
system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar_filt --sheep --out example_output/ROH/roh ",
              "--keep example_data/example_ids_plink.txt ",
              "--homozyg --homozyg-window-snp 50 --homozyg-snp 50 --homozyg-kb 1200 ",
              "--homozyg-gap 300 --homozyg-density 200 --homozyg-window-missing 2 ",
              "--homozyg-het 2 ",
              "--homozyg-window-het 2"))

# fitness and roh data
load("data/survival_mods_data.RData") 
fitness_data %>% 
        filter(id %in% example_ids$id) -> fitness_data
save(fitness_data, file = "example_data/survival_mods_data.RData") 

# make example plink files
system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar_filt ",
              "--sheep --keep example_data/example_ids_plink.txt ",
              "--make-bed --out example_data/sheep_geno_imputed_oar_filt"))
