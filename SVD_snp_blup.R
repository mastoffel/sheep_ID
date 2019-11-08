library(tidyverse)
library("corpcor")
library(snpStats)
library(tidyimpute)

chr <- 20

# general data
# data
load("data/fitness_roh_df.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# genotypes
# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019_pruned"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map_sub <- full_sample$map %>%  filter(chromosome == chr) 
# prepare additive genotypes subset
snps_sub <- full_sample$map %>% filter(chromosome == chr) %>%  .$snp.name
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_sub], Class = "numeric"), rownames = "id")

# fitness
# survival data
annual_survival <- fitness_data %>% 
        dplyr::rename(birth_year = BIRTHYEAR,
                      sheep_year = SheepYear,
                      age = Age,
                      id = ID,
                      twin = TWIN,
                      sex = SEX,
                      mum_id = MOTHER,
                      froh_short = FROH_short,
                      froh_medium = FROH_medium,
                      froh_long = FROH_long,
                      froh_all = FROH_all,
                      froh_not_roh = hom,
                      survival = Survival) %>% 
        # some individuals arent imputed well and should be discarded 
        filter(!(id %in% IDs_lots_missing$id)) %>% 
        #filter(age == 0) %>% 
        filter(!is.na(survival)) %>% 
        filter(!is.na(froh_all)) %>% 
        filter(!(is.na(birth_year) | is.na(sheep_year))) %>%  # no mum_id here
        mutate_at(c("id", "birth_year", "sex", "sheep_year", "survival"), as.factor) %>% 
        mutate(age2 = age^2) %>% 
        mutate(age_std = as.numeric(scale(age)),
               age2_std = as.numeric(scale(age2))) %>% 
        as.data.frame() 


# testset
annual_survival_sub <- annual_survival %>% 
                         sample_frac(0.05) 
# genos
geno_mat <- annual_survival_sub %>% 
                select(id) %>% 
                left_join(geno_sub) %>% 
                select(-id) %>% 
                map_df(function(x) {
                        if(length(x) <= 1) return(NULL) 
                        out <- (x-mean(x, na.rm = TRUE))/ sd(x, na.rm = TRUE)
                        }) %>% 
                impute_mean() %>% 
                as.matrix()

svd_geno <- svd(geno_mat)
S_mat <- diag(length(svd_geno$d)) * svd_geno$d
