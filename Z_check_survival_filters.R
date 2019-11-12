library(lme4)
library(tidyverse)
library(broom.mixed)
#source("theme_clean.R")
library(snpStats)
library(data.table)
library(furrr)

# best way to filter survival data

load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")
load("data/fitness_roh_df.RData")


# survival data
early_survival <- fitness_data %>% 
        dplyr::rename(birth_year = BIRTHYEAR,
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
        filter(!is.na(froh_all)) %>% 
        filter(!(is.na(birth_year) | is.na(mum_id))) %>% 
        mutate_at(c("id", "birth_year", "mum_id", "sex"), as.factor) 

# check nas
check_surv_na <- early_survival %>% filter(is.na(survival))
check_age_na <- early_survival %>% filter(is.na(age)) # none
check_surv_na_age_0 <- early_survival %>% filter(age == 0 & (is.na(survival))) %>% filter(str_detect(Comment, "aybe"))

early_survival_filter <- early_survival %>% 
        mutate(survived_1stwinter = case_when(
                # no NAs in age
                (age == 0) & (survival == 0) ~ 0,
                # if sheep was seen at any age after 0 
                (age > 0) | ((age == 0) & (survival == 1)) ~ 1,
                # deal with survival == NA
                #(age > 0) & (is.na(survival)) ~ 1
        )) %>% 
        group_by(id) %>% 
        arrange(desc(age)) %>% 
        filter(row_number() == 1)


# survival data
early_survival <- fitness_data %>% 
        dplyr::rename(birth_year = BIRTHYEAR,
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
        filter(age == 0) %>% 
        filter(!is.na(froh_all)) %>% 
        filter(!(is.na(birth_year) | is.na(mum_id))) %>% 
        mutate_at(c("id", "birth_year", "mum_id", "sex"), as.factor) %>% 
        as.data.frame() 

table(early_survival_filter$survived_1stwinter, useNA = "always")

test <- early_survival_filter %>% filter(is.na(survived_1stwinter))
        

# roh data
file_path <- "data/roh_nofilt_ram.hom"
roh_lengths <- fread(file_path)

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

snps_map_sub <- full_sample$map %>% 
        filter(chromosome == chr)  

