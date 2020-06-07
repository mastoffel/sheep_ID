# transform plinks PCA output to data files for GWAS
library(snpStats)
library(tidyverse)

# data
load("data/survival_mods_data.RData")
load("data/sheep_ped.RData")
#IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# survival data
annual_survival <- fitness_data %>% 
        # filter na rows
        filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
        mutate(age_cent = age - mean(age, na.rm = TRUE),
               age_cent2 = age_cent^2,
               age_std = as.numeric(scale(age)),
               age_std2 = age_std^2,
               # times 10 to estimate a 10% percent increase
               froh_all10 = froh_all * 10,
               froh_all10_cent = froh_all10 - mean(froh_all10, na.rm = TRUE),
               lamb = ifelse(age == 0, 1, 0),
               lamb_cent = lamb - mean(lamb, na.rm = TRUE),
               lamb = as.factor(lamb)) %>% 
        as.data.frame() 


# pcs made in script 1_ROH_calling.R
# check scree plot / looks like 7 is a good number
pca_eigenval <- read_lines("output/sheep_pca_oar.eigenval") %>% as.numeric() %>% plot()

# prepare eigenvectors to include in gwas
col_to <- paste0("pc", 1:20)
pca_eigenvec <- read_delim("output/sheep_pca_oar.eigenvec", delim = " ", col_names = FALSE) %>% 
                        dplyr::rename(id = X2) %>% 
                        dplyr::select(-X1) %>% 
                        #rename(across(starts_with("X")), ~col_to)
                        dplyr::rename_at(vars(X3:X22), ~col_to) %>% 
                        dplyr::select(id, pc1, pc2, pc3, pc4, pc5, pc6, pc7)

# plot
ggplot(pca_eigenvec, aes(pc1, pc7)) + geom_point()

# save
pca_eigenvec %>% 
        write_delim("data/ann_surv_pca.txt", " ")







