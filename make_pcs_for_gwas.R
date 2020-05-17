# pca
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









# make PCs for 80% of annual survival individuals
# extraxt all ids for gwas
set.seed(7336)
ids <- unique(as.character(annual_survival$id)) %>% 
        sample(length(.) * 0.8, replace = FALSE)
# make data.frame as input for plink 
id_df <- data.frame(1, ids) 
write.table(id_df, "data/ids_annual_survival_testset80.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# pca using only individuals in annual survival analysis
system(paste0("plink --bfile data/sheep_geno_imputed_ram_27092019_pruned --sheep ",
              "--pca 40 --keep data/ids_annual_survival_testset80.txt --nonfounders --out output/sheep_pca_ann_survival_testset80")) # --keep data/ind_testset_80.txt --keep data/ids_annual_survival.txt

# check scree plot / looks like 7 is a good number
pca_eigenval <- read_lines("output/sheep_pca_ann_survival_testset80.eigenval") %>% as.numeric() %>% plot()

# prepare eigenvectors to include in gwas
col_from <- paste0("X", 3:42)
col_to <- paste0("pc", 1:40)
pca_eigenvec <- read_delim("output/sheep_pca_ann_survival_testset80.eigenvec", delim = " ", col_names = FALSE) %>% 
        dplyr::rename(id = X2) %>% 
        select(-X1) %>% 
        dplyr::rename_at(vars(col_from), ~col_to) %>% 
        select(id, pc1, pc2, pc3, pc4, pc5, pc6, pc7)

# plot
ggplot(pca_eigenvec, aes(pc1, pc2)) + geom_point()

# save
pca_eigenvec %>% 
        write_delim("data/ann_surv_pca_testset80.txt", " ")
