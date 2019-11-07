# pca
library(snpStats)
library(tidyverse)

# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

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


# make PCs for all annual survival individuals

# extraxt all ids for gwas
ids <- as.character(unique(annual_survival$id))

# make data.frame as input for plink 
id_df <- data.frame(1, ids) 
write.table(id_df, "data/ids_annual_survival.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# pca using only individuals in annual survival analysis
system(paste0("~/programs/plink --bfile data/sheep_geno_imputed_ram_27092019_pruned --sheep ",
              "--pca 40 --keep data/ids_annual_survival.txt --nonfounders --out output/sheep_pca_ann_survival")) # --keep data/ind_testset_80.txt --keep data/ids_annual_survival.txt

# check scree plot / looks like 7 is a good number
pca_eigenval <- read_lines("output/sheep_pca_ann_survival.eigenval") %>% as.numeric() %>% plot()

# prepare eigenvectors to include in gwas
col_from <- paste0("X", 3:42)
col_to <- paste0("pc", 1:40)
pca_eigenvec <- read_delim("output/sheep_pca_ann_survival.eigenvec", delim = " ", col_names = FALSE) %>% 
                        dplyr::rename(id = X2) %>% 
                        select(-X1) %>% 
                        dplyr::rename_at(vars(col_from), ~col_to) %>% 
                        select(id, pc1, pc2, pc3, pc4, pc5, pc6, pc7)

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
