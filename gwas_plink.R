
# gwas
# potentially gridlmm ,saige

library(tidyverse)
library(snpStats)
library(data.table)
# calculate 
# --indep-pairwise 1000 200 0.5
system(paste0("~/programs/plink --bfile ../sheep/data/SNP_chip/sheep_geno_imputed_04092019 --sheep --out data/gwas/sheep_chr20 ",
              "--mind 0.1 --chr 20 --recode --make-bed --indep-pairwise 1000 500 0.7 ")) 

#9077 SNPs
# recode SNPs as hom/het
# cat sheep_chr20.ped | awk '{printf $1" "$2" "$3" "$4" "$5" "$6; for(i=7;i<6+9077*2;i=i+2) { j=i+1; if($i==$j) {printf " C C";} else {printf " C T";}} print""}' > sheep_chr20_het.ped


# fitness data
load("model_in/fitness_roh_df.RData")
load("model_in/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing_imputation.txt", delim = " ")

survival <- fitness_data %>% 
        filter(Age == 0) %>% 
        filter(!is.na(Survival)) %>% 
        filter(!is.na(FROH_all )) %>% 
        dplyr::select(ID, Survival, SEX, BIRTHYEAR) %>% 
        mutate(Survival = Survival + 1)

table(survival$Survival)

write_delim(tibble(FID= survival$ID, FAMILY = survival$ID, Survival = survival$Survival), "data/gwas/survival", col_names = TRUE)
write_delim(tibble(ids = survival$ID, ids2 = survival$ID), "data/gwas/survival_ids", col_names = FALSE)
# plink data
ped_data <- read_delim("data/gwas/sheep_chr20.ped", delim = " ", col_names = FALSE) 
ped_data_mod <- ped_data %>% rename(ID = X1) %>% 
                        left_join(survival[c(1,2)], by = "ID")
ped_data_mod2 <- ped_data_mod %>% 
                        mutate(X6 = Survival) %>% 
                        mutate(X6 = ifelse(is.na(X6), -9, X6)) %>% 
                        dplyr::select(-Survival)

write_delim(ped_data_mod2, "data/gwas/sheep_chr20_mod.ped", delim = " ", col_names = FALSE)
                     

system(paste0("~/programs/plink --ped data/gwas/sheep_chr20_mod.ped --map data/gwas/sheep_chr20.map --keep data/gwas/survival_ids ", # --covar quantitative.covars --linear 
              "--adjust --assoc --allow-no-sex --out data/gwas/inbreeding_analysis ",
              "--sheep --indep-pairwise 1000 500 0.8 ")) 

