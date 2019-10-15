# # run on server
library(lme4)
library(tidyverse)
library(broom.mixed)
#source("theme_clean.R")
library(snpStats)
library(data.table)
library(furrr)

#chr <-20

# SNP data
# which chromosome
#chr <- 25

# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# top snps from gwas
top_snps <- read_delim("output/top_snps_gwas.txt", delim = " ")

# roh data
file_path <- "data/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path)

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019_pruned"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

snps_map_sub <- full_sample$map %>% filter(chromosome %in% unique(top_snps$chromosome))

# survival data
early_survival <- fitness_data %>% 
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
        mutate_at(c("id", "birth_year", "sex", "sheep_year"), as.factor) %>% 
        mutate(age2 = age^2) %>% 
        mutate(age_std = as.numeric(scale(age)),
               age2_std = as.numeric(scale(age2))) %>% 
        as.data.frame() 

# prepare additive genotypes subset
snps_sub <- full_sample$map %>%
        filter(chromosome %in% unique(top_snps$chromosome)) %>%
        .$snp.name
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_sub], Class = "numeric"),
                      rownames = "id")

# subset roh 
roh_sub <- roh_lengths %>% filter(CHR %in% unique(top_snps$chromosome))

# define vectorized seq to work with mutate
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

# create indices for all rohs
roh_snps <- roh_sub %>% 
        as_tibble() %>% 
        # sample_frac(0.01) %>% 
        mutate(index1 = as.numeric(match(SNP1, names(geno_sub))),
               index2 = as.numeric(index1 + NSNP - 1)) %>% 
        mutate(all_snps = seq2(from = index1, to = index2)) %>% 
        group_by(IID) %>% 
        summarise(all_snps = list(all_snps)) %>% 
        mutate(all_snps = simplify_all(all_snps)) %>% 
        mutate(IID = as.character(IID)) %>% 
        rename(id = IID)
# join roh_snps
roh_snps_reord <- geno_sub %>% 
        dplyr::select(id) %>% 
        left_join(roh_snps, by = "id") %>% 
        dplyr::rename(id_roh = id)
# prepare roh yes/no matrix
roh_mat <- matrix(data = 0, nrow = nrow(geno_sub), ncol = ncol(geno_sub))
# set 1 where SNP is in an roh
roh_list <- pmap(roh_snps_reord, function(id_roh, all_snps) {
        df <- as.matrix(t(as.numeric(c(id_roh, rep(0, ncol(geno_sub) - 1)))))
        df[, all_snps] <- 1
        df
}) 

# make a tibble like for genotypes but with 0/1 for whether a SNP is in ROH or not
roh_df <- do.call(rbind, roh_list) %>% 
        as_tibble() %>% 
        mutate(V1 = as.character(V1))
names(roh_df) <- c("id", paste0("roh_", names(geno_sub)[-1]))


# make some space
rm(full_sample)
rm(roh_list)
rm(roh_mat)

# join additive and roh data to survival for gwas
early_survival_gwas <- early_survival %>% 
        dplyr::select(id, survival, sex, twin, birth_year, mum_id, sheep_year, age, age2) %>% 
        #left_join(geno_sub, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        as_tibble()

# take fewer loci
survival_svm_full <- early_survival_gwas %>% 
                        select(id:age2, paste0("roh_", top_snps$snp.name)) %>% 
                        mutate(age = as.numeric(scale(age)),
                               age2 = as.numeric(scale(age2))) %>% 
                        select(-mum_id) %>% 
                        mutate(survival = ifelse(survival == 1, "yes", "no")) %>% 
                        mutate(survival = factor(survival))

# do some machine learning
library(caret)
survival_svm <- survival_svm_full %>% 
                        select(-sheep_year) %>% 
                        select(-birth_year) %>% 
                        sample_frac(1) 
                        
set.seed(3456)
trainIndex <- createDataPartition(survival_svm$survival, p = .8, 
                                  list = FALSE, 
                                  times = 1)

survival_train <- survival_svm[trainIndex, ]
survival_test <- survival_svm[-trainIndex, ]

## 10-fold CV
# fitControl <- trainControl(
#         method = "repeatedcv",
#         number = 10,
#         ## repeated ten times
#         repeats = 10)


# Control params for SVM
ctrl <- trainControl(
        method = "cv", 
        number = 5, 
        classProbs = TRUE,                 
        summaryFunction = twoClassSummary  # also needed for AUC/ROC
)
# Tune an SVM
set.seed(5628)  # for reproducibility
surv_svm_auc <- train(
        survival ~ ., 
        data = survival_train, 
        method = "svmRadial",               
        preProcess = c("center", "scale"),  
        metric = "ROC",  # area under ROC curve (AUC)       
        trControl = ctrl,
        tuneLength = 10,
        na.action = "na.omit"
)
surv_svm_auc
trellis.par.set(caretTheme())
plot(surv_svm_auc)

predict(surv_svm_auc, survival_test, type = "prob")
head(surv_svm_auc$pred)
