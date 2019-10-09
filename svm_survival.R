# # run on server
library(lme4)
library(tidyverse)
library(broom.mixed)
#source("theme_clean.R")
library(snpStats)
library(data.table)
library(furrr)

chr <-20

# SNP data
# which chromosome
#chr <- 25

# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

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

# prepare additive genotypes subset
snps_sub <- full_sample$map %>% 
        filter(chromosome == chr) %>% 
        .$snp.name
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_sub], Class = "numeric"),
                      rownames = "id")

# subset roh 
roh_sub <- roh_lengths %>% filter(CHR == chr)

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
        dplyr::select(id, survival, sex, twin, birth_year, mum_id) %>% 
        #left_join(geno_sub, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        as_tibble()

# take fewer loci
survival_svm <- early_survival_gwas[c(2:4, sort(sample((7:ncol(early_survival_gwas)), 500)))] %>% 
                        mutate(survival = factor(survival, labels = c("no", "yes"))) %>% 
                        sample_frac(0.1) 


# do some machine learning
library(caret)
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
        number = 10, 
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
