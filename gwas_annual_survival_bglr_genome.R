# # run on server
library(tidyverse)
library(snpStats)
library(data.table)
library(furrr)
library(caret)
library(tidyimpute)
library(BGLR)

#~~~~~~~~~~~~~~~~~~~~~ GWAS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`#
# run gwas
input_folder <- "data/"
# prop_genome
#prop_geno <- 0.01
# useful to have _ at the end as some other stuff gets attached
run_name <- "first_run_"
# with / at end
output_folder <- "output/bglr/"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# snps
snp_map <- read_delim(paste0(input_folder, "snp_map.txt"), ' ')

# non genetic variables
annual_survival_gwas <- fread("data/annual_survival_gwas_vars.txt")
# response
y <- annual_survival_gwas$survival


# snps for gwas
snps_sub <- snp_map %>% 
       # sample_frac(prop_geno) %>% 
        .$snp.name

#2# Setting the linear predictor
ETA<-list(fixed = list(~factor(sex)+factor(twin)+age_std+age2_std,
                        data=annual_survival_gwas,model='FIXED', 
                        saveEffects=TRUE),
        random = list(~factor(id) + factor(sheep_year) + factor(birth_year), 
                        data=annual_survival_gwas, model='BRR', 
                        saveEffects=TRUE),
        roh = list(X_roh=fread(paste0(input_folder, "annual_survival_gwas_roh.txt"),  
                   select = paste0("roh_", snps_sub)) %>% as.matrix(), model='BayesC', 
                   saveEffects=TRUE),
        add = list(X_add=fread(paste0(input_folder, "annual_survival_gwas_snps.txt"), 
                select = snps_sub) %>% as.matrix(), model = 'BayesC',
                saveEffects=TRUE) # 
)

#3# Fitting the model
fm <- BGLR2(y=y,ETA=ETA, nIter=1000, burnIn=1000, thin = 10, 
        response_type = "ordinal",
        saveEnv=TRUE,
        # additional iterations with the following two lines
        #BGLR_ENV = paste0(output_folder, run_name, "BGLR_ENV.RData"), # default NULL
        #newChain = FALSE, # default TRUE
        # where to save
        saveAt = paste0(output_folder, run_name)) 

saveRDS(fm, file=paste0(output_folder, run_name, "mod.rds"))


# mod <- readRDS("output/bglr/first_run_mod.rds")
# summary(mod)
# mod$ETA$add

# eff_roh <- read_delim("output/bglr/first_run_ETA_roh_parBayesC.dat", " ", col_names = F)
# eff_roh2 <- readBinMat("output/bglr/first_run_ETA_roh_b.bin")
#save_samp_roh <- readBinMat("output/bglr/first_run_ETA_roh_b.bin")





