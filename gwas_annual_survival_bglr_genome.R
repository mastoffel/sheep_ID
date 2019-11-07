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
output_folder <- "/exports/eddie/scratch/mstoffel/bglr/"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# snps
snp_map <- read_delim(paste0(input_folder, "snp_map.txt"), ' ')

# non genetic variables
annual_survival_gwas <- fread("data/annual_survival_gwas_vars.txt")
# response
y <- annual_survival_gwas$survival


# snps for gwas
snps_sub <- snp_map %>% 
        #sample_frac(prop_geno) %>% 
        .$snp.name

#2# Setting the linear predictor
ETA<-list(fixed = list(~factor(sex)+factor(twin)+age_std+age2_std,
                        data=annual_survival_gwas,model='FIXED', 
                        saveEffects=TRUE),
             id = list(~factor(id), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
     sheep_year = list(~factor(sheep_year), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
     birth_year = list(~factor(birth_year), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
            roh = list(X_roh=fread(paste0(input_folder, "annual_survival_gwas_roh.txt"),  
                   select = paste0("roh_", snps_sub)) %>% as.matrix(), model='BayesC', 
                   saveEffects=TRUE),
            add = list(X_add=fread(paste0(input_folder, "annual_survival_gwas_snps.txt"), 
                  select = snps_sub) %>% as.matrix(), model = 'BayesC',
                  saveEffects=TRUE) # 
)

#3# Fitting the model
fm <- BGLR2(y=y,ETA=ETA, nIter=5000, burnIn=2000, thin = 50, 
        response_type = "ordinal",
        saveEnv=TRUE,
        # additional iterations with the following two lines
        #BGLR_ENV = paste0(output_folder, run_name, "BGLR_ENV.RData"), # default NULL
        #newChain = FALSE, # default TRUE
        # where to save
        saveAt = paste0(output_folder, run_name)) 


# save parts of the output
model_overview <- fm[-length(fm)]
eta_non_gen <- fm$ETA[1:4]
marker_effects <- tibble(snp_roh = names(fm$ETA$roh$b), b_roh = fm$ETA$roh$b, sd_b_roh = fm$ETA$roh$SD.b,
                 snp_add = names(fm$ETA$add$b), b_add = fm$ETA$add$b, sd_b_add = fm$ETA$add$SD.b)

saveRDS(list(model_overview, eta_non_gen), file=paste0(output_folder, run_name, "mod.rds"))
write_delim(marker_effects, path = paste0(output_folder,"marker_effects_", run_name, "mod.txt"), delim = " ")







