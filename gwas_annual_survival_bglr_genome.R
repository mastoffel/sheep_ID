# # run on server
library(tidyverse)
library(snpStats)
library(data.table)
library(BGLR)

#~~~~~~~~~~~~~~~~~~~~~ GWAS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# variable selection prior probIn=2/100,counts=100 (see https://github.com/gdlc/BGLR-R/issues/32)
# BGLR with variable selection?
var_sel <- TRUE
# run name
run_name <- "first_run"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# attach to run_name when var_sel selected
if(var_sel) run_name <- paste0(run_name, "_var_sel")

# different output folder depending on whether script runs on mac or server
if (Sys.info()["sysname"] == "Darwin") output_folder <- paste0("output/bglr/", run_name)
if (Sys.info()["sysname"] == "Linux") output_folder <- paste0("/exports/eddie/scratch/mstoffel/bglr/", run_name)

# create output folder if its not there
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# input folder (should be same on mac and server)
input_folder <- "data/"

# non genetic variables
annual_survival_gwas <- fread(paste0(input_folder, "annual_survival_gwas_vars.txt")) %>% as.data.frame()
# response
y <- as.numeric(annual_survival_gwas$survival)
# ids belonging to genetic matrices
ids <- fread(paste0(input_folder, "roh_ids.txt")) %>% as_tibble()

# US from SVD 
US_roh <- fread(paste0(input_folder, "roh_US.txt")) %>% 
        mutate(id = ids$id) %>% 
        right_join(data.frame(id = annual_survival_gwas[, "id"])) %>% 
        dplyr::select(-id) %>% 
        as.matrix() %>% 
        unname()
# US from SVD for additive effects
US_add <- fread(paste0(input_folder, "snps_US.txt")) %>%  
        mutate(id = ids$id) %>% 
        right_join(data.frame(id = annual_survival_gwas[, "id"])) %>% 
        dplyr::select(-id) %>% 
        as.matrix() %>% 
        unname()

#2# Setting the linear predictor
if (!var_sel) {
        ETA<-list(fixed = list(~factor(sex)+factor(twin)+age_std+age2_std,
                               data=annual_survival_gwas,model='FIXED',
                               saveEffects=TRUE),
                  id = list(~factor(id), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
                  sheep_year = list(~factor(sheep_year), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
                  birth_year = list(~factor(birth_year), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
                  # create genetic matrices with duplicated observations
                  roh = list(X_roh=US_roh, model='BayesC', saveEffects=TRUE),
                  add = list(X_add=US_add, model='BayesC', saveEffects=TRUE)
        )
}

if (var_sel) {
        ETA<-list(fixed = list(~factor(sex)+factor(twin)+age_std+age2_std,
                               data=annual_survival_gwas,model='FIXED',
                               saveEffects=TRUE),
                  id = list(~factor(id), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
                  sheep_year = list(~factor(sheep_year), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
                  birth_year = list(~factor(birth_year), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
                  # create genetic matrices with duplicated observations
                  roh = list(X_roh=US_roh, model='BayesC', saveEffects=TRUE, probIn=2/100,counts=100),
                  add = list(X_add=US_add, model='BayesC', saveEffects=TRUE, probIn=2/100,counts=100)
        )
}


#3# Fitting the model
fm <- BGLR2(y=y,ETA=ETA, nIter=100, burnIn=10, thin = 5, 
        response_type = "ordinal",
        #saveEnv=TRUE,
        # additional iterations with the following two lines
        BGLR_ENV = "output/bglr/var_sel/var_sel_svd_BGLR_ENV.RData", # default NULL
        newChain = FALSE, # default TRUE
        # where to save
        saveAt = paste0(output_folder, "/", run_name)) 


# save parts of the output
model_overview <- fm[-length(fm)]
eta_non_gen <- fm$ETA[1:4]
estimates <- tibble( b_roh = fm$ETA$roh$b, sd_b_roh = fm$ETA$roh$SD.b, b_add = fm$ETA$add$b, sd_b_add = fm$ETA$add$SD.b)

saveRDS(list(model_overview, eta_non_gen), file=paste0(output_folder, "/", run_name, "_mod.rds"))
write_delim(estimates, path = paste0(output_folder,"/estimates_", run_name, "mod.txt"), delim = " ")


# calculate marker effects
estimates <- estimates %>% 
        .$b_roh %>% 
        as.numeric()

V <- fread("data/gen_mats/roh_V.txt") %>% as.matrix()

effs <- V%*%estimates

write_delim(tibble(effs = effs), path = paste0(output_folder, "/marker_effects.txt"), delim = " ")




