# # run on server
library(tidyverse)
library(snpStats)
library(data.table)
library(BGLR)

#~~~~~~~~~~~~~~~~~~~~~ GWAS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`#
# run gwas
input_folder <- "data/"
# prop_genome
#prop_geno <- 0.01
# useful to have _ at the end as some other stuff gets attached
run_name <- "first_run_svd_"
# with / at end
output_folder <- "/exports/eddie/scratch/mstoffel/bglr/"
#output_folder <- "output/bglr/"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# non genetic variables
annual_survival_gwas <- fread("data/annual_survival_gwas_vars.txt") %>% as.data.frame()
# response
y <- as.numeric(annual_survival_gwas$survival)
# ids belonging to genetic matrices
ids <- fread("data/roh_ids.txt") %>% as_tibble()

# 
US_roh <- fread(paste0(input_folder, "roh_US.txt")) %>% 
        mutate(id = ids$id) %>% 
        right_join(data.frame(id = annual_survival_gwas[, "id"])) %>% 
        dplyr::select(-id) %>% 
        as.matrix() %>% 
        unname()
#US_roh <- cbind(US_roh, matrix(0, nrow = nrow(US_roh), ncol = abs(diff(dim(US_roh)))))

US_add <- fread(paste0(input_folder, "snps_US.txt")) %>%  
        mutate(id = ids$id) %>% 
        right_join(data.frame(id = annual_survival_gwas[, "id"])) %>% 
        dplyr::select(-id) %>% 
        as.matrix() %>% 
        unname()

#2# Setting the linear predictor
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

#3# Fitting the model
fm <- BGLR2(y=y,ETA=ETA, nIter=100000, burnIn=20000, thin = 80, 
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
estimates <- tibble( b_roh = fm$ETA$roh$b, sd_b_roh = fm$ETA$roh$SD.b, b_add = fm$ETA$add$b, sd_b_add = fm$ETA$add$SD.b)

saveRDS(list(model_overview, eta_non_gen), file=paste0(output_folder, run_name, "mod.rds"))
write_delim(estimates, path = paste0(output_folder,"estimates_", run_name, "mod.txt"), delim = " ")







