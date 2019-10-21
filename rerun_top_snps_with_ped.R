library(pedigreemm)
library(tidyverse)
library(INLA)
library(lme4qtl)
library(furrr)
# data
gwas_df <- read_delim("output/early_survival_top_snps_pca.txt", delim = " ") %>% 
                mutate(id = as.factor(id),
                       id2 = id,
                       sex = as.factor(sex),
                       twin = as.factor(twin),
                      # survival = as.factor(survival),
                       birth_year = as.factor(birth_year),
                       sheep_year = as.factor(sheep_year))
load("data/sheep_ped.RData")

snps <- read_delim("output/top_snps_gwas_pca.txt", delim = " ")
#~~~~~~~~~INLA~~~~~~~~~~~~~~~~#

# inla prep
sheep_ped_inla <- sheep_ped %>% 
        as_tibble() %>% 
        rename(id = ID,
               mother = MOTHER,
               father = FATHER) %>% 
        mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "F", "888")) %>% 
        mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "M", "999")) %>% 
        #filter(!is.na(id)) %>% 
        mutate(father = ifelse(is.na(father), 0, father)) %>% 
        mutate(mother = ifelse(is.na(mother), 0, mother)) %>% 
        mutate_if(is.character, list(as.numeric)) %>% 
        as.data.frame() 

comp_inv <- AnimalINLA::compute.Ainverse(sheep_ped_inla)
ainv <- comp_inv$Ainverse
ainv_map <- comp_inv$map
Cmatrix <- sparseMatrix(i=ainv[,1],j=ainv[,2],x=ainv[,3])

add_index_inla <- function(dat) {
        Ndata <- dim(dat)[1]
        dat$IndexA <- rep(0, times = Ndata)
        for(i in 1:Ndata) dat$IndexA[i] <- which(ainv_map[,1]==dat$id[i])
        dat
}
gwas_df <- add_index_inla(gwas_df)
gwas_df <- gwas_df %>% mutate(IndexA2 = IndexA)

# model
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))

# get snps
top_snps <- gwas_df %>% dplyr::select(contains("roh")) %>% names(.) %>% str_replace(., "roh_", "")

# inla model
run_mod <- function(snp) {
        formula_surv <- as.formula(paste('survival ~ 1 + sex + age_std + age2_std + twin', snp, paste0("roh_", snp),
                                         'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                         'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                         'f(IndexA2, model = "iid", hyper = list(prec = prec_prior))',
                                         'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))
        
        mod_inla <- inla(formula=formula_surv, family="binomial",
                         data=gwas_df, control.compute = list(dic = TRUE))
        
}

# set up plan
plan(multiprocess, workers = 8)
# run in parallel
all_top_snp_fits <- future_map(top_snps, run_mod)

saveRDS(all_top_snp_fits, file = "output/inla_fits_top_snps.rds")






# 
# 
# 
# 
# snps_top <- read_delim("output/top_snps_gwas_pca.txt", delim = " ")
# snps_top[snps_top$snp.name == "oar3_OAR9_4828463", ]
# 
# 
# 
# 
# 
# 
# 
# #~~~~~~~~~pedigreemm~~~~~~~~~~~~~~~~#
# 
# # time saver function for modeling
# nlopt <- function(par, fn, lower, upper, control) {
#         .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
#                                   opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
#                                               maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
#         list(par = res$solution,
#              fval = res$objective,
#              conv = if (res$status > 0) 0 else res$status,
#              message = res$message
#         )
# }
# 
# pedS4 <- pedigree(sire=as.character(sheep_ped$FATHER), dam=as.character(sheep_ped$MOTHER), label=as.character(sheep_ped$ID))  
# 
# # make second ID
# gwas_df$id2 <- gwas_df$id
# 
# gwas_df %>% 
#         sample_frac(0.1) -> mod_df
# 
# start_time <- Sys.time()
# mod <- pedigreemm(survival ~ sex + twin + age_std + age2_std + roh_oar3_OAR15_63750918 + oar3_OAR15_63750918 + (1|id) + (1|id2) + (1|birth_year) + (1|sheep_year),
#                   pedigree=list(id=pedS4), data=mod_df,
#                   family = "binomial", 
#                   control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
# end_time <- Sys.time()
# end_time - start_time
# summary(mod)
# 
# 
# 
# # try asreml
# library(asreml)
# library(MasterBayes)
# ##### Prepare pedigree for asreml #####
# # pedigree format: Needs to be sorted, and individual, father, mother
# sheep_ped <- read_delim("../sheep/data/SNP_chip/20190711_Soay_Pedigree.txt", delim = "\t")[c(1,3,2)] %>% 
#         filter(ID %in% gwas_df$id) %>% 
#         as.data.frame() %>% 
#         orderPed()
# # rename(Calf = ID,
# #        Sire = FATHER,
# #        Dam = MOTHER)
# # inverse relationship mat
# sheep_ainv <- asreml::ainverse(sheep_ped)
# gwas_df_mod <- gwas_df %>% sample_frac(1) 
# 
# mod <- asreml(fixed = survival ~ 1 + sex + twin + age_std + age2_std + roh_oar3_OAR9_4828463 + oar3_OAR9_4828463,
#               random = ~vm(id, sheep_ainv) + idv(id2) + idv(birth_year) + idv(sheep_year), 
#               data = gwas_df_mod, na.action = na.method(x=c("omit")),
#               family = asr_binomial(), maxit = 1000) # "omit" "include"
# summary(mod)
# 
# mod_sum <- summary(mod, coef = TRUE)
# mod_sum$varcomp
# mod_sum$coef.fixed
# 
# snps[snps$snp.name == "oar3_OAR9_4828463", ]
# 
