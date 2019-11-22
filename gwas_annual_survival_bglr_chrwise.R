# # run on server
library(tidyverse)
#source("theme_clean.R")
library(snpStats)
library(data.table)
library(furrr)
#library(caret)
#library(tidyimpute)
library(BGLR)

# for running on server, chromosome input through array job
chr_inp  <- commandArgs(trailingOnly=TRUE)
if (!(length(chr_inp) == 0)) {
        chr <- as.numeric(chr_inp[[1]])
} else {
        # SNP data
        # which chromosome
        chr <- 20
}

#~~~~~~~~~~~~~~~~~~~~~ GWAS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`#
# run gwas
snp_map <- read_delim("data/snp_map.txt", ' ')
annual_survival_gwas <- fread("data/annual_survival_gwas_vars.txt")
y <- annual_survival_gwas$survival

# with / at end
output_folder <- paste0("/exports/eddie/scratch/mstoffel/bglr_chr/chr_", chr, "/")
# output_folder <- ("output/bglr_chr/testruns/")
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

run_gwas_per_chr <- function(chr) {
        snps_sub <- snp_map %>% filter(chromosome == chr) %>% .$snp.name
        #2# Setting the linear predictor
        ETA<-list( fixed = list(~factor(sex)+factor(twin)+age_std+age2_std,
                data=annual_survival_gwas,model='FIXED', saveEffects=TRUE),
                #random = list(~factor(id) + factor(sheep_year) + factor(birth_year), data=annual_survival_gwas, model='BRR'),
                id = list(~factor(id), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
                sheep_year = list(~factor(sheep_year), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
                birth_year = list(~factor(birth_year), data=annual_survival_gwas, model='BRR', saveEffects=TRUE),
                roh = list(X=fread("data/annual_survival_gwas_roh.txt",  select = paste0("roh_", snps_sub)) %>% as.matrix(), model='BayesC', probIn=1/100,counts=100), # , saveEffects=TRUE
                add = list(X=fread("data/annual_survival_gwas_snps.txt", select = snps_sub) %>% as.matrix(), model = 'BayesC', probIn=1/100,counts=100) # , saveEffects=TRUE / BayesC
        )
        
        #3# Fitting the model / previous 10K / 5k burning
        fm <- BGLR2(y=y,ETA=ETA, nIter=10000, burnIn=1000, thin = 10, response_type = "ordinal",
                #saveEnv=TRUE,
                #newChain = TRUE,
                #BGLR_ENV = paste0(output_folder,"chr_", chr, "BGLR_ENV.RData"),
                saveAt = paste0(output_folder, "chr_", chr))
        
        marker_effects <- tibble(snp_roh = names(fm$ETA$roh$b), b_roh = fm$ETA$roh$b, sd_b_roh = fm$ETA$roh$SD.b,
                                 snp_add = names(fm$ETA$add$b), b_add = fm$ETA$add$b, sd_b_add = fm$ETA$add$SD.b)
        write_delim(marker_effects, path = paste0(output_folder,"marker_effects_", "chr_", chr, ".txt"), delim = " ")
        
        #save(fm, file=paste0(output_folder, "chr_", chr, ".rda"))
        #rm(ETA)
}

walk(chr, run_gwas_per_chr)
#out <- walk(20:26, run_gwas_per_chr)
#out <- walk(10:19, run_gwas_per_chr)
#out <- walk(1:9, run_gwas_per_chr)
# 
# load("output/bglr/chr_24.rda")
# out <- read_delim("output/bglr/chr24_ETA_add_parBayesC.dat", " ")
# 
# # 
# # make plot
# all_mods <- list.files("output/bglr", pattern = "*24.rda", full.names = TRUE)
# # load(all_mods[[1]])
# # 
# get_marker_effects <- function(chr) {
#         load(paste0("output/bglr/chr_", chr, ".rda"))
#         effs_add <- fm$ETA[[4]]$b
#         effs_roh <- fm$ETA[[3]]$b
#         out <- list(chr, effs_add, effs_roh)
# }
# # 
# all_marker_effs <- map(24, get_marker_effects)
# all_roh_effs <- map_df(all_marker_effs, function(x) {
#         tibble(snp = names(x[[3]]),
#                 eff = x[[2]],
#                 chr = x[[1]])
# })
# # 
# all_roh_effs %>%
#         mutate(num_snp = 1:nrow(.)) %>%
#         ggplot(aes(num_snp, eff, color = as.factor(chr))) + geom_point() + theme_classic()
# 
# 
# load("output/bglr/chr_20.rda")
# out <- read_delim("output/bglr/chr24_ETA_add_parBayesC.dat", " ")
# 
# fm$probs
# 
# #2# Predictions
# # Total prediction
# yHat<-fm$yHat
# tmp<-range(c(y,yHat))
# plot(yHat~y,xlab='Observed',ylab='Predicted',col=2,
#         xlim=tmp,ylim=tmp); abline(a=0,b=1,col=4,lwd=2)
# 
# #3# Godness of fit and related statistics
# fm$fit
# fm$varE # compare to var(y)
# 
# #4# Trace plots
# list.files()
# # Residual variance
# varE<-scan('output/bglr/chr24_varE.dat')
# plot(varE,type='o',col=2,cex=.5,ylab=expression(var[e]));
# abline(h=fm$varE,col=4,lwd=2);
# abline(v=fm$burnIn/fm$thin,col=4)
# # lambda (regularization parameter of the Bayesian Lasso)
# out <- read_delim("output/bglr/chr24_ETA_fixed_b.dat", " ")
# 
# 
# lambda<-scan('output/bglr/chr24_ETA_random_varB.dat')
# plot(lambda,type='o',col=2,cex=.5,ylab=expression(lambda));
# abline(h=fm$ETA[[3]]$lambda,col=4,lwd=2);
# abline(v=fm$burnIn/fm$thin,col=4)
# 
# 
# lambda<-scan('output/bglr/chr24_ETA_fixed_b.dat')
# plot(lambda,type='o',col=2,cex=.5,ylab=expression(lambda));
# abline(h=fm$ETA[[3]]$lambda,col=4,lwd=2);
# abline(v=fm$burnIn/fm$thin,col=4)
# 
# 
# 
# # load("output/bglr/chr_2.rda")
# bHat <- fm$ETA$roh$b
# SD.bHat<- fm$ETA[[3]]$SD.b
# test_stat <- bHat/SD.bHat
# plot(bHat^2, ylab='Estimated Squared-Marker Effect',
#         type='o',cex=.5,col=4,main='Marker Effects')
# # 
# bHat2 <- fm$ETA[[4]]$b
# SD.bHat<- fm$ETA[[4]]$SD.b
# ?p.value
# plot(bHat^2, ylab='Estimated Squared-Marker Effect',
#         type='o',cex=.5,col=4,main='Marker Effects')
# 
# df <- tibble(effs = c(bHat1, bHat2)^2, num_snp = rep(1:length(bHat1), 2),
#         type = rep(c("roh", "add"), each = length(bHat1)))
# ggplot(df, aes(num_snp, effs)) + geom_point() + facet_wrap(~type, nrow = 2) 
# # 
# # 
# # 
# # X_roh <- annual_survival_gwas[grep( "roh", names(annual_survival_gwas))] #%>% 
# #         # as.data.table() %>% 
# #         # impute_mean() %>% 
# #         #as.matrix()
# # 
# # seq(from = 1, to = ncol(x_roh), by = 1000)
# # 
# # library(ggplot2)
# # 
# # sub_X <- split(1:ncol(X_roh), cut_number(1:ncol(X_roh), 100))
# # 
# # impute_mean_across_mat <- function(sub_X_part) {
# #         out <- X_roh[, sub_X_part] %>% as_tibble() %>% impute_mean()
# #         out
# # }
# # 
# # X_roh_imp <- map_df(sub_X, impute_mean_across_mat)
# # 
# # annual_survival_gwas <- annual_survival_gwas[, -grep( "roh", names(annual_survival_gwas))]
# # 
# # fwrite(X_roh, file = "data/annual_survival_gwas_roh.txt")
# # 
# # 
# # ind1 <- which(names(annual_survival_gwas) == "oar3_OAR20_121503")
# # X_add <- annual_survival_gwas[, ind1:(ind1 + ncol(X_roh)-1)] %>% 
# #         impute_mean() %>% 
# #         as.matrix()
# # 
# # #2# Setting the linear predictor
# # ETA<-list( list(~factor(sex)+factor(twin)+factor(age_std)+factor(age2_std),
# #                 data=annual_survival_gwas,model='FIXED'),
# #            list(~factor(id) + factor(sheep_year) + factor(birth_year), data=annual_survival_gwas, model='BRR'),
# #            list(X_roh=X_roh, model='BayesC'),
# #            list(X_add=X_add, model = 'BayesC')
# # )
# # #3# Fitting the model
# # fm<-BGLR(y=y,ETA=ETA, nIter=5000, burnIn=1000, thin = 10, response_type = "ordinal",
# #          saveAt = "output/bglr")
# # save(fm,file='fm.rda')
# # 
# # #1# Estimated Marker Effects & posterior SDs
# # bHat<- fm$ETA[[3]]$b
# # SD.bHat<- fm$ETA[[3]]$SD.b
# # plot(bHat^2, ylab='Estimated Squared-Marker Effect',
# #      type='o',cex=.5,col=4,main='Marker Effects')
# # 
# # bHat<- fm$ETA[[4]]$b
# # SD.bHat<- fm$ETA[[4]]$SD.b
# # plot(bHat^2, ylab='Estimated Squared-Marker Effect',
# #      type='o',cex=.5,col=4,main='Marker Effects')
