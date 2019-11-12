library(tidyverse)
library(lme4)
library(performance)

# get df
surv_df <- read_delim("output/annual_survival_top_snps_pca_testset.txt", delim = " ") %>% 
                                mutate(age_std = as.numeric(scale(age)),
                                       age2_std = as.numeric(scale(age2))) %>% 
                                mutate(sex = as.factor(sex), twin = as.factor(twin),
                                       birth_year = as.factor(birth_year),
                                       sheep_year = as.factor(sheep_year),
                                       id = as.factor(id)) %>% 
                                filter(!is.na(sex) & !is.na(twin))

ids_test <- read_delim("data/ids_annual_survival_testset80.txt", " ", col_names = FALSE) %>% 
                rename(id = X2)

surv_df <- surv_df %>% mutate(trainset = ifelse(id %in% ids_test$id, "yes", "no"))



# split 
library(caret)
library(lme4)
library(BGLR)
library(tidyimpute)
survival_train <- surv_df %>% filter(trainset == "yes")
survival_test <- surv_df %>% filter(trainset == "no")

roh_snps <- names(surv_df)[grep("roh", names(surv_df))[-1]]
add_snps <- str_replace(roh_snps, "roh_", "")

# BGLR prediction
X_roh <- surv_df[roh_snps] %>% 
        impute_mean() %>% 
        as.matrix()
X_add <- surv_df[add_snps] %>% 
        impute_mean() %>% 
        as.matrix()

# set y to NA for testset predictions
surv_df$survival[surv_df$trainset == "no"] <- NA

#2# Setting the linear predictor
ETA<-list(fixed = list(~factor(sex)+factor(twin)+age_std+age2_std,
                       data=surv_df,model='FIXED'),
          id = list(~factor(id), data=surv_df, model='BRR'),
          sheep_year = list(~factor(sheep_year), data=surv_df, model='BRR'),
          birth_year = list(~factor(birth_year), data=surv_df, model='BRR'),
          roh = list(X_roh=X_roh, model='BayesC'),
          add = list(X_add=X_add, model = 'BayesC') # 
)
y <- surv_df$survival
#3# Fitting the model
fm <- BGLR(y=y,ETA=ETA, nIter=5000, burnIn=1000, thin = 50, 
            response_type = "ordinal",
            # additional iterations with the following two lines
            #BGLR_ENV = paste0(output_folder, run_name, "BGLR_ENV.RData"), # default NULL
            #newChain = FALSE, # default TRUE
            # where to save
            saveAt = "output/bglr_pred/try") 

plot(fm$ETA$roh$b^2)
hist(fm$yHat[surv_df$trainset == "no"])

roh_snps <- names((fm$ETA$roh$b^2)[fm$ETA$roh$b^2 > 0.001])
add_snps <- str_replace(roh_snps, "roh_", "")
        
lme_form <- reformulate(c(roh_snps, add_snps, paste0("pc", 1:7), "sex", "twin", "age_std", "age2_std",  "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
mod1 <- glmer(lme_form, data = survival_train, family = "binomial")
summary(mod1)
?predict


pred_mod <- predict(mod1, newdata = survival_test, allow.new.levels = TRUE, type = "response")

#%>% 
        #         mutate(surv = ifelse(survival > 0.5, 1, 0)) %>% 
        # mutate(surv = as.factor(surv),
        #        survival_org = as.factor(survival_org))

confusionMatrix(data = newdat$surv, reference = newdat$survival_org)

newdat <- survival_test %>% rename(survival_org = survival) %>% 
        mutate(survival = predict(mod1, newdat, 
                                  allow.new.levels = TRUE, type = "response")) %>% 
        mutate(surv = ifelse(survival > 0.5, 1, 0)) %>% 
        mutate(surv = as.factor(surv),
               survival_org = as.factor(survival_org))

confusionMatrix(data = newdat$surv, reference = newdat$survival_org)

library(merTools)
predictInterval(mod1, newdata = survival_test,
                level = 0.95, n.sims = 5)

mm <- model.matrix(terms(mod1),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod1),mm))









library(BGLR)
y <- surv_df$survival
X <- surv_df[, grep( "roh", names(surv_df))[-1]] %>% 
        impute_mean()
X_add <- surv_df %>% dplyr::select(oar3_OAR15_63743364:oar3_OAR15_64194091) %>% 
        impute_mean()
#2# Setting the linear predictor
ETA<-list( list(~factor(sex)+factor(twin)+factor(age_std)+factor(age2_std),
                data=surv_df,model='FIXED'),
           list(~factor(id) + factor(sheep_year) + factor(birth_year), data=surv_df, model='BRR'),
           list(X=X, model='BL'),
           list(X_add=X_add, model = 'BL')
)
#3# Fitting the model
fm<-BGLR(y=y,ETA=ETA, nIter=1000, burnIn=200)
save(fm,file='fm.rda')

#1# Estimated Marker Effects & posterior SDs
bHat<- fm$ETA[[3]]$b
SD.bHat<- fm$ETA[[3]]$SD.b
plot(bHat^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')

bHat<- fm$ETA[[4]]$b
SD.bHat<- fm$ETA[[4]]$SD.b
plot(bHat^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')

#3# Godness of fit and related statistics
fm$fit
fm$varE # compare to var(y)
#4# Trace plots
list.files()
# Residual variance
varE<-scan('varE.dat')
plot(varE,type='o',col=2,cex=.5,ylab=expression(var[e]));
abline(h=fm$varE,col=4,lwd=2);
abline(v=fm$burnIn/fm$thin,col=4)
# lambda (regularization parameter of the Bayesian Lasso)
lambda<-scan('ETA_3_lambda.dat')
plot(lambda,type='o',col=2,cex=.5,ylab=expression(lambda));
abline(h=fm$ETA[[3]]$lambda,col=4,lwd=2);
abline(v=fm$burnIn/fm$thin,col=4)




# how much variation do all ROH snps explain?
# time saver function for modeling
nlopt <- function(par, fn, lower, upper, control) {
        .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
                                  opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                              maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
        list(par = res$solution,
             fval = res$objective,
             conv = if (res$status > 0) 0 else res$status,
             message = res$message
        )
}

# get roh snps
roh_snps <- names(early_survival_top_snps)[str_detect(names(early_survival_top_snps), "roh")][-1]

# check
test <- survival_top_snps %>% 
        filter(!is.na(survival_top_snps$survival)) %>% 
        summarise_at(vars(contains("roh")), function(x) length(table(x)))

# with roh_snps
glme_form <- reformulate(c("sex", "twin", "age_std", "age2_std", roh_snps, "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
mod1 <- glmer(glme_form, data = survival_top_snps, family = "binomial",
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod1)
r2_nakagawa(mod1)

glme_form2 <- reformulate(c("sex", "twin", "age_std", "age2_std",  "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
mod2 <- glmer(glme_form2, data = survival_top_snps, family = "binomial",
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
summary(mod2)
library(performance)
r2_nakagawa(mod2)


# get all inla mods
top_snps_gwas <- read_delim("output/top_snps_gwas_pca.txt", delim = " ")
mods <- list.files("output/inla_mods_gwas_top/", full.names = TRUE)
all_inla <- map(mods, read_rds)
all_eff <- map_df(all_inla, function(mod) {
        mod <- as_tibble(mod$summary.fixed, rownames = "pred")
        mod[c(7), ]
        }) %>% 
        mutate(pred = str_replace(pred, "roh_", ""))

top_snps_gwas %>% 
        select(snp.name, estimate) %>% 
        rename(pred=snp.name) %>% 
        inner_join(all_eff, by = "pred") %>% 
        ggplot(aes(estimate, mean)) + geom_point() + geom_smooth(method = "lm")




















library(INLA)
# inla
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

add_index_inla <- function(traits) {
        Ndata <- dim(traits)[1]
        traits$IndexA <- rep(0, times = Ndata)
        for(i in 1:Ndata) {
                traits$IndexA[i] <- which(ainv_map[,1]==traits$id[i])
        }
        traits
}

early_survival_gwas <- add_index_inla(early_survival_gwas)
# for all other random effects too
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))

formula_surv <- as.formula(paste('survival ~ oar3_OAR15_63743364 + roh_oar3_OAR15_63743364 + sex + age_std + age2_std + twin + 1', 
                                 'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                 'f(id, model = "iid", hyper = list(prec = prec_prior))',
                                 # 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                                 'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))

mod_inla <- inla(formula=formula_surv, family="binomial",
                 data=early_survival_gwas, 
                 control.compute = list(dic = TRUE))
summary(mod_inla)

formula_surv_norel <- as.formula(paste('survival ~ oar3_OAR15_63743364 + roh_oar3_OAR15_63743364 + sex + age_std + age2_std + twin + 1', 
                                       'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                                       'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                                       'f(id, model = "iid", hyper = list(prec = prec_prior))', sep = " + "))
# 'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
#'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + "))

mod_inla2 <- inla(formula=formula_surv_norel, family="binomial",
                  data=early_survival_gwas, 
                  control.compute = list(dic = TRUE))
summary(mod_inla2)
