# evaluate whether additive genetic variation changes things
library(tidyverse)
library(data.table)
library(snpStats)
top_snps <- read_delim("output/top_snps_gwas_train08.txt", delim = " ")
chrs <- unique(top_snps$chromosome)

# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing_imputation.txt", delim = " ")

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

snps_map_sub <- full_sample$map #%>% filter(chromosome %in% chrs)

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
        filter(chromosome %in% chrs) %>% 
        .$snp.name
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_sub], Class = "numeric"),
                      rownames = "id")

# subset roh 
roh_sub <- roh_lengths %>% filter(KB > 1000) %>% filter(CHR %in% chrs)

# define vectorized seq to work with mutate
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

# create indices for all rohs
roh_snps <- roh_sub %>% 
        as_tibble() %>% 
        # sample_frac(0.01) %>% 
        mutate(index1 = as.numeric(match(SNP1, names(geno_sub))),
               index2 = as.numeric(index1 + NSNP - 1)) %>% 
        mutate(all_snps = seq2(from = index1, to = index2)) %>% 
        mutate(snps_roh = map(all_snps, function(x) names(geno_sub)[x])) %>% 
        # it might be that this results in two snps for a roh
        mutate(gwas_snp = map_chr(snps_roh, function(x) {
                snp_in_roh <- top_snps$snp.name %in% x
                out <- ifelse(any(snp_in_roh), top_snps$snp.name[snp_in_roh], NA)
                out
        })) %>% 
        filter(!is.na(gwas_snp)) %>% 
        group_by(IID) %>% 
        summarise(all_snps = list(gwas_snp)) %>% 
        mutate(all_snps = simplify_all(all_snps)) %>% 
        mutate(IID = as.character(IID)) %>% 
        rename(id = IID)

# filter geno_sub
geno_sub <- geno_sub %>% 
        dplyr::select(id, top_snps$snp.name) 

# join roh_snps
roh_snps_reord <- geno_sub %>% 
       # dplyr::select(id, top_snps$snp.name) %>% 
        left_join(roh_snps, by = "id") %>% 
        dplyr::rename(id_roh = id)

# prepare roh yes/no matrix
roh_mat <- matrix(data = 0, nrow = nrow(geno_sub), ncol = ncol(geno_sub))
colnames(roh_mat) <- c(id, top_snps$snp.name)

# set 1 where SNP is in an roh
# make a tibble like for genotypes but with 0/1 for whether a SNP is in ROH or not
roh_df <- pmap_df(roh_snps_reord[c("id_roh", "all_snps")], function(id_roh, all_snps) {
        df <- as_tibble(t(as.numeric(c(id_roh, rep(0, ncol(geno_sub) - 1)))))
        names(df) <- c("id_roh", top_snps$snp.name)
        if (!is.null(all_snps)){
                df[, all_snps] <- 1
        } else {
                return()
        }
        df
        #as_tibble(df)
}) 

names(roh_df) <- c("id", paste0("roh_", names(geno_sub)[-1]))
roh_df <- mutate(roh_df, id = as.character(id))

# join additive and roh data to survival for gwas
early_survival_top_snps <- early_survival %>% 
        dplyr::select(id, survival, sex, twin, birth_year, sheep_year, mum_id, age_std, age2_std) %>% 
        left_join(geno_sub, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        as_tibble()

# make some space
rm(full_sample)
rm(roh_list)
rm(roh_mat)

write_delim(early_survival_top_snps, "output/early_survival_top_snps_train08.txt", delim = " ")






# get df
early_survival_top_snps <- read_delim("output/early_survival_top_snps.txt", delim = " ") %>% 
                                        mutate(roh_all = rowSums(select_at(., vars(contains("roh"))), na.rm = TRUE))
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

library(lme4)
library(performance)
roh_snps <- names(early_survival_top_snps)[str_detect(names(early_survival_top_snps), "roh")]
glme_form <- reformulate(c("sex", "twin", "age_std", "age2_std", roh_snps, "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
mod1 <- glmer(glme_form, data = early_survival_top_snps, family = "binomial",
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))




# split 
library(caret)
library(lme4)
roh_snps <- paste0("roh_", top_snps$snp.name)
trainIndex <- createDataPartition(early_survival_top_snps$survival, p = .7, 
                                  list = FALSE, 
                                  times = 1)

survival_train <- early_survival_top_snps[trainIndex, ]
survival_test <- early_survival_top_snps[-trainIndex, ]
survival_test %>% filter((birth_year%in%survival_train$birth_year)) %>% 
                  filter((sheep_year %in% survival_train$sheep_year))


lme_form <- reformulate(c(roh_snps, "sex", "twin", "age_std", "age2_std",  "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
mod1 <- glmer(lme_form, data = survival_train, family = "binomial")
summary(mod1)
?predict

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
