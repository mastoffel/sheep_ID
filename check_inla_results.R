# evaluate whether additive genetic variation changes things

library(tidyverse)
library(data.table)
library(snpStats)

top_snps <- read_delim("output/top_snps_gwas.txt", delim = " ")
snp <- "oar3_OAR11_8762236"

chr <- 15

# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

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

snps_map_sub <- full_sample$map %>% 
        filter(chromosome == chr) 

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
        filter(chromosome == chr) %>% 
        .$snp.name
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_sub], Class = "numeric"),
                      rownames = "id")

# subset roh 
roh_sub <- roh_lengths %>% filter(CHR == chr) %>% filter(KB > 1000)

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
        dplyr::select(id, survival, sex, twin, birth_year, sheep_year, mum_id, age_std, age2_std) %>% 
        left_join(geno_sub, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        as_tibble()



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
