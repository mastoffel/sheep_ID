library(Matrix)
library(magrittr)
library(pedigreemm)
library(MasterBayes)
library(lme4qtl)
data(milk)

milk <- within(milk, {
        id <- as.character(id)
        sdMilk <- milk / sd(milk)
})

ids <- with(milk, id[sire %in% 1:3]) # for more ids: c(1:3, 319-321)
milk_subset <- subset(milk, id %in% ids)

milk_subset <- within(milk_subset, {
        herd <- droplevels(herd)
        herd_id <- paste0("herd", seq(1, length(herd)))
})


A_herd <- with(milk_subset, model.matrix(~ herd - 1)) %>% 
        tcrossprod %>% Matrix

rownames(A_herd) <- milk_subset$herd_id
colnames(A_herd) <- milk_subset$herd_id
A_gen <- getA(pedCowsR)

stopifnot(all(ids %in% rownames(A_gen)))
ind <- rownames(A_gen) %in% ids

A_gen <- A_gen[ind, ind]
image(A_herd, main = "A_herd")
image(A_gen, main = "A_gen")

library(nadiv)
load("model_in/sheep_ped.RData")
sheep_ped <- sheep_ped %>% 
        rename(id = ID, dam = MOTHER, sire = FATHER) %>% 
        mutate(dam = str_replace(dam, "F", "999"),
               sire = str_replace(sire, "M", "888")) %>% 
        mutate(id = str_replace(id, "F", "999")) %>% 
        mutate(id = str_replace(id, "M", "888")) %>% 
        mutate_if(is.character, as.integer)

sheep_ped1 <- pedigreemm::pedigree(sire = sheep_ped$sire, dam = sheep_ped$dam, label = as.character(sheep_ped$id))
A_gen <- getA(sheep_ped1)
ind <- rownames(A_gen) %in% early_survival_gwas$id
A_gen <- A_gen[ind, ind]
image(A_gen[1:500, 1:500], main = "A_gen")

dat <- early_survival_gwas %>% mutate(id = as.integer(id))
early_survival_gwas %>% 
        sample_n(500) -> dat
mod1 <- relmatGlmer(survival ~ sex + twin + (1|mum_id) + (1|birth_year) + (1|id), 
                    relmat = list(id = A_gen),
                    data = dat, 
                    family = "binomial")
VarCorr(mod1)
summary(mod1)
