# # run on server
library(lme4)
library(tidyverse)
library(broom.mixed)
source("theme_clean.R")
library(snpStats)
library(data.table)
library(furrr)

# time saver function
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

# data
load("model_in/fitness_roh_df.RData")
load("model_in/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing_imputation.txt", delim = " ")

# survival data
early_survival <- fitness_data %>% 
        dplyr::rename(birth_year = BIRTHYEAR,
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
        filter(!(id %in% IDs_lots_missing)) %>% 
        filter(age == 0) %>% 
        filter(!is.na(froh_all)) %>% 
        filter(!(is.na(birth_year) | is.na(mum_id))) %>% 
        mutate_at(c("id", "birth_year", "mum_id", "sex"), as.factor) %>% 
        as.data.frame() 

# roh data
file_path <- "output/ROH/roh_nofilt_ram.hom"
roh_lengths <- fread(file_path)

# SNP data
# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/ramb_mapping/sheep_geno_imputed_ram_27092019"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# which chromosome
chr <- 20


# genotypes
snps_sub <- full_sample$map %>% 
                filter(chromosome == chr) %>% 
                .$snp.name
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_sub], Class = "numeric"),
                      rownames = "id")

# roh
roh_sub <- roh_lengths %>% 
        filter(CHR == chr)

# define vectorized seq to work with mutate
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))

# create indices for all rohs
roh_snps <- roh_sub %>% 
        as_tibble() %>% 
        mutate(index1 = as.numeric(match(SNP1, names(geno_sub))),
               index2 = as.numeric(index1 + NSNP - 1)) %>% 
        mutate(all_snps = seq2(from = index1, to = index2)) %>% 
        group_by(IID) %>% 
        summarise(all_snps = list(all_snps)) %>% 
        mutate(all_snps = simplify_all(all_snps)) %>% 
        mutate(IID = as.character(IID)) %>% 
        rename(id = IID)

roh_snps_reord <- geno_sub %>% 
        select(id) %>% 
        left_join(roh_snps) %>% 
        rename(id_roh = id)

roh_mat <- matrix(data = 0, nrow = nrow(geno_sub), ncol = ncol(geno_sub))

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


rm(full_sample)
rm(roh_list)
rm(ro_mat)

# tibble for gwas
early_survival_gwas <- early_survival %>% 
                        dplyr::select(id, survival, sex, twin, birth_year, mum_id) %>% 
                        left_join(geno_sub, by = "id") %>% 
                        left_join(roh_df, by = "id") %>% 
                        as_tibble()


# GWAS
run_gwas <- function(snp, data) {
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + ", snp, "+ ", paste0("roh_", snp), "+ (1|birth_year) + (1|mum_id)"))
        mod <- glmer(formula = formula_snp,
                     data = data, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- tidy(mod)
}

start_time <- Sys.time()
out <- map(snps_sub[1:1000], run_gwas, early_survival_gwas)
end_time <- Sys.time()

end_time-start_time

microbenchmark::microbenchmark(map(snps_sub1, run_gwas, early_survival_gwas))

gwas_p <- map_df(out, function(x) tibble(snp = tidy(x)$term[4], pval = tidy(x)$p.value[5]))

gwas_p_full <- gwas_p %>% 
        dplyr::left_join(snp_maps, by = "snp") %>% 
        dplyr::mutate(p_adj = p.adjust(pval, method = "fdr")) 





# overall model
mod1 <- glmer(survival ~ 1 + froh_all + sex + twin + (1|birth_year) + (1|mum_id),
              data = early_survival, family = "binomial",
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))

confint(mod1)
summary(mod1)
VarCorr(mod1)
tidy(test)
dotplot(ranef(mod1,condVar=TRUE))
methods(class="merMod")
plot(mod1)

# make plink file
system(paste0("~/programs/plink --bfile ../sheep/data/SNP_chip/ramb_mapping/sheep_geno_imputed_ram_27092019 --sheep --out data/gwas/sheep_chr20 ",
              "--mind 0.1 --chr 20 --recode --make-bed --indep-pairwise 55000 55000 0.95 ")) 

# keep snps
snps_keep <- readLines("data/gwas/sheep_chr20.prune.in")

# plink name
sheep_plink_name <- "data/gwas/sheep_chr20"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

sheep_geno <- as(full_sample$genotypes, Class = "numeric") %>% 
                as_tibble(rownames = "id") %>% 
                dplyr::select(id, snps_keep)

sheep_geno[sheep_geno == 2] <- 0 # make homozygous 0 and het 1

gwas_dat <- early_survival %>% 
        left_join(sheep_geno, by = "id")

names(gwas_dat)[1:100]
snps_gwas <- names(gwas_dat)[59:length(names(gwas_dat))]

# get mapping positions
# snp_maps <- full_sample$map %>% 
#         dplyr::select(snp.name, position) %>% 
#         dplyr::rename(snp = snp.name) %>% 
#         filter(snp %in% snps_gwas) %>% 
#         as_tibble() %>% 
#         mutate(snp = str_replace(snp, "\\.", ""))

run_gwas <- function(snp, data) {
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + ",snp, "+ (1|birth_year) + (1|mum_id)"))
        out <- glmer(formula = formula_snp,
              data = data, family = "binomial",
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out
}
#snps_gwas_sub <- sample(snps_gwas, 10)
gwas_out <- map(snps_gwas_sub, run_gwas, gwas_dat)

tidy(gwas_out[[1]])

gwas_p <- map_df(gwas_out, function(x) tibble(snp = tidy(x)$term[4], pval = tidy(x)$p.value[4]))

gwas_p_full <- gwas_p %>% 
        dplyr::left_join(snp_maps, by = "snp") %>% 
        dplyr::mutate(p_adj = p.adjust(pval, method = "fdr")) 

        

hist(gwas_p$p_adj)
sum(gwas_p$p_adj < 0.05, na.rm = TRUE)
