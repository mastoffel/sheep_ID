# make SVDs for annual survival
# # run on server
library(tidyverse)
library(broom.mixed)
library(tidyimpute)
library(snpStats)
library(data.table)
library(furrr)

output_folder <- "/exports/eddie/scratch/mstoffel/svd/"

# data
load("data/fitness_roh_df.RData")
load("data/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# pcas 
pcs <- read_delim("data/ann_surv_pca.txt", " ", col_names = TRUE) %>% 
        mutate(id = as.character(id))

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

snps_map_sub <- full_sample$map 

# survival data
annual_survival <- fitness_data %>% 
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
        mutate_at(c("id", "birth_year", "sex", "sheep_year", "survival"), as.factor) %>% 
        mutate(age2 = age^2) %>% 
        mutate(age_std = as.numeric(scale(age)),
               age2_std = as.numeric(scale(age2))) %>% 
        as.data.frame() 



# prepare additive genotypes subset
snps_sub <- full_sample$map %>% 
        #filter(chromosome == chr) %>% 
        .$snp.name
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_sub], Class = "numeric"),
                      rownames = "id")

# subset roh 
roh_sub <- roh_lengths %>% 
        #filter(CHR == chr)  %>% 
        filter(KB > 1000)

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
annual_survival_gwas <- annual_survival %>% 
        dplyr::select(id) %>% 
        left_join(geno_sub, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        # remove all duplicates
        dplyr::distinct(id, .keep_all = TRUE) %>% 
        as_tibble()

# vector of ids
ids <- annual_survival_gwas$id
# now make svds
roh_snps <- grep("roh", names(annual_survival_gwas), value = TRUE)
add_snps <- stringr::str_replace(snp_names_roh, "roh_", "")

run_svd <- function(snps) {
        
        gen_mat <- annual_survival_gwas[snps] %>% 
                impute_mean() %>% 
                as.matrix() %>% 
                scale()
        # svd
        svd_gen_mat <- svd(gen_mat)
        S_mat <- diag(svd_gen_mat$d)
        
        # name
        out_name <- ifelse(any(str_detect(roh_snps, "roh")), "roh", "add")
        
        # save
        fwrite(svd_gen_mat$u, paste0(output_folder, out_name, "_U.txt"), nThread = 1, col.names = FALSE, sep = " ")
        fwrite(svd_gen_mat$v, paste0(output_folder, out_name, "_V.txt"), nThread = 1, col.names = FALSE, sep = " ")
        fwrite(S_mat, paste0(output_folder, out_name, "_S.txt"), nThread = 1, col.names = FALSE, sep = " ")
        
        # US mat
        gen_US <- svd_gen_mat$u %*% S_mat
        fwrite(gen_US, paste0(output_folder, out_name, "_US.txt"), nThread = 1, col.names = FALSE, sep = " ")
        
}

# set up plan
plan(multiprocess, workers = 2)
future_map(list(roh_snps, add_snps), run_svd)





