# evaluate whether additive genetic variation changes things
library(tidyverse)
library(data.table)
library(snpStats)
top_snps <- read_delim("output/top_roh_snps_gwas_testset.txt", delim = " ") %>% 
                distinct(snp.name, .keep_all = TRUE)
chrs <- unique(top_snps$chromosome)

# data
load("data/fitness_roh_df.RData")
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

snps_map_sub <- full_sample$map #%>% filter(chromosome %in% chrs)

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
        mutate(gwas_snp = map(snps_roh, function(x) {
                snp_in_roh <- top_snps$snp.name %in% x
                out <- if (any(snp_in_roh)) {
                    out <- top_snps$snp.name[snp_in_roh]
                } else {
                    out <- NA
                }
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
colnames(roh_mat) <- c("id", top_snps$snp.name)

# set 1 where SNP is in an roh
# make a tibble like for genotypes but with 0/1 for whether a SNP is in ROH or not
roh_df <- pmap_df(roh_snps_reord[c("id_roh", "all_snps")], function(id_roh, all_snps) {
        
        df <- as_tibble(t(as.numeric(c(id_roh, rep(0, ncol(geno_sub) - 1)))))
        names(df) <- c("id_roh", unique(top_snps$snp.name))
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
ann_survival_top_snps <- annual_survival %>% 
        dplyr::select(id, survival, sex, twin, birth_year, sheep_year, mum_id, age, age2, froh_all) %>% 
        left_join(geno_sub, by = "id") %>% 
        left_join(roh_df, by = "id") %>% 
        as_tibble()

# make some space
rm(full_sample)
rm(roh_list)
rm(roh_mat)

# pcas 
pcs <- read_delim("data/ann_surv_pca_testset80.txt", " ", col_names = TRUE) %>% 
        mutate(id = as.character(id))

out <- ann_survival_top_snps %>% 
        left_join(pcs)

write_delim(out, "output/annual_survival_top_snps_pca_testset.txt", delim = " ")








