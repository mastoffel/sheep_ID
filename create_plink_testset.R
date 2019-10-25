# creates a testset containing 80% of individuals for gwas
library(tidyverse)
library(snpStats)
# data
load("data/fitness_roh_df.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")
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
        mutate_at(c("id", "birth_year", "sex", "sheep_year", "survival"), as.factor) %>% 
        mutate(age2 = age^2) %>% 
        mutate(age_std = as.numeric(scale(age)),
               age2_std = as.numeric(scale(age2))) %>% 
        as.data.frame() 


# use only training set (80% of individuals)
sample_frac_groups = function(tbl, size, replace = FALSE, weight = NULL) {
        # regroup when done
        grps = tbl %>% groups %>% lapply(as.character) %>% unlist
        # check length of groups non-zero
        keep = tbl %>% summarise() %>% ungroup() %>% sample_frac(size, replace, weight)
        # keep only selected groups, regroup because joins change count.
        # regrouping may be unnecessary but joins do something funky to grouping variable
        tbl %>% right_join(keep, by=grps) %>% group_by(.dots = grps)
}

set.seed(7336)
all_ids <- unique(as.character(early_survival$id)) %>% 
                sample(length(.) * 0.8, replace = FALSE)

# write_lines(early_survival$index, "output/ind_index_test.txt")
#write_lines(unique(early_survival$id), "data/ind_testset_80.txt")

write_delim(tibble(1, all_ids), path = "data/ind_testset_80.txt", 
            delim = " ", col_names = FALSE)

system(paste0("~/programs/plink --bfile data/sheep_geno_imputed_ram_27092019_pruned ",
              "--keep data/ind_testset_80.txt --make-bed --out data/sheep_geno_imputed_ram_27092019_pruned_cv80")) # --keep data/ind_testset_80.txt

# repeat pca
# PCA for GWAS
system(paste0("~/programs/plink --bfile data/sheep_geno_imputed_ram_27092019_pruned --keep data/ind_testset_80 ",
              "--sheep --make-grm-bin --pca --out output/sheep_pca_cv80"))

# gcta64  --grm test --keep test.indi.list  --pca 20  --out test
# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019_pruned"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
chr <- 20
inds <- full_sample$fam

snps_sub <- full_sample$map %>% 
        filter(chromosome == chr) %>% 
        .$snp.name

geno_sub <- as_tibble(as(full_sample$genotypes[, snps_sub], Class = "numeric"),
                      rownames = "id")

nrow(full_sample$fam)
str(full_sample$genotypes)
nrow(full_sample$map)
