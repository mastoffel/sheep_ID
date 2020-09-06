library(tidyverse)
library(snpStats)

# fitness data 
load("data/survival_mods_data.RData") 
# gwas
gwas_full <- readRDS("output/gwas_res_oar_roh_fac_full.rds")

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map <- full_sample$map 

pos_ests <- gwas_full %>% 
        filter(p.value < 0.05/(39149*2)) %>% 
        filter(estimate >0)
neg_ests <- gwas_full %>% 
        filter(p.value < 0.05/(39149*2)) %>% 
        filter(estimate < 0)
table(as(full_sample$genotypes[, "oar3_OAR19_36946184"], "numeric"))

birth_dates <- fitness_data %>% 
                        group_by(id) %>% 
                        mutate(birth_year = as.numeric(as.character(birth_year))) %>% 
                        summarise(birth_year = mean(birth_year))


get_freqs_pos <- function(year) {
        ids <- birth_dates %>% filter(birth_year == year) %>% .$id
        ids_ind <- rownames(full_sample$genotypes) %in% ids
        snp_stats <- col.summary(full_sample$genotypes[ids_ind, pos_ests$snp.name]) %>% 
                        as_tibble(rownames = "snp") %>% 
                        mutate(birth_year = year)
}

all_freqs_pos <- map(sort(unique(birth_dates$birth_year)), get_freqs_pos) %>% 
                bind_rows()
all_freqs_pos %>% 
        filter(birth_year > 1990) %>% 
ggplot(aes(birth_year,  MAF, color = snp)) +
        geom_line() 



# check wether there is a general trewnd

# check maf of positives
geno_sub <- as_tibble(as(full_sample$genotypes[, pos_ests$snp.name], Class = "numeric"),
                                  rownames = "id")
map(geno_sub[2:ncol(geno_sub)], function(x) as.data.frame(rbind(table(x)))) %>% bind_rows() -> test
snp_stats_sub <- col.summary(full_sample$genotypes[, pos_ests$snp.name])

get_freqs <- function(year, snps_ind) {
        ids <- birth_dates %>% filter(birth_year == year) %>% .$id
        ids_ind <- rownames(full_sample$genotypes) %in% ids
        snp_stats <- col.summary(full_sample$genotypes[ids_ind, snps_ind]) %>% 
                as_tibble(rownames = "snp") %>% 
                mutate(birth_year = year)
}

all_snp_stats <- col.summary(full_sample$genotypes) 
snps_ind <- which(all_snp_stats$MAF < 0.15 & all_snp_stats$MAF > 0.05) %>% 
                sample(10)

#snp_ind <- sample(1:417373, 10)
all_freqs <- map(sort(unique(birth_dates$birth_year)), get_freqs, snps_ind ) %>% 
        bind_rows()

all_freqs %>% 
        bind_rows(all_freqs_pos, .id = "type") %>% 
        filter(birth_year > 1990) %>%
        ggplot(aes(birth_year,  MAF, color = snp, alpha = type)) +
        geom_line(size = 1) +
        theme(legend.position = "none") +
        scale_alpha_manual(values = c(0.2, 1))
               
               
               
           
               
