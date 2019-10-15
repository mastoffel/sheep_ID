library(tidyverse)
library(snpStats)
source("theme_clean.R")
gwas_files <- list.files("output/gwas_full_roh", pattern = "*.rds", full.names = TRUE)

# extract results
all_gwas <- purrr::map(gwas_files, readRDS) %>% 
            purrr::flatten() %>% 
            # only extract results
            .[seq(1,length(.),by=2)] %>% 
           # remove snps that didnt work
            purrr::compact()

# check errors
# all_gwas[seq(2,length(all_gwas),by=2)] %>% 
#         compact()

# get roh pval
gwas_res <- map_df(all_gwas, function(x) x %>% .[c(6,7), ] %>% 
                           dplyr::select(term, estimate, p.value))

# check whether roh didnt work anywhere
gwas_res[which(str_detect(gwas_res$term, "sd")), ]

# remove those 
gwas_res <- gwas_res[-which(str_detect(gwas_res$term, "sd")), ]

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019_pruned"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map <- full_sample$map 
table(full_sample$map$chromosome, useNA = "always")


# gwas_res_roh <- map_df(all_gwas, function(x) x %>% .[c(5), ] %>% dplyr::select(term,estimate, p.value))
# gwas_res_snp <- map_df(all_gwas, function(x) x %>% .[c(4), ] %>% dplyr::select(term,estimate, p.value))

# check how often roh var couldnt be estimated
# gwas_roh <- gwas_res %>% filter(str_detect(term, "roh"))
# gwas_roh[which(!(gwas_snp$term %in% full_sample$map$snp.name)), ]

all_gwas[[11842]]
# put into df
gwas_full <- gwas_res %>%
        rename(snp.name = term) %>%
        # roh variable couldnt be estimated in some models, so the next
        # row was extracted from the tidy output
        filter(!str_detect(snp.name, "sd")) %>% 
        mutate(groups = ifelse(str_detect(snp.name, "roh"), "roh", "add")) %>% 
        mutate(snp.name = str_replace(snp.name, "roh_", "")) %>%
        left_join(snps_map) 
# 
qqman::qq(gwas_full[gwas_full$groups == "add", ]$p.value)
qqman::qq(gwas_full[gwas_full$groups == "roh", ]$p.value)

# manhattan
## computing new x axis
gwas_roh <- gwas_full %>% 
                        #filter(chromosome %in% c(4:26)) %>% 
                        group_by(groups) %>% 
                        arrange(chromosome, position) %>% 
                        dplyr::mutate(tmp = 1, cumsum.tmp = cumsum(tmp))
## calculating x axis location for chromosome label
med.dat <- gwas_roh %>% dplyr::group_by(groups, chromosome) %>% 
                dplyr::summarise(median.x = median(cumsum.tmp))


ggplot(data = gwas_roh) + 
        geom_point(aes(x = cumsum.tmp, y = -log10(p.value), color = chromosome %%2 == 0),
                   size = 1) + ## create alternate coloring
        geom_hline(yintercept = -log10(0.05/28946)) + ## add horizontal line
        scale_x_continuous(breaks = med.dat$median.x, labels = med.dat$chromosome) + ## add new x labels 
        guides(colour=FALSE) +  ## remove legend
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
        scale_color_manual(values = c(gray(0.5), gray(0))) +## instead of colors, go for gray
        theme_clean() +
        facet_wrap(groups~., nrow = 2)


gwas_roh %>% arrange(p.value) %>% filter(groups == "roh")



# prepare genotypes for simpleM
snps_geno <- full_sample$map 
sheep_geno <- as(full_sample$genotypes[, snps_geno$snp.name], Class = "numeric")
missings <- rowSums(is.na(sheep_geno))
sheep_geno <- sheep_geno[missings < (0.01*ncol(sheep_geno)), ] # remove inds with more than 1% missing
sheep_geno[is.na(sheep_geno)] <- sample(c(0,1,2), 1)
sheep_geno_t <- t(sheep_geno)

library(data.table)
fwrite(sheep_geno_t, file = "data/geno_mat_simpleM_allchr.txt", col.names = FALSE, row.names = FALSE)

# library(synbreed)
# ?codeGeno
# impute_matrix <- function(sheep_geno) {
#         col_means <- colMeans(sheep_geno, na.rm = TRUE)
# }

gwas_roh %>% arrange(p.value) %>% 
        filter(groups == "roh") %>% 
        .[1:50, ] %>% write_delim("output/top_snps_gwas.txt")







