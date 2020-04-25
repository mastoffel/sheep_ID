library(tidyverse)

hd <- read_delim("data/wrongly_mapped_hd.txt", " ") %>% 
        select(chr_new, snp.name, pos_new, allele.1.x, allele.2.x, chr_old, pos_old)
ld <- read_delim("data/wrongly_mapped_ld.txt", " ") %>% 
        select(chr_new, snp.name, pos_new, allele.1.x, allele.2.x, chr_old, pos_old)

all <- hd %>% 
        bind_rows(ld, .id = "chip") %>% 
        mutate(chip = ifelse(chip == 1, "HD", "50K")) %>% 
        select(chip, snp.name, chr_new, chr_old, pos_new, pos_old, everything()) %>% 
        setNames(c("snp_chip", "snp_name", "chr_rambV1", "chr_oarV3_1", "pos_rambV1", "pos_oarV3_1", "allele1", "allele2")) %>% 
        write_delim("data/RambV1_problematic_SNPs.txt")
