library(data.table)
library(tidyverse)
mismatches <- fread("../sheep_ID_paper/submission_NatComm/final_submission/SNP_data_submission_old/mismatches.txt",
                    header = FALSE, sep = " ") %>% 
              as_tibble() %>% 
              select(V4, V6, V9, V17) %>% 
              setNames(c("chr", "pos", "ref", "exp")) %>% 
              mutate(chr = str_replace(chr, ",", ""),
                     pos = as.numeric(str_replace(pos, ",", "")),
                     ref = str_replace_all(ref, "'", ""),
                     exp = str_replace_all(exp, "'", ""),
                     chr = as.numeric(str_replace(chr, "OAR", "")))

# only ld chip
snps_ld <- data.table::fread("../sheep/data/SNP_chip/oar31_mapping/20200323_Plates_1to91.bim") %>% 
                setNames(c("chr", "snp", "cM", "pos", "alt_ld", "ref_ld")) %>% 
                select(chr, pos, ref_ld, alt_ld) %>% 
                filter(chr > 0)

snps_hd <- data.table::fread("../sheep/data/SNP_chip/oar31_mapping/Plates_1-2_HD_QC2_snpfile_for_merge.bim") %>% 
        setNames(c("chr", "snp", "cM", "pos", "alt_hd", "ref_hd")) %>% 
        select(chr, pos, ref_hd, alt_hd) %>% 
        filter(chr > 0)

# read genotypes, data.table does this correctly
# i.e. jumping over the first 5 lines
gt <- data.table::fread("../sheep_ID_paper/submission_NatComm/final_submission/SNP_data_submission_old/sheep_geno_oar31.vcf", 
                        skip = 32, select = 1:5)
all <- gt %>% 
        as_tibble() %>% 
        setNames(c("chr", "pos", "id", "ref_vcf", "alt_vcf")) %>% 
        left_join(mismatches, by = c("chr", "pos")) %>% 
        left_join(snps_ld) %>% 
        left_join(snps_hd) %>% 
        filter(!is.na(ref))



all_sub <- all %>% 
   drop_na() 

sum(all_sub$alt_vcf == all_sub$exp)/nrow(all_sub)

sum(all$alt_ld == all$exp, na.rm = TRUE)/sum(!is.na(all$alt_ld))

all %>% 
        filter(!is.na(ref_ld),
               ref_ld == exp)
