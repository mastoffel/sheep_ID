# create data for EVA submission

# to make vcfs
# plink --bfile merged_sheep_geno_oar31 --recode vcf --out sheep_geno_oar31 --sheep --snps-only just-acgt
# ~/Downloads/tabix-0.2.6/bgzip sheep_geno_oar31.vcf
# ~/Downloads/tabix-0.2.6/tabix -p vcf sheep_geno_oar31.vcf.gz
# check
# ~/Downloads/vcf_validator_macos -i sheep_geno_oar31.vcf


library(snpStats)
library(tibble)
library(writexl)

### merged genos
# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/oar31_mapping/merged_sheep_geno_oar31"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)


df <- tibble(id = full_sample$fam$member, title = paste0("Soah_sheep_sample", 1:nrow(full_sample$fam)))
df %>% 
  write_xlsx(path = "output/ids_for_eva.xlsx")



### imputed genos
# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/oar31_mapping/sheep_geno_imputed_oar31_17052020"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)


df <- tibble(id = full_sample$fam$member, title = paste0("Soah_sheep_sample", 1:nrow(full_sample$fam)))
df %>% 
        write_xlsx(path = "output/ids_for_eva_imputed.xlsx")
