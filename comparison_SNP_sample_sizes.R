library(snpStats)

# Oar3.1 mapping
# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/oar31_mapping/Plates_1to87_QC3"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
ncol(full_sample$genotypes) # 39368
nrow(full_sample$genotypes) # 7700

# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/oar31_mapping/Plates_1-2_HD_QC2"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
ncol(full_sample$genotypes) # 430702
nrow(full_sample$genotypes) # 7700


# Ramb mapping
# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/ramb_mapping/Plates_1to87_QC4_ram"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample_ld <- read.plink(sheep_bed, sheep_bim, sheep_fam)
ncol(full_sample$genotypes) # 38300
nrow(full_sample$genotypes) # 7700

# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
#full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
full_sample_hd <- read.plink(sheep_bed, sheep_bim, sheep_fam)
ncol(full_sample$genotypes) # 411800
nrow(full_sample$genotypes) # 7700


# merged RAM
# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/ramb_mapping/merged_sheep_geno_ram_snpfilt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
ncol(full_sample$genotypes) # 411800
nrow(full_sample$genotypes) # 7700

ld <- colnames(full_sample_ld$genotypes)
hd <- colnames(full_sample_hd$genotypes)
merged <- colnames(full_sample$genotypes)

# genotyped on both chips and in merged
sum(merged %in% ld[(ld %in% hd)])

# imputed
sheep_plink_name <- "data/sheep_geno_imputed_ram_400k_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
ncol(full_sample$genotypes) # 411800
nrow(full_sample$genotypes) # 7700

summary(full_sample$genotypes)
