library(windowscanr)
library(data.table)
library(tidyverse)

#~~~ ROH density
hom_sum <- fread("output/ROH/roh_ram.hom.summary") # ROH_surv_subset/
#hom_sum_sub <- hom_sum %>% filter(UNAFF == 0)

hom_sum <- hom_sum %>%
        mutate(MB = BP / 1000000,
               KB = BP / 1000,
               index = 1:nrow(.))

#check at every SNP that calling an ROH would be possible
roh_possible <- function(snp_ind, hom_sum_df, roh_snps = 30, roh_kb = 600) {

 # hom_sum_df <- hom_sum_df[hom_sum_df$CHR == chr, ]
  snp_min <- snp_ind-roh_snps+1
  snp_max <- snp_ind+roh_snps-1
  if (snp_min < 1) snp_min <- 1
  if (snp_max > nrow(hom_sum_df)) snp_max <- nrow(hom_sum_df)
  out <- winScan(x = hom_sum_df[snp_min:snp_max, ],
                 values = "KB",
                 win_size = roh_snps,
                 win_step = 1,
                 funs = c("max", "min"))
  out$KB_diff <- out$KB_max - out$KB_min
  res <- ifelse(any(out[1:roh_snps, ]$KB_diff <= roh_kb), "yes", "no")

}

# check at which SNPs there is 
pot_snps <- hom_sum %>% filter(UNAFF == 0) %>% .$SNP
pot_snps_ind <- which(hom_sum$SNP %in% pot_snps)
snps_ok <- data.frame(pot_snps, roh_possible = map_chr(pot_snps_ind, roh_possible, hom_sum))

write_delim(snps_ok_df, "output/snps_that_can_have_roh_400k")

