# make SVDs for annual survival
# # run on server
library(tidyverse)
library(data.table)

output_folder <- "/exports/eddie/scratch/mstoffel/svd/"
#output_folder <- "output/svd/"
# name
out_name <- "snps" # snps roh
out_name <- "roh"
# non genetic variables
annual_survival_gwas <- fread("data/annual_survival_gwas_vars.txt")
gen_dat <- fread(paste0("data/annual_survival_gwas_", out_name, ".txt"))
#add <- fread("data/annual_survival_gwas_snps.txt", select = 1:10)

# get unique ids
ids <- annual_survival_gwas$id[!duplicated(annual_survival_gwas$id)]

# extract unique ids from geno mat
gen_dat <- gen_dat[!duplicated(annual_survival_gwas$id), ]
gen_dat_scaled <- scale(as.matrix(gen_dat))

# get snps without any variation and remove from data
snps_without_var <- which(attr(gen_dat_scaled, "scaled:scale") == 0)

# remove non-variant snps
if (length(snps_without_var) > 0) {
        gen_dat_scaled <- gen_dat_scaled[, -snps_without_var]
        
}

# which were removed 
fwrite(tibble(snps_removed = names(snps_without_var), snps_removed_index = snps_without_var), 
       paste0(output_folder, out_name, "_snps_removed.txt"), nThread = 1, sep = " ")

# svd
svd_gen_mat <- svd(gen_dat_scaled)
S_mat <- diag(svd_gen_mat$d)

# write_ids
fwrite(tibble(id = ids), paste0(output_folder, out_name, "_ids.txt"), nThread = 4)
# save
fwrite(svd_gen_mat$u, paste0(output_folder, out_name, "_U.txt"), nThread = 4, col.names = FALSE, sep = " ")
fwrite(svd_gen_mat$v, paste0(output_folder, out_name, "_V.txt"), nThread = 4, col.names = FALSE, sep = " ")
fwrite(S_mat, paste0(output_folder, out_name, "_S.txt"), nThread = 4, col.names = FALSE, sep = " ")
        
# US mat
gen_US <- svd_gen_mat$u %*% S_mat
fwrite(gen_US, paste0(output_folder, out_name, "_US.txt"), nThread = 4, col.names = FALSE, sep = " ")
        




