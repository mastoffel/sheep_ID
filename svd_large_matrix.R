library(data.table)
library(matrixStats)
library(RSpectra)
library(MASS)

# change here
out_name <- "roh_mat"
gen_mat <- fread("data/annual_survival_gwas_roh.txt", nThread = 4) %>% as.matrix()
#gen_mat <- fread("data/annual_survival_gwas_snps.txt", nThread = 4) %>% as.matrix()
#####

output_folder <- "output/"

#gen_mat <- matrix(sample(x = c(1, 0), size = 2000, replace = TRUE), nrow = 20)

# standardise
cm <- colMeans(gen_mat, na.rm = TRUE)
csd <- colSds(gen_mat, center = cm)
gen_mat_stand <- t( (t(gen_mat) - cm) / csd )

# svd
svd_gen_mat <- svd(gen_mat)
S_mat <- diag(svd_gen_mat$d)

# save
fwrite(svd_gen_mat$u, paste0(output_folder, out_name, "_U.txt"), nThread = 4, col.names = FALSE, sep = " ")
fwrite(svd_gen_mat$v, paste0(output_folder, out_name, "_V.txt"), nThread = 4, col.names = FALSE, sep = " ")
fwrite(S_mat, paste0(output_folder, out_name, "_S.txt"), nThread = 4, col.names = FALSE, sep = " ")

# US mat
gen_US <- svd_gen_mat$u %*% S_mat
fwrite(gen_US, paste0(output_folder, out_name, "_US.txt"), nThread = 4, col.names = FALSE, sep = " ")

