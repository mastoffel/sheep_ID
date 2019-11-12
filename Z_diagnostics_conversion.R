# diagnostics and conversion

# missing genotypes
system("~/programs/plink --bfile data/sheep_geno_imputed --missing --out diagnostics/missing/miss_stat --sheep")
# allele frequencies
system("~/programs/plink --bfile data/sheep_geno_imputed --freq --out diagnostics/AF/freq_stat --sheep")

# convert to ped to add family id
system("~/programs/plink --bfile data/sheep_geno_imputed --recode --tab --out data/sheep_imp --sheep --no-fid")
# convert back to bed
system("~/programs/plink --file data/sheep_imp --make-bed --out data/sheep_imp --sheep")
