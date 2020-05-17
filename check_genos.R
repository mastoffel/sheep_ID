# check susies suggestions

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_398k_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map <- as_tibble(full_sample$map)

sheep_map_oar <- read_delim("../sheep/data/SNP_chip/oar31_mapping/merged_sheep_geno_oar31.bim", 
                            "\t", col_names = FALSE) %>% 
                            filter(X1 %in% c(1:26)) %>% 
                            setNames(names(snps_map))

map_all <- snps_map %>% 
                left_join(sheep_map_oar, by = "snp.name") %>% 
                mutate(chr_diff = chromosome.x != chromosome.y) %>% 
                filter(chromosome.x == 1) %>% 
                mutate(gwas_hit =  as.factor(ifelse(position.x == 142390230, 1, 0)))

ggplot(map_all, aes(position.x, position.y)) +
        geom_point(aes(color = gwas_hit)) +
        geom_vline(aes(xintercept = 142390230))


gwas_res <- read_rds("output/gwas_res_397k.rds")
# put into df
gwas_full <- gwas_res %>%
        rename(snp.name = term) %>%
        filter(!str_detect(snp.name, "sd")) %>% 
        mutate(groups = ifelse(str_detect(snp.name, "roh"), "roh", "add")) %>% 
        mutate(snp.name = str_replace(snp.name, "roh_", "")) %>%
        left_join(snps_map) 

# extract roh
gwas_roh <- gwas_full %>% filter(groups == "roh") 
gwas_roh %>% arrange(desc(-log10(p.value)))





