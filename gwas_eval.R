library(tidyverse)
library(snpStats)
source("theme_clean.R")

chr_info <- read_delim("../sheep/data/sheep_genome/chromosome_info_ram.txt", "\t") %>% 
                .[-1, ] %>% 
                rename(chromosome = Part) %>% 
                mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
                mutate(chromosome = as.integer(chromosome)) %>% 
                filter(!is.na(chromosome))

gwas_files <- list.files("output/gwas_full_roh_pca", pattern = "*.rds", full.names = TRUE)
# extract results
all_gwas <- purrr::map(gwas_files, readRDS) %>% 
            purrr::flatten() %>% 
            # only extract results
            .[seq(1,length(.),by=2)] %>% 
           # remove snps that didnt work
            purrr::compact()

# check errors
# all_gwas[seq(2,length(all_gwas),by=2)] %>% 
#       compact()

# get roh pval
gwas_res <- map_df(all_gwas, function(x) x %>% .[c(13,14), ] %>% 
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

full_sample$map %>% 
        group_by(chromosome) %>% 
        tally() %>% 
        write_delim("data/nsnps_pruned.txt")

# gwas_res_roh <- map_df(all_gwas, function(x) x %>% .[c(5), ] %>% dplyr::select(term,estimate, p.value))
# gwas_res_snp <- map_df(all_gwas, function(x) x %>% .[c(4), ] %>% dplyr::select(term,estimate, p.value))

# check how often roh var couldnt be estimated
# gwas_roh <- gwas_res %>% filter(str_detect(term, "roh"))
# gwas_roh[which(!(gwas_snp$term %in% full_sample$map$snp.name)), ]

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

qualityTools::qqPlot(gwas_full[gwas_full$groups == "add", ]$p.value)
GWASTools::qqPlot(gwas_full[gwas_full$groups == "roh", ]$p.value)
# manhattan
## computing new x axis
gwas_roh <- gwas_full %>% 
                        filter(groups == "roh") %>% 
                        filter(chromosome == 20)
                        #group_by(groups) %>% 
                        #arrange(chromosome, position) %>% 
                        #dplyr::mutate(tmp = 1, cumsum.tmp = cumsum(tmp))
## calculating x axis location for chromosome label
med.dat <- gwas_roh %>% dplyr::group_by(groups, chromosome) %>% 
                dplyr::summarise(median.x = median(cumsum.tmp))

chr_labels <- c(c(1:18),"","20","",  "22","", "24","", "26")
#chr_labels <- unique(gwas_roh$chromosome)
#chr_labels <- med.dat$chromosome
library(viridis)
library(magrittr)
cols <- c("#336B87", "#2A3132")

# get cumsums
chr_info %<>% 
        mutate(tot=cumsum(Length)-Length) %>% 
        dplyr::select(chromosome, tot)

gwas_p <- gwas_roh %>% 
        left_join(chr_info) %>% 
        # Add a cumulative position of each SNP
        arrange(chromosome, position) %>%
        mutate(pos_cum = position + tot) 

axisdf <- gwas_p %>% group_by(chromosome) %>% 
                summarize(center = (max(pos_cum) + min(pos_cum)) / 2 )


pgwas <- ggplot(gwas_p, aes(x=pos_cum, y=-log10(p.value))) +
        # Show all points
        geom_point(aes(color=as.factor(chromosome), fill = chromosome %%2 == 0),  
                   size = 2,color = "black", alpha = 0.7, shape = 21, stroke = 0.01) +
        #scale_color_manual(values = rep(cols, 26 )) +
        geom_hline(yintercept = -log10(0.05/28946), linetype="dashed", color = "grey") +
        # custom X axis:
        scale_x_continuous(labels = chr_labels, breaks= axisdf$center ) +
        #scale_x_continuous(breaks= axisdf$center ) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,8)) +
        # Add label using ggrepel to avoid overlapping
       # geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
        scale_fill_manual(values = cols) +## instead of colors, go for gray
        theme_clean() +
        theme(axis.text.x = element_text(size = 10),
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE) 
        #facet_wrap(groups~., nrow = 2)

pgwas
ggsave( "figs/survival_gwas_roh_pca_20.jpg",pgwas, height = 3, width = 12)

gwas_roh %>% arrange(p.value) %>% filter(groups == "roh") %>% filter(p.value < 0.05/28946)


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
        filter(groups == "roh") %>%  # filter(p.value < 0.05/28946) 
        .[1:100, ] %>% 
        write_delim("output/top_snps_gwas_pca.txt")


snps <- gwas_roh %>% arrange(p.value) %>% 
        filter(groups == "roh") %>% 
        filter(p.value < 0.05/28946) %>% 
        write_delim("output/top_snps_gwas_pca.txt")


# some other plots
gwas_roh_mod <- gwas_roh %>% filter(chromosome == 10)
ggplot(data = gwas_roh_mod) + 
        geom_point(aes(x = cumsum.tmp, y = -log10(p.value), fill = chromosome %%2 == 0),
                   size = 2.5,color = "black", alpha = 0.7, shape = 21, stroke = 0.01) + ## create alternate coloring
        geom_hline(yintercept = -log10(0.05/28946), linetype="dashed", color = "grey") + ## add horizontal line
       # scale_x_continuous(breaks = med.dat$median.x, labels = chr_labels) + ##  labels = chr_labels
        # guides(colour=FALSE) +  ## remove legend
        guides(fill=FALSE) +
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
        scale_fill_manual(values = cols) +## instead of colors, go for gray
        theme_clean() +
        theme(axis.text.x = element_text(size = 10),
              axis.ticks = element_line(size = 0.1))+
        facet_wrap(groups~., nrow = 2)
