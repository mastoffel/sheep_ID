library(tidyverse)
library(snpStats)
source("theme_clean.R")
library(viridis)
library(magrittr)
chr_info <- read_delim("../sheep/data/sheep_genome/chromosome_info_ram.txt", "\t") %>% 
                .[-1, ] %>% 
                rename(chromosome = Part) %>% 
                mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
                mutate(chromosome = as.integer(chromosome)) %>% 
                filter(!is.na(chromosome))

gwas_files <- list.files("output/gwas_full_pca_with_f", pattern = "*.rds", full.names = TRUE)
#gwas_files <- list.files("output/gwas_survival_gaussian/", pattern = "*.rds", full.names = TRUE)

# extract results
all_gwas <- purrr::map(gwas_files, readRDS) %>% 
            purrr::flatten() %>% 
            # only extract results
            .[seq(1,length(.),by=2)] %>% 
           # remove snps that didnt work
            purrr::compact()

# get roh/add pval
gwas_res <- map_df(all_gwas, function(x) x %>% .[c(14,15), ] %>% 
                           dplyr::select(term, estimate, p.value))

# check whether snp roh didnt work in some models anywhere
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
                        filter(groups == "roh") #%>% 
                       # filter(chromosome == 20)
                        #group_by(groups) %>% 
                        #arrange(chromosome, position) %>% 
                        #dplyr::mutate(tmp = 1, cumsum.tmp = cumsum(tmp))
## calculating x axis location for chromosome label
# med.dat <- gwas_roh %>% dplyr::group_by(groups, chromosome) %>% 
#                 dplyr::summarise(median.x = median(cumsum.tmp))
library(QQperm) # lambda roh = 1.602886, lambda add = 1.570673
exp_p <- runif(nrow(gwas_roh), min = 0, max = 1)
lambda <- estlambda2(gwas_roh$p.value, exp_p)$estimate
gwas_roh <- gwas_roh %>% mutate(p_cor = p.value/lambda)


chr_labels <- c(c(1:18),"","20","",  "22","", "24","", "26")
#chr_labels <- unique(gwas_roh$chromosome)
#chr_labels <- med.dat$chromosome

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


# roh info
#devtools::install_github('tavareshugo/windowscanr')
library(windowscanr)
library(data.table)
hom_sum <- fread("output/ROH/roh_nofilt_ram.hom.summary") %>% 
                rename(snp.name = SNP, roh_count = UNAFF) %>% 
                dplyr::select(snp.name, roh_count) 

quants <- quantile(hom_sum$roh_count, probs = c(0.1, 0.9))
hom_sum <- hom_sum %>% mutate(roh_prevalence = case_when(
                                roh_count < quants[1] ~ "<10% of inds",
                                roh_count > quants[2] ~ ">90% of inds",
                                TRUE ~ "in between")
                              )

gwas_plot <- gwas_p %>% 
                left_join(hom_sum)

cols <- c("#4393c3", "#2A3132")
pgwas <- ggplot(gwas_plot, aes(x=pos_cum, y=-log10(p.value))) +
        # Show all points
 geom_point(aes(fill = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
                   size = 3, alpha = 0.5, shape = 21, stroke = 0.1) +
        # geom_point(aes(fill = chromosome %%2 == 0,  shape = roh_prevalence), #  shape = roh_prevalence  #fill = chromosome %%2 == 0
        #            size = 2.5, alpha = 0.5, stroke = 0.2) +
        #geom_point(size = 2.5, alpha = 1, shape = 21, stroke = 0.1) + 
        #scale_color_manual(values = rep(cols, 26 )) +
        geom_hline(yintercept = -log10(0.05/28946), linetype="dashed", color = "grey") +
        # for top snps
        #geom_hline(yintercept = -log10(0.05/1000), linetype="dashed", color = "grey") +
        #scale_shape_manual(values = c(25,24,21)) +
   
        scale_x_continuous(labels = chr_labels, breaks= axisdf$center ) +

        scale_y_continuous(expand = c(0, 0), limits = c(0,8)) +
        # Add label using ggrepel to avoid overlapping
       # geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
        scale_fill_manual(values = cols) +##values = cols
        theme_clean() +
        theme(axis.text.x = element_text(size = 10),
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE) 
        #facet_wrap(groups~., nrow = 2)

pgwas

ggsave( "figs/survival_gwas_roh_pca_with_f_no_shape.jpg",pgwas, height = 3, width = 15)

# save top snps
top_roh_snps <- gwas_full %>% 
        arrange(p.value) %>% 
        filter(groups == "roh") %>%  # filter(p.value < 0.05/28946) 
        filter(p.value < 0.00005) %>% 
        # only take top snp per peak
        mutate(pos_round = round(position/1000000)) %>% 
        arrange(chromosome, pos_round, p.value) %>% 
        group_by(chromosome, pos_round) %>% 
        top_n(-1, p.value)
      
gwas_full %>% 
        filter(snp.name %in% top_roh_snps$snp.name) %>% 
        write_delim("output/top_roh_snps_gwas_testset.txt")

gwas_roh %>% arrange(p.value) %>% filter(groups == "roh") %>% filter(p.value < 0.05/28946)



# effect size
gwas_plot %>% 
        #filter(chromosome == 20) %>% 
       # filter(p.value < 0.05/90000) %>% 
ggplot(aes(x=pos_cum, y=estimate)) + # inv.logit(estimate)
        geom_point(aes(fill = chromosome %%2 == 0, shape = roh_prevalence),  #fill = chromosome %%2 == 0
                   size = 2.5, alpha = 0.5, stroke = 0.2) +
        #geom_hline(yintercept = -log10(0.05/28946), linetype="dashed", color = "grey") +
        scale_shape_manual(values = c(25,24,21)) +
        scale_x_continuous(labels = chr_labels, breaks= axisdf$center ) +
        #scale_x_continuous(breaks= axisdf$center ) +
        #scale_y_continuous(limits = c(-1,1)) +
        #scale_y_sqrt() + 
        xlab("Chromosome") + 
        ylab("estimate") + ## y label from qqman::qq
        scale_fill_manual(values = cols) +##values = cols
        theme_clean() +
        theme(axis.text.x = element_text(size = 10),
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE) 



# additive plot
pgwas <- ggplot(gwas_plot, aes(x=pos_cum, y=-log10(p.value))) +
        geom_point(aes(fill = chromosome %%2 == 0),  #fill = chromosome %%2 == 0
                   size = 2.5, alpha = 0.5, stroke = 0.2) +
        geom_hline(yintercept = -log10(0.05/28946), linetype="dashed", color = "grey") +
        scale_shape_manual(values = c(25,24,21)) +
        scale_x_continuous(labels = chr_labels, breaks= axisdf$center ) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,8)) +
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
        scale_fill_manual(values = cols) +##values = cols
        theme_clean() +
        theme(axis.text.x = element_text(size = 10),
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE) 

pgwas







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

gwas_full %>% 
        arrange(p.value) %>% 
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
