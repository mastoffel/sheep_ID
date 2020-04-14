library(tidyverse)
library(snpStats)
source("theme_simple.R")
library(viridis)
library(magrittr)
library(patchwork)
library(furrr)
chr_info <- read_delim("../sheep/data/sheep_genome/chromosome_info_ram.txt", "\t") %>% 
                .[-1, ] %>% 
                rename(chromosome = Part) %>% 
                mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
                mutate(chromosome = as.integer(chromosome)) %>% 
                filter(!is.na(chromosome))

#gwas_files <- list.files("output/gwas_full_pca_with_f", pattern = "*.rds", full.names = TRUE)
gwas_files <- list.files("output/gwas_new/second_run_age_nolamb/", pattern = "*.rds", full.names = TRUE)
#gwas_files <- list.files("output/gwas_survival_gaussian/", pattern = "*.rds", full.names = TRUE)

# extract results
all_gwas <- purrr::map(gwas_files, readRDS) %>% 
            purrr::flatten() %>% 
            purrr::flatten() %>% 
            .[seq(1,length(.),by=2)] #%>% 
           # remove snps that didnt work
           # purrr::compact()

# get roh/add pval
not_working <- map(all_gwas, is.null)
which(unlist(not_working))
all_gwas[[1364]]
gwas_res0 <- map(all_gwas, function(x) nrow(x) == 19) %>% unlist()
which(!gwas_res0)

all_gwas <- all_gwas[-which(unlist(not_working))]

plan(multiprocess, workers = 8)
gwas_res <- future_map_dfr(all_gwas, function(x) x %>% .[c(14,15), ] %>% 
                           dplyr::select(term, estimate, p.value))

# check whether snp roh didnt work in some models anywhere
gwas_res[which(str_detect(gwas_res$term, "sd")), ]

# remove those 
gwas_res <- gwas_res[-which(str_detect(gwas_res$term, "sd")), ]

# plink name
sheep_plink_name <- "output/plink_files/sheep_geno_imputed_ram_pruned"
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
        #filter(estimate>-0.1) %>% 
        # roh variable couldnt be estimated in some models, so the next
        # row was extracted from the tidy output
        filter(!str_detect(snp.name, "sd")) %>% 
        mutate(groups = ifelse(str_detect(snp.name, "roh"), "roh", "add")) %>% 
        mutate(snp.name = str_replace(snp.name, "roh_", "")) %>%
        left_join(snps_map) 

# qqplots
qqman::qq(gwas_full[gwas_full$groups == "add", ]$p.value)
qqman::qq(gwas_full[gwas_full$groups == "roh", ]$p.value)
qualityTools::qqPlot(gwas_full[gwas_full$groups == "roh", ]$p.value)
GWASTools::qqPlot(gwas_full[gwas_full$groups == "roh", ]$p.value)

# manhattan
## computing new x axis
gwas_roh <- gwas_full %>% filter(groups == "roh") 

# distribution of effect sizes
cols <- viridis(2)
# alternative
# cols <- c("#4393c3", "#2A3132")
options(scipen = 999)
p1 <- gwas_roh %>% 
  mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
  mutate(estimate = abs(estimate)) %>% 
  ggplot(aes(estimate, fill = direction)) +
  geom_histogram(bins = 1000) +
  theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
  theme(axis.line.y = element_blank()) +
  scale_fill_manual("Direction of\neffect on survival", values = cols) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1), labels = c(0.001, 0.01, 0.1, 1),
                limits = c(0.00001, 2)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Estimate (log-odds of survival)") +
  ylab("SNPs")
p1

p2 <- gwas_roh %>% 
  mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
  mutate(estimate = abs(estimate)) %>% 
  ggplot(aes(estimate, fill = direction)) +
  geom_histogram(bins = 10) +
  theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 7) +
  theme(axis.line.y = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(margin=margin(t=1)),
        axis.title.y = element_text(margin=margin(r=1)),
        axis.text.x=element_text(margin=margin(t=1)),
        axis.text.y=element_text(margin=margin(r=1))) +
  scale_fill_manual("Direction of\neffect on survival", values = cols) +
  scale_x_continuous(limits = c(0.7, 1.6)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Estimate") +
  ylab("SNPs")
p2

p3 <- gwas_roh %>% 
  mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
  ggplot(aes(p.value, fill = direction)) +
  geom_histogram(bins = 1000) +
  theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
  scale_fill_manual("Direction of\neffect on survival", values = cols) +
  theme(axis.line.y = element_blank()) +
  scale_x_continuous(limits = c(0, 1),  expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("SNPs") +
  xlab("p-value")
p3

p4 <- gwas_roh %>% 
  mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
  ggplot(aes(p.value, fill = direction)) +
  geom_histogram(bins = 30) +
  theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 7) +
  scale_fill_manual("Direction of\neffect on survival", values = cols) +
  theme(axis.line.y = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(margin=margin(t=1)),
        axis.title.y = element_text(margin=margin(r=1)),
        axis.text.x=element_text(margin=margin(t=1)),
        axis.text.y=element_text(margin=margin(r=1)))+
  scale_x_continuous(limits = c(0, 0.001),  expand = c(0,0),
                     breaks = c(0, 0.001)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("SNPs") +
  xlab("p-value") 
  
p4

layout <- c(
  area(t = 1, l = 1, b = 12, r = 16),
  area(t = 2, l = 2, b = 8, r = 7),
  area(t = 13, l = 1, b = 24, r = 16),
  area(t = 14, l = 2, b = 20, r = 7)
)
p_effsize <- p1 + p2 + p3 + p4 +
  plot_layout(design = layout, guides = 'collect')

p_effsize
ggsave("figs/eff_size_dist_roh.jpg", p_effsize, width = 4.5, height = 4.5)




# manhattan plots
chr_labels <- c(c(1:18),"","20","",  "22","", "24","", "26")
cols <- c("#336B87", "#2A3132")

# get cumsums
chr_info2 <- chr_info %>% 
        mutate(tot=cumsum(Length)-Length) %>% 
        dplyr::select(chromosome, tot)

gwas_p <- gwas_roh %>% 
        left_join(chr_info2) %>% 
        # Add a cumulative position of each SNP
        arrange(chromosome, position) %>%
        mutate(positive_cum = position + tot) 

axisdf <- gwas_p %>% group_by(chromosome) %>% 
                summarize(center = (max(positive_cum) + min(positive_cum)) / 2 )


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
gwas_plot <- gwas_p %>% left_join(hom_sum)

# check whether effects sizes and roh_count correlate
p_rohvsgwas <- ggplot(gwas_plot, aes(roh_count, estimate)) + 
      geom_point(shape = 21, alpha = 1, stroke = 0.1) +
      geom_smooth(se = FALSE, method = "lm") +
      theme_simple(axis_lines = TRUE, 
                   grid_lines = FALSE,
                   base_size = 8) + 
      facet_wrap(~chromosome, scales = "free") +
      theme(axis.title = element_text(size = 13)) +
      ylab("Estimate (log-odds)") +
      xlab("ROH count")
p_rohvsgwas
#ggsave("figs/roh_vs_gwas.jpg", width = 10, height = 6)

cols <- c("#4393c3", "#2A3132")
cols <- viridis(2)
#cols <- c(viridis(2)[1], "grey", viridis(2)[2], "lightgrey")
gwas_plot <- gwas_plot %>% 
  mutate(direction = ifelse(estimate < 0, "negative", "positive"))
  
pgwas <- ggplot(gwas_plot, aes(x=positive_cum, y=-log10(p.value))) +
        # Show all points
        #geom_point(aes(fill = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
        #           size = 2, shape = 21, alpha = 0.7, stroke = 0.02, color = "black") +
        geom_point(aes(fill = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
                   size = 2, shape = 21, alpha = 1, stroke = 0, color = "black") +
        geom_point(data = gwas_plot %>% filter(-log10(p.value) > -log10(0.05/31705)), # aes(fill = direction),  0.00001
                  fill=cols[[1]], size = 2, shape = 21, alpha = 0.5, stroke = 0.3, color = "black") +
       # gghighlight(-log10(p.value) > 2,
        #            unhighlighted_params = list(aes(fill = chromosome %%2 == 0))) +
        geom_hline(yintercept = -log10(0.05/31705), linetype="dashed", color = "grey") +
        scale_x_continuous(labels = chr_labels, breaks= axisdf$center ) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,8)) +
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
        #scale_fill_manual(values = c("#d8dee9", viridis(2),"#ECEFF4")) + 
        scale_fill_manual(values = c("#d8dee9","#ECEFF4")) + 
        #scale_fill_manual(values = c("#d8dee9","#4393c3","#2A3132","#ECEFF4")) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        theme(axis.text.x = element_text(size = 8),
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE) 

pgwas
#ggsave( "figs/survival_gwas_roh_pca.jpg",pgwas, height = 3, width = 15)

layout2 <- c(
  area(t = 1, l = 1, b = 12, r = 16),
  area(t = 2, l = 2, b = 8, r = 7),
  area(t = 1, l = 17, b = 12, r = 33),
  area(t = 2, l = 18, b = 8, r = 23),
  area(t = 13, l = 1, b = 34, r = 33)
)
p_effsize_h <- p1 + p2 + p3 + p4 + pgwas +
  plot_layout(design = layout2, guides = 'collect') +
  plot_annotation(tag_levels = 'A')
p_effsize_h 

ggsave(filename = "figs/gwas_plot_viridis.jpg",p_effsize_h, width = 9, height = 5)

cols <- c("#440154FF", "#FDE725FF")
# add effect size of significant snps
pgwas2 <- gwas_plot %>% 
  #filter(-log10(p.value) > -log10(0.05/28946)) %>% 
  filter(-log10(p.value) > -log10(0.05/31705)) %>% 
  mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
ggplot(aes(estimate, fill = direction)) +
  geom_histogram(bins = 7) +
  scale_fill_manual("Direction of\neffect on survival", values = cols) +
  theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 7) +
  #scale_x_continuous(breaks = c(-1.5, -1, -0.5, 0, 0.5), limits = c(-1.6, 0.7)) +
  scale_x_continuous(limits = c(-0.9, -0.55)) +
  theme(axis.line.y = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(margin=margin(t=1)),
        axis.title.y = element_text(margin=margin(r=1)),
        axis.text.x=element_text(margin=margin(t=1)),
        axis.text.y=element_text(margin=margin(r=1))) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("SNPs") +
  xlab("estimate") 
pgwas2 

layout3 <- c(
  area(t = 1, l = 1, b = 12, r = 16),
  area(t = 2, l = 2, b = 8, r = 7),
  area(t = 1, l = 17, b = 12, r = 33),
  area(t = 2, l = 18, b = 8, r = 23),
  area(t = 14, l = 1, b = 35, r = 33),
  area(t = 13, l = 12, b = 16, r = 17)
)
p_effsize_h2 <- p1 + p2 + p3 + p4 + pgwas + pgwas2 +
  plot_layout(design = layout3, guides = 'collect') +
  plot_annotation(tag_levels = 'A')
p_effsize_h2

ggsave(filename = "figs/gwas_plot_viridis_190k.jpg", p_effsize_h2, width = 9, height = 5)


# detailed look at significant regions
gwas_plot %>% arrange(p.value)

get_genome_region <- function(chr, pos) {
  gwas_tbl <- gwas_plot %>% filter(
    (chromosome == chr) & ((position > (pos - plusminus  * 1000000)) & (position < (pos +  plusminus * 1000000)))
  )
  out <- gwas_tbl %>% left_join(hom_sum) %>%
      mutate(log_p = -log10(p.value), estimate = abs(estimate)) %>%
      pivot_longer(names_to = "var", values_to = "vals",
                                cols = c("log_p", "roh_count", "estimate"))
  
}
plusminus <- 1.5 # Mb

top_snps <- gwas_plot %>% filter(-log10(p.value) > 5.6) %>% 
                      group_by(chromosome) %>% 
                      top_n(-1, wt = p.value)

out <- map2_df(top_snps$chromosome, top_snps$position, get_genome_region, .id = "snp")
df_plot <- out %>% 
            mutate(top_snp = ifelse(-log10(p.value) > 5.6, 1, 0)) %>% 
            mutate(top_snp = as.factor(top_snp)) %>% 
            mutate(pos_Mb = position/1000000) 

ggplot() +
  geom_point(data= df_plot %>% filter(top_snp == 0), aes(pos_Mb, y = vals), color = "lightgrey") +
  geom_point(data=df_plot %>% filter(top_snp == 1), aes(pos_Mb, y = vals), color = "blue") +
  facet_wrap(var~chromosome, scales = "free", nrow = 3,ncol = 6) +
 # scale_color_manual(values = viridis(2)) +
  #facet_wrap(~var, nrow = 3, scales = "free") + 
  theme_simple(grid_lines = FALSE, axis_lines = TRUE) + 
  theme(axis.line.x=element_line())

gwas_chr11 <- ggplot(gwas_plot_chr11, aes(x=positive_cum, y=-log10(p.value))) +
  geom_point()
  geom_point(size = 2, shape = 21, alpha = 1, stroke = 0, color = "black") +
  geom_point(data = gwas_plot_chr11 %>% filter(-log10(p.value) > -log10(0.05/31705)), # aes(fill = direction),  0.00001
             fill=cols[[1]], size = 2, shape = 21, alpha = 0.5, stroke = 0.3, color = "black") +
  # gghighlight(-log10(p.value) > 2,
  #            unhighlighted_params = list(aes(fill = chromosome %%2 == 0))) +
  geom_hline(yintercept = -log10(0.05/31705), linetype="dashed", color = "grey") +
  scale_x_continuous(labels = chr_labels, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,8)) +
  xlab("Chromosome") + 
  ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
  #scale_fill_manual(values = c("#d8dee9", viridis(2),"#ECEFF4")) + 
  scale_fill_manual(values = c("#d8dee9","#ECEFF4")) + 
  #scale_fill_manual(values = c("#d8dee9","#4393c3","#2A3132","#ECEFF4")) +
  theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
  theme(axis.text.x = element_text(size = 8),
        axis.ticks = element_line(size = 0.1)) +
  guides(fill=FALSE) 




pgwas_effs <- ggplot(gwas_plot, aes(x=positive_cum, y=exp(estimate))) +
  # Show all points
  geom_hline(yintercept = 0, size = 0.2) +
  geom_point(aes(fill = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
             size = 3, alpha = 0.5, shape = 21, stroke = 0.1) +
  scale_x_continuous(labels = chr_labels, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0,8)) +
  xlab("Chromosome") + 
  ylab("") + ## y label from qqman::qq
  scale_fill_manual(values = cols) + ##values = cols
  theme_simple(axis_lines = TRUE) +
  theme(axis.text.x = element_text(size = 10),
        axis.ticks = element_line(size = 0.1)) +

  guides(fill=FALSE) 
pgwas_effs 

pgwas / pgwas_effs

ggplot(gwas_plot, aes(estimate)) +
  geom_histogram(bins = 1000) +
  scale_y_sqrt()

plot(gwas_plot$p.value, abs(gwas_plot$estimate))
# save top snps
top_roh_snps <- gwas_full %>% 
        arrange(p.value) %>% 
        filter(groups == "roh") %>%  # filter(p.value < 0.05/28946) 
        filter(p.value < 0.00005) %>% 
        # only take top snp per peak
        mutate(positive_round = round(position/1000000)) %>% 
        arrange(chromosome, positive_round, p.value) %>% 
        group_by(chromosome, positive_round) %>% 
        top_n(-1, p.value)
      
gwas_full %>% 
        filter(snp.name %in% top_roh_snps$snp.name) %>% 
        write_delim("output/top_roh_snps_gwas_testset.txt")

gwas_roh %>% arrange(p.value) %>% filter(groups == "roh") %>% filter(p.value < 0.05/28946)



# effect size
gwas_plot %>% 
        #filter(chromosome == 20) %>% 
       # filter(p.value < 0.05/90000) %>% 
ggplot(aes(x=positive_cum, y=estimate)) + # inv.logit(estimate)
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
pgwas <- ggplot(gwas_plot, aes(x=positive_cum, y=-log10(p.value))) +
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
