library(tidyverse)
library(snpStats)
source("theme_simple.R")
library(viridis)
library(magrittr)
library(patchwork)
library(furrr)

chr_info <- read_delim("../sheep/data/sheep_genome/chromosome_info_oar31.txt", "\t") %>% 
                .[-1, ] %>% 
                rename(chromosome = Part) %>% 
                mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
                mutate(chromosome = as.integer(chromosome)) %>% 
                filter(!is.na(chromosome))

#gwas_files <- list.files("output/gwas_full_pca_with_f", pattern = "*.rds", full.names = TRUE)
gwas_files <- list.files("output/gwas_new/", pattern = "*.rds", full.names = TRUE) # oar31_roh_long/
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
all_gwas[[67507]]
gwas_res0 <- map(all_gwas, function(x) nrow(x) == 19) %>% unlist()
which(!gwas_res0)

# remove models which did not work
all_gwas <- all_gwas[-which(unlist(not_working))]

# extract additive and roh
plan(multiprocess, workers = 8)
gwas_res <- future_map_dfr(all_gwas, function(x) x %>% .[c(14,15), ] %>% 
                           dplyr::select(term, estimate, p.value))

# check whether snp roh didnt work in some models anywhere
gwas_res[which(str_detect(gwas_res$term, "sd")), ]

# remove those 
gwas_res <- gwas_res[-which(str_detect(gwas_res$term, "sd")), ]

saveRDS(gwas_res, file = "output/gwas_res_oar.rds")
gwas_res <- read_rds("output/gwas_res_oar.rds")
#gwas_res <- read_rds("output/gwas_res_oar.rds")

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_oar_filt"
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

# how many negative and positive effects?
roh_neg_pos <- gwas_full %>% filter(groups == "roh") %>% 
                    mutate(eff = ifelse(estimate <0, "neg", "pos")) %>% 
                    group_by(eff) %>% tally()

233920/(233920 + 163913)
binom.test(x= 233920, n = (233920 + 163913),
           p = 0.5, alternative = "greater",
           conf.level = 0.95)
formatC(0.00000000000000022, format = "e", digits = 2)

add_neg_pos <- gwas_full %>% filter(groups == "add") %>% 
  mutate(eff = ifelse(estimate <0, "neg", "pos")) %>% 
  group_by(eff) %>% tally()

binom.test(x= 203712, n = (203712 + 194122),
           p = 0.5, alternative = "two.sided",
           conf.level = 0.95)

# test whether negative estimates are more often larger estimates or
# smaller p-values
gwas_full %>% filter(groups == "roh") %>% 
  mutate(p_log = -log10(p.value)) %>% 
  mutate(p_log = as.numeric(scale(p_log)),
         estimate = as.numeric(scale(estimate))) %>% 
  mutate(eff = ifelse(estimate <0, "neg", "pos")) %>% 
  mutate(eff = factor(eff, levels = c("pos", "neg"))) -> gwas_roh_mod
gwas_roh_mod

mod1 <- glm(eff ~ p_log, data = gwas_roh_mod,
      family = "binomial")
tidy(mod1, conf.int = TRUE)
mod2 <- glm(eff ~ abs(estimate), data = gwas_roh_mod,
            family = "binomial")
tidy(mod2, conf.int = TRUE)

# manhattan
## computing new x axis
gwas_roh <- gwas_full %>% filter(groups == "roh") 
gwas_roh %>% arrange(desc(abs(estimate)))
# distribution of effect sizes
cols <- viridis(2)
# alternative
#cols <- c("#4393c3", "#2A3132")
options(scipen = 999)
  
p1 <- gwas_roh %>% 
  mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
  mutate(estimate = abs(estimate)) %>% 
  ggplot(aes(estimate, fill = direction)) +
  geom_histogram(bins = 500, position="identity") +
  theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
  theme(axis.line.y = element_blank()) +
  scale_fill_manual("Direction of\neffect on survival", values = cols) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 3), labels = c(0.001, 0.01, 0.1, 1, 3),
                limits = c(0.00005, 5.8)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Estimate (log-odds of survival)") +
  ylab("SNPs")
p1

p2 <- gwas_roh %>% 
  mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
  mutate(estimate = abs(estimate)) %>% 
  ggplot(aes(estimate, fill = direction)) +
  geom_histogram(bins = 30, position="identity") +
  theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 7) +
  theme(axis.line.y = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(margin=margin(t=1)),
        axis.title.y = element_text(margin=margin(r=1)),
        axis.text.x=element_text(margin=margin(t=1)),
        axis.text.y=element_text(margin=margin(r=1))) +
  scale_fill_manual("Direction of\neffect on survival", values = cols) +
  scale_x_continuous(limits = c(0.7, 3)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Estimate") +
  ylab("SNPs")
p2

p3 <- gwas_roh %>% 
  mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
  ggplot(aes(p.value, fill = direction)) +
  geom_histogram(bins = 500, position="identity") +
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
  geom_histogram(bins = 30,  position="identity") +
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
chr_labels <- c(c(1:10),"","12","", "14","", "16","", "18", "", "20","", "22","", "24","", "26")
chr_labels_full <- as.character(1:26)
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
hom_sum <- fread("output/ROH/roh_ram.hom.summary") %>% 
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
        geom_hline(yintercept = -log10(0.05/39149), linetype="dashed", color = "grey") +
        geom_point(data = gwas_plot %>% filter(-log10(p.value) <= -log10(0.05/39149)),
                   aes(fill = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
                   size = 2, shape = 21, alpha = 1, stroke = 0, color = "black") +
        geom_point(data = gwas_plot %>% filter(-log10(p.value) > -log10(0.05/39149)), # aes(fill = direction),  0.00001
                  fill=alpha(cols[[1]], 0.5), size = 2.5, shape = 21, stroke = 0.1, color = "black") + # "#d8dee9"
       # gghighlight(-log10(p.value) > 2,
        #            unhighlighted_params = list(aes(fill = chromosome %%2 == 0))) +
       # scale_x_continuous(labels = chr_labels, breaks= axisdf$center) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,9), labels = as.character(0:8), breaks = 0:8) +
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
 
        scale_fill_manual(values = c("#d8dee9","#ECEFF4")) + 
       # scale_fill_manual(values = cols) +

        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        theme(axis.text.x = element_text(size = 8),
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE) 

pgwas
ggsave( "figs/survival_gwas_oar.jpg",pgwas, height = 3, width = 15)

# simple plot without subplots
p_gwas_simple <- (p1 + p3) / pgwas +
  plot_layout(guides = 'collect', height = c(1, 1.5)) +
  plot_annotation(tag_levels = 'A')
p_gwas_simple

ggsave(filename = "figs/gwas_plot_viridis_397k.jpg", 
       p_gwas_simple, width = 8, height = 4)



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

#ggsave(filename = "figs/gwas_plot_viridis.jpg",p_effsize_h, width = 9, height = 5)

cols <- c("#440154FF", "#FDE725FF")
# add effect size of significant snps
pgwas2 <- gwas_plot %>% 
  #filter(-log10(p.value) > -log10(0.05/28946)) %>% 
  filter(-log10(p.value) > -log10(0.05/39149)) %>% 
  mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
ggplot(aes(abs(estimate), fill = direction)) +
  geom_histogram(bins = 25,  position="identity") +
  scale_fill_manual("Direction of\neffect on survival", values = cols) +
  theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 7) +
  #scale_x_continuous(breaks = c(-1.5, -1, -0.5, 0, 0.5), limits = c(-1.6, 0.7)) +
  #scale_x_continuous(limits = c(-0.85, -0.65)) +
  theme(axis.line.y = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(margin=margin(t=1)),
        axis.title.y = element_text(margin=margin(r=1)),
        axis.text.x=element_text(margin=margin(t=1)),
        axis.text.y=element_text(margin=margin(r=1))) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("SNPs") +
  xlab("Estimate") 
pgwas2  

layout3 <- c(
  area(t = 1, l = 1, b = 12, r = 16),
  area(t = 2, l = 2, b = 8, r = 7),
  area(t = 1, l = 17, b = 12, r = 33),
  area(t = 2, l = 18, b = 8, r = 23),
  area(t = 15, l = 1, b = 35, r  = 34),
  area(t = 15, l = 7, b = 17, r = 16)
)
p_effsize_h2 <- p1 + p2 + p3 + p4 + pgwas + pgwas2 +
  plot_layout(design = layout3, guides = 'collect') +
  plot_annotation(tag_levels = 'A')
p_effsize_h2

ggsave(filename = "figs/gwas_plot_viridis_400k.jpg", p_effsize_h2, width = 9, height = 5)

 




