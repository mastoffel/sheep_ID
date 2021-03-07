# Evaluate and plot GWAS output
# when running with example data, load libraries and start from line 50 #

library(tidyverse)
library(snpStats)
source("theme_simple.R")
library(viridis)
library(magrittr)
library(patchwork)
library(furrr)
library(windowscanr)
library(data.table)
library(broom)
library(broom.mixed)
# convert GWAS results into tidy data.frame ------------------------------------
# GWAS results are models tidied by broom.mixed::tidy() and saved as .rds
# Results from script 6_alt_gwas_annual_survival.R
gwas_files <- list.files("output/gwas_bothA_sep", pattern = "*.rds", full.names = TRUE) # oar31_roh_long/
# extract results
all_gwas <- purrr::map(gwas_files, readRDS) %>% 
        purrr::flatten() %>% 
        purrr::flatten() %>% 
        .[seq(1,length(.),by=2)] #%>% 

# remove snps that didnt work
not_working <- map(all_gwas, is.null)
which(unlist(not_working))

# remove models which did not work
all_gwas <- all_gwas[-which(unlist(not_working))]

# bind together and filter ROH predictors
gwas_res <- bind_rows(all_gwas) %>% 
            #filter(str_detect(term, "^roh")) %>% 
            filter(!(term %in% c("(Intercept)", "sexM", "twin1", "age_std",
                                 "age_std2", "sd__(Intercept)", paste0("pc", 1:7)))) %>% 
            filter(!str_detect(term, "froh_no_chr")) %>% 
            select(term, estimate, p.value) %>% 
            mutate(state = case_when(
                    str_detect(term, "roh_0") ~ "roh_0",
                    str_detect(term, "roh_2") ~ "roh_2",
                    TRUE ~ "add"
            )) %>% 
            mutate(term = str_remove(term, pattern = "roh_[0-9]_")) %>% 
            rename(snp = term)
       
gwas_res %>% 
        filter(state %in% c("roh_0", "roh_2")) %>% 
        pivot_wider(names_from = state, values_from = c("estimate", "p.value")) %>% 
        mutate(interesting = ifelse((p.value_roh_0 < 0.000005) | (p.value_roh_2 < 0.000005), 1, 0)) %>% 
        filter(interesting == 1) %>% 
        #filter(estimate_1 < 0 & estimate_2 <0) %>% 
        #filter(estimate_1 > -10) %>% 
        ggplot(aes(estimate_roh_0, estimate_roh_2)) +
        geom_point() +
        geom_smooth(method = "lm")

saveRDS(gwas_res , file = "output/gwas_res_oar_roh_sep.rds")

# start here when running example data -----------------------------------------
# get SNP map
sheep_plink_name <- "data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# geno_sub <- as_tibble(as(full_sample$genotypes[, 100000:101000], Class = "numeric"),
#                      rownames = "id")

# get snp map and add summary statistics
snps_map <- full_sample$map %>% as_tibble()
snps_stats <- col.summary(full_sample$genotypes)
plot(snps_stats$P.AA, snps_stats$P.BB)
snps_stats <- snps_stats %>% as_tibble(rownames = "snp.name")
snps_map <- snps_map %>% 
        left_join(snps_stats, by = "snp.name")
table(full_sample$map$chromosome, useNA = "always")

# chromosome info from assembly
chr_info <- read_delim("data/chromosome_info_oar31.txt", "\t") %>% 
        .[-1, ] %>% 
        rename(chromosome = Part) %>% 
        mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
        mutate(chromosome = as.integer(chromosome)) %>% 
        filter(!is.na(chromosome))

# load gwas results 
#gwas_res <- read_rds("output/gwas_res_oar_long.rds")

# put into df
gwas_full <- gwas_res %>%
        rename(snp.name = snp) %>%
        left_join(snps_map) 

# check which alleles
gwas_full %>% 
        filter(p.value < 0.05/(39149*2)) %>% 
        mutate(direction = ifelse(estimate < 0, "neg", "pos")) %>% 
        mutate(AF = ifelse(allele == 2, 1-MAF, MAF)) %>% 
       # filter(estimate < 0) %>% 
        ggplot(aes(as.factor(direction), AF)) +
                geom_boxplot() +
        geom_jitter(size = 3) 
        
# how many negative and positive effects? --------------------------------------
# roh effects
roh_neg_pos <- gwas_full %>%
        filter(state != "add") %>% 
        mutate(eff = ifelse(estimate <0, "neg", "pos")) %>% 
        group_by(eff) %>% tally()

# binomial test
binom.test(x= 465227, n = (465227 + 354279),
           p = 0.5, alternative = "greater",
           conf.level = 0.95)
formatC(0.00000000000000022, format = "e", digits = 2)

# add effects
add_neg_pos <- gwas_full %>%
        filter(state == "add") %>% 
        mutate(eff = ifelse(estimate <0, "neg", "pos")) %>% 
        group_by(eff) %>% tally()

binom.test(x= 207013, n = (207013 + 210346),
           p = 0.5, alternative = "two.sided",
           conf.level = 0.95)

# 39184

# test whether negative estimates are more often larger estimates or
# smaller p-values
gwas_full %>%
        filter(state != "add") %>% 
        mutate(p_log = -log10(p.value)) %>% 
        mutate(p_log = as.numeric(scale(p_log)),
               estimate = as.numeric(scale(estimate))) %>% 
        mutate(eff = ifelse(estimate <0, "neg", "pos")) %>% 
        mutate(eff = factor(eff, levels = c("pos", "neg"))) -> gwas_roh_mod
gwas_roh_mod


mod1 <- glm(eff ~ p.value , data = gwas_roh_mod,
            family = "binomial")
tidy(mod1, conf.int = TRUE)

mod2 <- glm(eff ~ abs(estimate), data = gwas_roh_mod,
            family = "binomial")
tidy(mod2, conf.int = TRUE)

mod1 <- glm(eff ~ p.value + abs(estimate) , data = gwas_roh_mod,
            family = "binomial")
tidy(mod1, conf.int = TRUE)

# manhattan plot ---------------------------------------------------------------
## computing new x axis
gwas_roh <- gwas_full %>% filter(state != "add")
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
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_family = "Helvetica") +
        theme(axis.line.y = element_blank()) +
        scale_fill_manual("Direction of\neffect on survival", values = cols) +
        #scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1), labels = c(0.001, 0.01, 0.1, 1),
        #              limits = c(0.00005, 1.1)) +
        scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10), limits = c(0.00005, 40),
                      labels = c(0.001, 0.01, 0.1, 1, 10)) +
        scale_y_continuous(expand = c(0,0)) +
        xlab("Estimate (log-odds of survival)") +
        ylab("SNPs")
p1

# p2 <- gwas_roh %>% 
#   mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
#   mutate(estimate = abs(estimate)) %>% 
#   ggplot(aes(estimate, fill = direction)) +
#   geom_histogram(bins = 30, position="identity") +
#   theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 7) +
#   theme(axis.line.y = element_blank(),
#         legend.position = "none",
#         axis.title.x = element_text(margin=margin(t=1)),
#         axis.title.y = element_text(margin=margin(r=1)),
#         axis.text.x=element_text(margin=margin(t=1)),
#         axis.text.y=element_text(margin=margin(r=1))) +
#   scale_fill_manual("Direction of\neffect on survival", values = cols) +
#   scale_x_continuous(limits = c(0.7, 3)) +
#   scale_y_continuous(expand = c(0,0)) +
#   xlab("Estimate") +
#   ylab("SNPs")

p3 <- gwas_roh %>% 
        mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
        ggplot(aes(p.value, fill = direction)) +
        geom_histogram(bins = 400, position="identity") +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_family = "Helvetica") +
        scale_fill_manual("Direction of\neffect on survival", values = cols) +
        theme(axis.line.y = element_blank()) +
        scale_x_continuous(limits = c(0, 1),  expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        ylab("SNPs") +
        xlab("p-value")
p3

# p4 <- gwas_roh %>% 
#   mutate(direction = ifelse(estimate < 0, "negative", "positive")) %>% 
#   ggplot(aes(p.value, fill = direction)) +
#   geom_histogram(bins = 30,  position="identity") +
#   theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 7) +
#   scale_fill_manual("Direction of\neffect on survival", values = cols) +
#   theme(axis.line.y = element_blank(),
#         legend.position = "none",
#         axis.title.x = element_text(margin=margin(t=1)),
#         axis.title.y = element_text(margin=margin(r=1)),
#         axis.text.x=element_text(margin=margin(t=1)),
#         axis.text.y=element_text(margin=margin(r=1)))+
#   scale_x_continuous(limits = c(0, 0.001),  expand = c(0,0),
#                      breaks = c(0, 0.001)) +
#   scale_y_continuous(expand = c(0,0)) +
#   ylab("SNPs") +
#   xlab("p-value") 
#   
# p4
# 
# layout <- c(
#   area(t = 1, l = 1, b = 12, r = 16),
#   area(t = 2, l = 2, b = 8, r = 7),
#   area(t = 13, l = 1, b = 24, r = 16),
#   area(t = 14, l = 2, b = 20, r = 7)
# )
# p_effsize <- p1 + p2 + p3 + p4 +
#   plot_layout(design = layout, guides = 'collect')
# 
# p_effsize
# ggsave("figs/eff_size_dist_roh.jpg", p_effsize, width = 4.5, height = 4.5)



# manhattan plots
#chr_labels <- c(c(1:18),"","20","",  "22","", "24","", "26")
chr_labels <- c(c(1:10),"","12","", "14","", "16","", "18", "", "20","", "22","", "24","", "26")
chr_labels_full <- as.character(1:26)
cols <- c("#336B87", "#2A3132")

# get cumsums
chr_info2 <- chr_info %>% 
        mutate(tot=cumsum(Length)-Length + cumsum(rep(35e6, 26))) %>% 
        dplyr::select(chromosome, tot)

gwas_plot <- gwas_roh %>% 
        left_join(chr_info2) %>% 
        # Add a cumulative position of each SNP
        arrange(chromosome, position) %>%
        mutate(positive_cum = position + tot) 
        

axisdf <- gwas_plot %>% group_by(chromosome) %>% 
        summarise(center = (max(positive_cum) + min(positive_cum)) / 2 )

#ggsave("figs/roh_vs_gwas.jpg", width = 10, height = 6)

#cols <- c("#4393c3", "#2A3132")
cols <- viridis(2)
gwas_plot_roh <- gwas_plot %>% 
        mutate(direction = ifelse(estimate < 0, "negative", "positive"))
# 39149
# 
eff_tests <- 2*39149
gwas_plot_roh_sub <- gwas_plot_roh #%>% 
        #sample_n(100000)

pgwas <- ggplot(gwas_plot_roh_sub, aes(x=positive_cum, y=-log10(p.value))) +
        geom_hline(yintercept = -log10(0.05/(eff_tests)), linetype="dashed", color = "grey") +
        #geom_point(data = gwas_plot_roh %>% filter(-log10(p.value) <= -log10(0.05/(eff_tests))),
        #           aes(fill = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
        #           size = 1.2, shape = 21, alpha = 1, stroke = 0, color = "transparent") +
        geom_point(data = gwas_plot_roh_sub %>% filter(-log10(p.value) <= -log10(0.05/(eff_tests))),
                   aes(color = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
                   size = 0.8) +
        # geom_point(data = gwas_plot %>% filter(-log10(p.value) > -log10(0.05/(39149*2))), # aes(fill = direction),  0.00001
        #            fill=alpha(cols[[1]], 0.8), size = 2.5, shape = 21, stroke = 0.1, color = "black") + # "#d8dee9
        geom_point(data = gwas_plot_roh_sub %>% filter(-log10(p.value) > -log10(0.05/(eff_tests))), aes(fill = direction), # aes(fill = direction),  0.00001
                   size = 2, shape = 21, stroke = 0.1, color = "black") + # "#d8dee9"
        scale_x_continuous(labels = chr_labels, breaks= axisdf$center) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,9), labels = as.character(0:8), breaks = 0:8) +
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
       # scale_fill_manual(values = c("#ECEFF4", cols[[1]], cols[[2]],"#d8dee9")) + # #dbe1eb #d1d8e5  "#ECEFF4" #d8dee9
        scale_fill_manual(values = c(cols[[1]], cols[[2]])) +
        scale_color_manual(values = c("#ECEFF4","#d8dee9")) + # #dbe1eb #d1d8e5  "#ECEFF4" #d8dee9
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_family = "Helvetica") +
        theme(axis.text.x = element_text(size = 8),
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE, color = FALSE) 
pgwas 


# pgwas <- ggplot(gwas_plot, aes(x=positive_cum, y=-log10(p.value))) +
#         geom_hline(yintercept = -log10(0.05/(39149*2)), linetype="dashed", color = "grey") +
#         geom_point(data = gwas_plot %>% filter(-log10(p.value) <= -log10(0.05/(39184*2))),
#                    aes(fill = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
#                    size = 2, shape = 21, alpha = 1, stroke = 0, color = "black") +
#         geom_point(data = gwas_plot %>% filter(-log10(p.value) > -log10(0.05/(39184*2))), # aes(fill = direction),  0.00001
#                    fill=alpha(cols[[1]], 0.8), size = 2.5, shape = 21, stroke = 0.1, color = "black") + # "#d8dee9"
#         scale_x_continuous(labels = chr_labels, breaks= axisdf$center) +
#         scale_y_continuous(expand = c(0, 0), limits = c(0,9), labels = as.character(0:8), breaks = 0:8) +
#         xlab("Chromosome") + 
#         ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
#         scale_fill_manual(values = c("#d8dee9","#ECEFF4")) + 
#         theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
#         theme(axis.text.x = element_text(size = 8),
#               axis.ticks = element_line(size = 0.1)) +
#         guides(fill=FALSE) 

#ggsave( "figs/survival_gwas_oar_shortnew.jpg",pgwas, height = 3, width = 15)

# simple plot without subplots
p_gwas_simple <- (p1 + p3) / pgwas +
        plot_layout(guides = 'collect', height = c(1, 1.5)) +
        plot_annotation(tag_levels = 'a') &
        theme(plot.tag = element_text(face = "bold"))
p_gwas_simple

ggsave(filename = "figs/gwas_plot_viridis_sep.jpg", 
       p_gwas_simple, width = 8, height = 3.7)
ggsave(filename = "figs/gwas_plot_viridis_sep.pdf", 
       p_gwas_simple, width = 8, height = 3.7)




# additive ---------------------------------------------------------------------
gwas_add <- gwas_full %>% filter(state == "add")
# manhattan plots
#chr_labels <- c(c(1:18),"","20","",  "22","", "24","", "26")
chr_labels <- c(c(1:10),"","12","", "14","", "16","", "18", "", "20","", "22","", "24","", "26")
chr_labels_full <- as.character(1:26)
cols <- c("#336B87", "#2A3132")

# get cumsums
chr_info2 <- chr_info %>% 
        mutate(tot=cumsum(Length)-Length) %>% 
        dplyr::select(chromosome, tot)

gwas_plot <- gwas_add %>% 
        left_join(chr_info2) %>% 
        # Add a cumulative position of each SNP
        arrange(chromosome, position) %>%
        mutate(positive_cum = position + tot) 

axisdf <- gwas_plot %>% group_by(chromosome) %>% 
        summarize(center = (max(positive_cum) + min(positive_cum)) / 2 )

#ggsave("figs/roh_vs_gwas.jpg", width = 10, height = 6)

#cols <- c("#4393c3", "#2A3132")
cols <- viridis(2)
gwas_plot_add <- gwas_plot %>% 
        mutate(direction = ifelse(estimate < 0, "negative", "positive"))
# 39149
# 
eff_tests <- 2*39149
pgwas_add <- ggplot(gwas_plot_add, aes(x=positive_cum, y=-log10(p.value))) +
        geom_hline(yintercept = -log10(0.05/(eff_tests)), linetype="dashed", color = "grey") +
        geom_point(data = gwas_plot_add %>% filter(-log10(p.value) <= -log10(0.05/(eff_tests))),
                   aes(fill = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
                   size = 2, shape = 21, alpha = 1, stroke = 0, color = "black") +
        # geom_point(data = gwas_plot %>% filter(-log10(p.value) > -log10(0.05/(39149*2))), # aes(fill = direction),  0.00001
        #            fill=alpha(cols[[1]], 0.8), size = 2.5, shape = 21, stroke = 0.1, color = "black") + # "#d8dee9
        geom_point(data = gwas_plot_roh %>% filter(-log10(p.value) > -log10(0.05/(eff_tests))), aes(fill = direction), # aes(fill = direction),  0.00001
                   size = 2.5, shape = 21, stroke = 0.1, color = "black") + # "#d8dee9"
        scale_x_continuous(labels = chr_labels, breaks= axisdf$center) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,9), labels = as.character(0:8), breaks = 0:8) +
        xlab("Chromosome") + 
        ylab(expression(-log[10](italic(p)))) + ## y label from qqman::qq
        scale_fill_manual(values = c("#ECEFF4", cols[[1]], cols[[2]],"#d8dee9")) + 
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        theme(axis.text.x = element_text(size = 8),
              axis.ticks = element_line(size = 0.1)) +
        guides(fill=FALSE) 
pgwas_add



