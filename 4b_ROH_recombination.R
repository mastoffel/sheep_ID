# ROH and recombination rate variation analysis
# Most of this is done in running windows across the genome.

library(data.table)
library(tidyverse)
source("theme_simple.R")
library(snpStats)
library(windowscanr)
library(patchwork)
library(gghalves)
library(lme4)
library(broom.mixed)
library(partR2)
#library(DHARMa)


# check ROH and recombination rate variation for shorter ROH
# system(paste0("/usr/local/bin/plink --bfile data/sheep_geno_imputed_oar_filt --sheep --out output/ROH/ROH_over_3Mb/roh ",
#               # "--keep output/ROH/ids_surv.txt ",
#               "--homozyg --homozyg-window-snp 50 --homozyg-snp 50 --homozyg-kb 3000 ",
#               "--homozyg-gap 300 --homozyg-density 200 --homozyg-window-missing 2 ",
#               "--homozyg-het 2 ",
#               "--homozyg-window-het 2"))


# ROH and recombination rate variation analysis --------------------------------

# 1) SNP statistics, heterozygosity
# plink name
sheep_plink_name <- "data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

snps_stats <- col.summary(full_sample$genotypes)
snps_stats <- snps_stats %>% as_tibble(rownames = "snp_name")
snps_stats_full <- full_sample$map %>% 
        as_tibble(rownames = "snp_name") %>% 
        mutate(KB = position / 1000) %>% 
        left_join(snps_stats, by = "snp_name") %>% 
        rename(het = P.AB)

# 2) roh prevalence 
snp_roh <- fread("output/ROH/roh.hom.summary") %>% 
        mutate(MB = BP / 1000000,
               KB = BP / 1000,
               index = 1:nrow(.))

# 3) linkage map from Johnston et al. (2020)
lmap <- read_delim("data/7_20200504_Full_Linkage_Map.txt", "\t") %>% 
        rename(SNP = SNP.Name)

# inner join subsets for shared snps
snp_df <- snp_roh %>% 
        inner_join(lmap, by = "SNP") %>% 
        mutate(KB = BP/1000)

# calculate running windows
# summed recombination fraction in 500 Kb non-overlapping running windows
running_rec <- winScan(x = snp_df,
                       groups = "CHR",
                       position = "KB",
                       values = c("cMdiff", "r"),
                       win_size = 500,
                       win_step = 500,
                       funs = c("sum"),
                       cores = 8)

running_rec <- running_rec %>% mutate(cM_Mb = (cMdiff_sum) / 0.5) 
running_rec %>% 
        ggplot(aes(r_sum, cM_Mb)) +
        geom_point() +
        facet_wrap(~CHR)

# mean ROH in 500 Kb non-overlapping running windows
running_roh <- winScan(x = snp_roh,
                       groups = "CHR",
                       position = "KB",
                       values = "UNAFF",
                       win_size = 500,
                       win_step = 500,
                       funs = c("mean"),
                       cores = 8)

# mean SNP heterozygosity in 500 Kb non-overlapping running windows.
running_het <- winScan(x = snps_stats_full,
                       groups = "chromosome",
                       position = "KB",
                       values = c("het"),
                       win_size = 500,
                       win_step = 500,
                       funs = c("mean"),
                       cores = 8)

running_het <- running_het %>% 
        rename(CHR = chromosome)

# recombination fraction vs cM position
ggplot(snp_df, aes(r, cMdiff)) +
        geom_point() +
        geom_smooth(method = "lm") +
        facet_wrap(~CHR)

# 4) islands and deserts
roh_isl_des <- read_delim(file = "output/roh_islands_deserts.txt", " ") %>% 
        rename(roh_n = UNAFF_n, roh_mean = UNAFF_mean)

# combine roh density in windows calculated from imputed data 
# and SNP heterozygosity in windows calculated from imputed data 
# and recombinatio fraction in windows calculated from 50k data
run_roh_rec_het <- running_roh %>% 
        as_tibble() %>% 
        rename(roh_n = UNAFF_n,
               roh_mean = UNAFF_mean) %>% 
        # filter windows with too few snps (lowest percentile)
        filter(roh_n > 35) %>% 
        mutate(prop_ROH = (roh_mean / 5952) * 100) %>% 
        left_join(running_rec, by = c("CHR", "win_start", "win_end", "win_mid")) %>% 
        left_join(running_het, by = c("CHR", "win_start", "win_end", "win_mid")) %>% 
        left_join(roh_isl_des) %>% 
        select(-prop_roh) %>% 
        mutate(extreme = ifelse(is.na(extreme), "average", extreme))


run_roh_rec_het

# create source data file for natcomms
library(writexl)
#writexl::write_xlsx(run_roh_rec_het, path = "output/source_data_fig2.xlsx")

roh_mod <- run_roh_rec_het %>% 
        # filter(row_number() %% 5 == 1) %>% 
        mutate(cM_Mb_std = as.numeric(scale(cM_Mb)),
               het_mean_std = as.numeric(scale(het_mean))) %>% 
        drop_na()

# model
mod <- lmer(prop_ROH ~ cM_Mb_std + het_mean_std + (1|CHR), data = roh_mod)
tidy(mod, conf.int = TRUE)
plan(multisession, workers = 8)
modR2 <- partR2(mod, partvars = c(" cM_Mb_std", "het_mean_std"), data = roh_mod,
                nboot = 1000, parallel = TRUE)
modrpt <- rptR::rptGaussian(prop_ROH ~ cM_Mb_std + het_mean_std + (1|CHR), data = roh_mod,
                            grname = "CHR")
summary(modR2)

sum_mod <- tidy(mod, conf.int = TRUE, conf.method = "boot")
sum_mod[4, "std.error"] <- (sum_mod[4, "conf.high"] - sum_mod[4, "conf.low"])/4
sum_mod[5, "std.error"] <- (sum_mod[5, "conf.high"] - sum_mod[5, "conf.low"])/4



library(gt)
# make table
sum_mod %>% 
        mutate(Term = c("Intercept", "Recombination rate (cM/Mb)", "Heterozygosity", "Chromosome", "Residual")) %>% 
        mutate(effect = c(rep("Fixed effects", 3), 
                          rep("Random effects (variances)", 2))) %>% 
        select(c(1,9,4,5,7,8,9)) %>% 
        setNames(c("effect", "Term", "Estimate", "Std.Error", "CI (2.5%)", "CI (97.5%)")) %>% 
        mutate(Standardization = c("", "(x-mean(x))/sd(x)", "(x-mean(x))/sd(x)", "", ""),
               Info = c("", "continuous", "continuous", "n = 26", ""),
               R2 = c("", "0.08, 95%CI [0.06, 0.11]", #0.42; 95%CI [0.40, 0.44]", 
                      "0.17, 95%CI [0.15, 0.19]", "", "")) %>% 
               
               # R2 = c("", "0.04, 95%CI [0.02, 0.07]", #0.42; 95%CI [0.40, 0.44]", 
               #        "0.38, 95%CI [0.36, 0.40]", "", "")) %>% 
        mutate(across(is.numeric, round, 3)) %>% 
        gt(
                rowname_col = "term",
                groupname_col = "effect"
        ) %>% 
        tab_style(
                style = cell_text( weight = "bold"),
                locations = cells_column_labels(columns = TRUE)
        ) %>% 
        fmt_markdown(columns = TRUE) %>% 
        gtsave("Rec_model_table_roh_over_3MB.png", path = "figs/tables/")


# only deserts/islands
roh_mod_extreme <- roh_mod %>% filter(extreme %in% c("island", "desert")) %>% 
        mutate(extreme_num = ifelse(extreme == "island", 0, 1))
mod2 <- lmer(prop_ROH ~ cM_Mb_std + het_mean_std + (1|CHR), data = roh_mod_extreme)
modR22 <- partR2(mod2, partvars = c("cM_Mb_std", "het_mean_std"), data = roh_mod_extreme,
                 nboot = 1000)

#all in one plot
coefs <- coef(lm(prop_ROH ~ cM_Mb, data = run_roh_rec_het))
library(viridis)
viridis(3)
p_roh_rec <- run_roh_rec_het %>% 
        ggplot(aes(cM_Mb, prop_ROH, fill = extreme)) +
        geom_point(shape = 21, stroke = 0.1, alpha = 0.7, size = 2,
                   color = "darkgrey") + #  fill = "#eceff4",
        ylab("% of sheep with ROH") +
        xlab("Recombination rate (cM/Mb)") +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_family = "Helvetica") +
        scale_fill_manual("ROH", values = c("#eceff4", viridis(3)[c(3, 1)]), 
                          breaks = c( "island", "desert")) +
        scale_x_continuous(breaks = seq(from=0, to=9, by = 2), limits = c(0, 8.8)) +
        geom_abline(intercept = coefs[1], slope = coefs[2], color = "#2e3440", size = 0.5) +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y=element_text(margin=margin(r=3)))

p_roh_rec
#ggsave("figs/roh_rec.jpg", plot = p_roh_rec, width = 5, height = 3.7)

coefs2 <- coef(lm(prop_ROH ~ het_mean, data = run_roh_rec_het))
p_roh_het <- run_roh_rec_het %>% 
        ggplot(aes(het_mean, prop_ROH,  fill = extreme)) +
        geom_point(shape = 21, stroke = 0.1, alpha = 0.7, size = 2,
                   color = "darkgrey") + #  fill = "#eceff4",
        ylab("% of sheep with ROH") +
        xlab("Heterozygosity") +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_family = "Helvetica") +
        scale_fill_manual("ROH", values = c("#eceff4", viridis(3)[c(3, 1)]), 
                          breaks = c( "island", "desert")) +
        scale_x_continuous(breaks = seq(from=0.1, to=0.4, by = 0.1), limits = c(0.025, 0.46)) +
        geom_abline(intercept = coefs2[1], slope = coefs2[2], color = "#2e3440", size = 0.5)  +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
p_roh_het

p_roh_rec / p_roh_het 

library(gghalves)
p_ext_rec <- run_roh_rec_het %>% 
        filter(extreme != "average") %>% 
        ggplot(aes(extreme, cM_Mb, fill = extreme)) +
        geom_hline(yintercept = mean(run_roh_rec_het$cM_Mb, na.rm = TRUE),
                   linetype = "dashed", color = "#4c566a", size = 0.5) +
        geom_half_point(side = "r", shape = 21, alpha = 0.7, stroke = 0.1, size = 2,
                        transformation_params = list(height = 0, width = 1.3, seed = 1),
                        color = "darkgrey") +
        geom_half_boxplot(side = "l", outlier.color = NA,
                          width = 0.6, lwd = 0.3, color = "black",
                          alpha = 0.8) + 
        scale_fill_manual(values = viridis(3)[c(3, 1)]) +
        scale_y_continuous(breaks = seq(from=0, to=9, by = 2), limits = c(0, 8.8)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_family = "Helvetica") +
        xlab("ROH region") +
        ylab("Recombination rate (cM/Mb)") +
        coord_flip() +
        theme(legend.position = "none",
              axis.title.y=element_text(margin=margin(r=3)))

p_ext_het <- run_roh_rec_het %>% 
        filter(extreme != "average") %>% 
        ggplot(aes(extreme, het_mean, fill = extreme)) +
        geom_hline(yintercept = mean(run_roh_rec_het$het_mean, na.rm = TRUE),
                   linetype = "dashed", color = "#4c566a", size = 0.5) +
        geom_half_point(side = "r", shape = 21, alpha = 0.7, stroke = 0.1, size = 2,
                        transformation_params = list(height = 0, width = 1.3, seed = 1),
                        color = "darkgrey") +
        geom_half_boxplot(side = "l", outlier.color = NA,
                          width = 0.6, lwd = 0.3, color = "black",
                          alpha = 0.8) + 
        scale_fill_manual(values = viridis(3)[c(3, 1)]) +
        scale_y_continuous(breaks = seq(from=0.1, to=0.4, by = 0.1), limits = c(0.025, 0.46)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_family = "Helvetica") +
        xlab("ROH region") +
        ylab("Heterozygosity") + 
        coord_flip() +
        theme(legend.position = "none",
              axis.title.y=element_blank())

p_full <- p_roh_rec + p_roh_het + p_ext_rec +  p_ext_het + 
        plot_layout(heights = c(1.41803398875, 1),
                    widths = c(0.5, 0.5)) +
        plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(face = "bold", vjust = 4),
              axis.title = element_text(size = 10),
              axis.text = element_text(colour = "black", size = 8))
p_full 
ggsave("figs/Fig2_roh_rec.jpg", width = 5.7, height = 4.3)
ggsave("figs/Fig2_roh_rec.pdf", width = 5.7, height = 4.3, device = cairo_pdf)


# supplementary plot
ggplot(run_roh_rec_het, aes(cM_Mb, prop_ROH)) +
        geom_point(shape = 21, fill = "#eceff4", stroke = 0.1, alpha = 1) +
        ylab("% of sheep with ROH") +
        xlab("Recombination rate (cM/Mb)") +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        geom_smooth(method = "lm", se = FALSE, color = "#5e81ac", size = 1) +
        scale_x_continuous( breaks = scales::pretty_breaks(2)) +
        scale_y_continuous( breaks = scales::pretty_breaks(4)) +
        theme(axis.title = element_text(size = 16)) +
        facet_wrap(~CHR, scales = "free") -> p_rec
p_rec
ggsave("figs/Sup_ROH_Rec.jpg", width = 8, height = 7)



# model df
running_snp <- run_roh_snp %>% 
        left_join(chr_sizes) %>% 
        mutate(size_MB = size_KB/1000) %>% 
        mutate(r_sum_std = as.numeric(scale(r_sum)),
               size_MB_std = as.numeric(scale(size_MB)),
               ROH_prop = roh_mean/5952)

# linear mixed model. r_sum is the recombination fraction in a given window
# size_MB is the chromosome size, CHR is the chromosome id
# predictors were standardised (scale()). 
mod <- lmer(ROH_prop ~ r_sum_std + size_MB_std + (1|CHR), data = running_snp)
summary(mod)

sum_mod <- tidy(mod, conf.int = TRUE, conf.method = "boot")
sum_mod[4, "std.error"] <- (sum_mod[4, "conf.high"] - sum_mod[4, "conf.low"])/4
sum_mod[5, "std.error"] <- (sum_mod[5, "conf.high"] - sum_mod[5, "conf.low"])/4

# make table
sum_mod %>% 
        mutate(Term = c("Intercept", "Recombination rate", "Chromosome size (Mb)", "Chromosome", "Residual")) %>% 
        mutate(effect = c(rep("Fixed effects", 3), 
                          rep("Random effects (variances)", 2))) %>% 
        select(c(1,9,4,5,7,8,9)) %>% 
        setNames(c("effect", "Term", "Estimate", "Std.Error", "CI (2.5%)", "CI (97.5%)")) %>% 
        mutate(Standardization = c("", "(x-mean(x))/sd(x)", "(x-mean(x))/sd(x)", "", ""),
               Info = c("", "continuous", "continuous", "n = 26", "")) %>% 
        mutate(across(is.numeric, round, 3)) %>% 
        gt(
                rowname_col = "term",
                groupname_col = "effect"
        ) %>% 
        tab_style(
                style = cell_text( weight = "bold"),
                locations = cells_column_labels(columns = TRUE)
        ) %>% 
        fmt_markdown(columns = TRUE) #%>% 
#gtsave("Rec_model_table.png", path = "figs/tables/")
#







