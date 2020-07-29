# Patterns of ROH analysis
# This script contains the analyses for the first part of the paper,
# exploring ROH variation across the genome, ROH islands and deserts,
# and the association between ROH prevalence and recombination rate variation.

library(data.table)
library(tidyverse)
source("theme_simple.R")
library(windowscanr)
library(cowplot)
library(gt)
library(grid)
library(ggplotify)
library(patchwork)
library(viridis)
library(scales)
library(ggridges)
library(gt)
library(lme4)
library(broom)
options(scipen=999)

# Chr lengths
chr_data <- read_delim("data/chromosome_info_oar31.txt", delim = "\t") %>% 
        rename(size_BP = Length,
               CHR = Part) %>% 
        mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
        .[2:27, ] %>% 
        summarise(sum_KB = sum(size_KB)) %>% 
        as.numeric()

#~~ ROH for survival data subset
file_path <- "output/ROH/roh.hom"
roh_lengths <- fread(file_path)

# max ROH length
roh_lengths[which.max(roh_lengths$KB), ]

# ROH overview -----------------------------------------------------------------
# descriptive ROH statistics
num_roh_per_ind <- roh_lengths %>% group_by(IID) %>% tally() 
summary(num_roh_per_ind$n)
sd(num_roh_per_ind$n)

# inbreeding coefficients
froh <- roh_lengths %>%
        dplyr::group_by(IID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/autosomal_genome_size) %>% 
        mutate(FROH_cent = FROH - mean(FROH))
mean(froh$FROH)
range(froh$FROH)

# longest ROH
roh_lengths %>% arrange(-KB)

# longest ROH proportional to chr size
chr_sizes <- chr_data %>% .[-1, ] %>% 
        mutate(CHR = str_replace(CHR, "Chromosome ", "")) %>% 
        mutate(CHR = as.integer(CHR))

roh_lengths %>% 
        arrange(desc(KB)) %>% 
        left_join(chr_sizes, by = "CHR") %>% 
        mutate(prop_chr = KB / size_KB) %>% 
        arrange(desc(prop_chr))
#plot(froh$KBSUM, num_roh_per_ind$n)

# ROH length and abundance in the 1% least and most inbred individuals 
num_roh_per_ind %>% 
        left_join(froh) %>% 
       # top_frac(-0.01, FROH) %>%    # top 1% least inbred individuals
        top_frac(0.01, FROH) %>%      # top 1% most inbred individuals
        summarise(mean(n), mean(KBAVG))

# Supplementary Figure FROH / ROH across individuals ---------------------------
p_froh <- ggplot(froh, aes(FROH)) +
        geom_histogram(bins = 100,  color = "white", size = 0.1, position = "identity",
                       alpha = 1, fill = "#4c566a") +
        ylab("individuals") +
        xlab(expression(inbreeding~coefficient~F[ROH])) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12, 
                     base_family = "Lato") 
p_froh

p_roh <- ggplot(num_roh_per_ind, aes(n)) +
        #geom_histogram(binwidth = 1,  fill = "#E5E9F0", color = "black",  size = 0.1) +
        geom_histogram(binwidth = 1, color = "white", size = 0.1, position = "identity",
                       alpha = 1, fill = "#4c566a") +
        ylab("individuals") +
        xlab("ROH per genome") +
        scale_y_continuous(expand = c(0, 0)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12, 
                     base_family = "Lato") 
p_roh

p_roh_dist <- p_froh + p_roh + plot_annotation(tag_levels = 'A')
p_roh_dist
ggsave("figs/Sup_ROH_dist.jpg", p_roh_dist, width = 7, height = 2.5)

# Supplementary Figure HD vs. imputed individuals ------------------------------
# this plot does not really work with example data, so please skip
# HD inds 
hd_inds <- read_delim("../sheep/data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram.fam", 
                      " ", col_names = FALSE)[[2]]
froh_plot <- froh %>% mutate(hd = ifelse(IID %in% hd_inds, "HD", "imputed")) %>% 
        mutate(hd = as.factor(hd)) %>% mutate(hd = relevel(hd, "imputed"))

# medians
froh_plot %>% group_by(hd) %>% summarise(median(FROH))

#
froh_imp_vs_hd <- ggplot(froh_plot, aes(x=FROH, fill = relevel(hd, "imputed"))) +
        geom_histogram(bins = 500,  color = "black", size = 0.1, position = "identity",
                       alpha = 0.8) +
        ylab("individuals") +
        xlab(expression(inbreeding~coefficient~F[ROH])) +
        scale_fill_manual(name = "Genotypes", values = c("#E5E9F0", "black"),
                          labels = c("partially imputed", "high-density")) +
       # scale_y_log10(expand = c(0, 0)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12,  # "#E5E9F0"
                     base_family = "Lato") +
        theme(legend.position="right")
froh_imp_vs_hd
#ggsave("figs/Sup_imp_vs_hd.jpg", froh_imp_vs_hd, width = 6, height = 2.7)

# alternative plot
p_imp <- froh_plot %>% 
        filter(hd == "imputed") %>% 
        ggplot(aes(x=FROH)) +
        geom_histogram(bins = 100,  color = "white", size = 0.1, position = "identity",
                       alpha = 1, fill = "#4c566a") + # #E5E9F0
        ylab("individuals") +
        scale_x_continuous(limits = c(0.18, 0.53)) +
        xlab(expression(inbreeding~coefficient~F[ROH])) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12,  # "#E5E9F0"
                     base_family = "Lato") +
        geom_vline(aes(xintercept = 0.241), size = 0.8, linetype = "dashed", color = "#d08770") +
        theme(legend.position="right",
              axis.line = element_blank(),
              plot.title = element_text(size = 12, face = "bold")) +
        ggtitle("imputed genotypes")

p_hd <- froh_plot %>% 
        filter(hd == "HD") %>% 
        ggplot(aes(x=FROH)) +
        geom_histogram(bins = 100,  color = "white", size = 0.1, position = "identity",
                       alpha = 1, fill = "#4c566a") +
        ylab("individuals") +
        scale_x_continuous(limits = c(0.18, 0.53)) +
        xlab(expression(inbreeding~coefficient~F[ROH])) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12,  # "#E5E9F0"
                     base_family = "Lato") +
        geom_vline(aes(xintercept = 0.239), size = 0.8, linetype = "dashed", color = "#d08770") +
        theme(legend.position="right",
              axis.line = element_blank(),
              plot.title = element_text(size = 12, face = "bold")) +
        scale_y_continuous(expand = c(0, 0)) +
        ggtitle("high-density genotypes")

p_out <- p_hd/p_imp + plot_annotation(tag_levels = "A")
p_out 
ggsave("figs/Sup_imp_vs_hd.jpg", p_out, width = 5.5, height = 4.6)

#~~ ROH for some indiividuals --------------------------------------------------
# Figure 1A, ROH for most and least inbred individuals

all_roh <- roh_lengths %>% 
        group_by(IID) %>% 
        summarise(sum_roh = sum(KB)) %>% 
        ungroup() %>% 
        arrange(desc(sum_roh))

longest_roh <- all_roh %>% 
        top_n(7)
shortest_roh <- all_roh %>% 
        top_n(-7)
num_ind <- 14

extreme_roh <- rbind(longest_roh, shortest_roh)

df <- roh_lengths %>%
        mutate(POS1 = POS1 / 1e+6,
               POS2 = POS2 / 1e+6,
               MB = KB / 1000)

df <- df %>% filter(IID %in% extreme_roh$IID) %>% 
        mutate(IID = factor(IID, levels = extreme_roh$IID))

yax <- data.frame(IID = fct_inorder(levels(df$IID))) %>%
        mutate(yax = seq(from = 2,
                         to = 2*length(unique(df$IID)),
                         by = 2)) 

df <- left_join(df, yax, by = "IID")

shade <- df %>%
        group_by(CHR) %>%
        summarise(min = min(POS1), max = max(POS2)) %>%
        mutate(min = case_when(CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10 |
                                       CHR == 12 | CHR == 14 | CHR == 16 | CHR == 18 | CHR == 20 |
                                       CHR == 22 | CHR == 24 | CHR == 26 ~ 0,
                               TRUE ~ min)) %>%
        mutate(max = case_when(CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10 |
                                       CHR == 12 | CHR == 14 | CHR == 16 | CHR == 18 | CHR == 20 |
                                       CHR == 22 | CHR == 24 | CHR == 26 ~ 0,
                               TRUE ~ max))

col <- c("#4c566a", "#d8dee9")
chr_names <- as.character(1:26)
names(chr_names) <- as.character(1:26)
chr_names[c(11, 13, 15, 17, 19, 21, 23, 25)] <- ""

df %>% 
        filter(MB > 5) %>% 
        filter(CHR %in% 1:26) %>% 
        ggplot() +
        geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
                  alpha=0.5, fill = "#eceff4") + # "#f7f7f7" "#eceff4"
        geom_hline(data = yax, aes(yintercept = yax), color = "#d8dee9", size = 0.4) +
        geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.9, 
                      fill = as.factor(CHR)),  col = "black", size = 0.1, alpha = 1) + 
        scale_fill_manual(values = rep(col, 18)) + 
        scale_color_manual(values = rep(col, 18)) +
        scale_y_reverse(expand = c(0, 0)) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13) +
        facet_grid(~CHR,scales = 'free_x', space = 'free_x', switch = 'x',
                   labeller = as_labeller(chr_names)) +
        theme(#strip.placement = 'outside',
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                panel.spacing = unit(0, "lines"),
                plot.margin = margin(r = 0.5, b = 0.5, t = 0.5, unit = "cm"),
                axis.line.x = element_blank(),
                legend.position="none",
                axis.title.x = element_text(margin=margin(t=0)),
                axis.title.y = element_text(margin=margin(r=0)),
                axis.text.y = element_text(colour = "white"),
                axis.line.y = element_blank()) +
        coord_cartesian(clip = 'off') +
        xlab("Chromosome") +
        ylab("Individuals") -> ROH_per_ind
ROH_per_ind
#ggsave("figs/fig1b_roh_per_ind_5Mb.jpg", ROH_per_ind, width = 7, height = 3)

# needs to be a grob to display axis labels correctly
pg <- ggplotGrob(ROH_per_ind)

for(i in which(grepl("strip-b", pg$layout$name))){
        pg$grobs[[i]]$layout$clip <- "off"
}

ROH_per_ind_grob <- as.ggplot(pg)
#ggsave("figs/fig1b_roh_per_ind_5Mb.jpg", ROH_per_ind_grob, width = 6, height = 3.5)

#p1 / p_roh_classes + plot_layout(heights = c(1, 0.8))

#~~~ ROH density ---------------------------------------------------------------
hom_sum <- fread("output/ROH/roh.hom.summary")

hom_sum <- hom_sum %>%
        mutate(MB = BP / 1000000,
               KB = BP / 1000,
               index = 1:nrow(.))

running_roh <- winScan(x = hom_sum,
                       groups = "CHR",
                       position = "KB",
                       values = "UNAFF",
                       win_size = 500,
                       win_step = 500,
                       funs = c("mean"))

# ROH sharing 1
running_roh %>% 
        mutate(UNAFF_mean = UNAFF_mean/5952) %>% 
        filter(UNAFF_n > 0) -> running_roh_p

#cor(running_roh$UNAFF_mean, running_roh$UNAFF_n, use = "complete.obs")
#plot(running_roh$UNAFF_mean, running_roh$UNAFF_n)

# remove windows without snps
running_roh %>% 
        mutate(UNAFF_mean = UNAFF_mean/5952) %>% 
        filter(UNAFF_n > 0) -> running_roh_p

# check distribution of snps
#ggplot(running_roh, aes(win_start, UNAFF_n)) + 
#  geom_point() + geom_smooth(span = 0.1) +
#  facet_wrap(~CHR, scales = "free_x")

# ROH across the genome plot
fill_cols <- viridis(20, option = "D")
qn <- scales::rescale(quantile(running_roh_p$UNAFF_mean,
                               probs=seq(0, 1, length.out=length(fill_cols))))

p1 <- ggplot(running_roh_p, aes(x = win_start, y = 0.5, fill = UNAFF_mean)) + 
        geom_tile(color = "black", size = 0) +
        theme_simple(base_size = 13, grid_lines = FALSE) + 
        scale_y_continuous(expand = c(0,0))+
        scale_x_continuous(expand = c(0,0), 
                           breaks = seq(0, 300000, by = 50000),
                           labels = as.character(seq(0, 300, 50)))+
        ylab("Chromosome") +
        scale_fill_gradientn("% of sheep with ROH",
                             colors = rev(fill_cols), values = qn,
                             breaks = c(0.1,0.3, 0.5, 0.7, 0.9),
                             labels = c(10, 30, 50 , 70, 90)) +
        facet_grid(CHR~., switch="both") +
        xlab("Position in Mb") +

        theme(panel.spacing.y=unit(0.1, "lines"),
              axis.title.x = element_text(margin=margin(t=5)),
              axis.title.y = element_text(margin=margin(r=5)),
              axis.text.y = element_blank(),
              axis.text.x = element_text(color = "black"),
              axis.ticks.x = element_line(size = 0.3),
              #axis.title.y = element_blank(),
              plot.margin = margin(r = 0.5, b = 0.5, unit = "cm"),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = c(0.805,0.13),
              legend.direction = "horizontal",
              strip.text.y.left = element_text(size = 10, angle = 0),
            #  axis.line.x = element_line(size = 0.3)
              ) +
        guides(fill = guide_colourbar(title.position = "bottom" ,
                                      barwidth = 10.95, barheight = 0.5))
p1
p2 <- ggplot(running_roh_p, aes(UNAFF_mean, "test", fill = ..x..)) +
        geom_density_ridges_gradient(scale = 2.5, lwd = 0.1) +
        scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) +
        theme_void() +
        scale_fill_gradientn("Proportion of Sheep with ROH",
                             colors = rev(fill_cols), values = qn,
                             breaks = c(0.1,0.3, 0.5, 0.7, 0.9)) +
        theme(legend.position = "none",
              plot.margin = margin(1, 1, 1, 1, unit = "cm"))

p2
ggsave("figs/roh_genome_legend.jpg", p2, width = 5, height = 2.5)

# try simple combined plot
p_roh_comb_simple <- plot_grid(ROH_per_ind_grob, p1, nrow = 2, 
                               rel_heights = c(0.658, 1), label_size = 15, 
                               labels = c("A", "B"), align = "v")
p_roh_comb_simple
ggsave("figs/roh_patterns_simple.jpg", p_roh_comb_simple, width = 7, height = 6.5)

# combine legend density and plot in keynote or somewhere else for full plot.

# ROH ISLANDS AND DESERTS ------------------------------------------------------

#~~~ ROH density
hom_sum <- fread("output/ROH/roh.hom.summary") # ROH_surv_subset/
head(hom_sum)

hom_sum <- hom_sum %>%
        mutate(MB = BP / 1000000,
               KB = BP / 1000,
               index = 1:nrow(.))# %>% 

# overview over SNPs with most and least ROH
hom_sum %>% 
        filter(UNAFF > 0) %>% 
        mutate(prop_roh = UNAFF/5952) %>% 
        arrange(prop_roh) %>% 
        .[1:50, ]

hom_sum %>% 
        filter(UNAFF > 0) %>% 
        mutate(prop_roh = UNAFF/5952) %>% 
        arrange(desc(prop_roh)) %>% 
        .[1:50, ]

# count ROH in running windows of 500 Kb
# UNAFF_n is number of SNPs in a window
# UNAFF_mean is mean ROH prevalence in a window
running_roh <- winScan(x = hom_sum,
                       groups = "CHR",
                       position = "KB",
                       values = "UNAFF",
                       win_size = 500,
                       win_step = 500,
                       funs = c("mean"),
                       cores = 8)
head(running_roh)
plot(running_roh$UNAFF_mean, running_roh$UNAFF_n)
cor(running_roh$UNAFF_mean, running_roh$UNAFF_n, use = "complete")

# check dist
hist(running_roh$UNAFF_mean, xlim = c(0,6000), breaks = 1000)


# roh desert
# filter windows with too few snps
quantile(running_roh$UNAFF_n, probs = c(0.01))

roh_deserts <- running_roh %>% 
        filter(UNAFF_n > 35) %>% 
        as_tibble() %>% 
        mutate(prop_roh = UNAFF_mean/5952) %>%  # 7691 5952
        arrange(prop_roh) %>% 
        # top 0.5% of windows
        .[1:24, ]
mean(roh_deserts$prop_roh)

# make table for supplementary
roh_deserts %>% 
        mutate(win_start = round(win_start/1000, 2), win_end = round(win_end/1000, 2), 
               prop_roh = round(prop_roh * 100, 2)) %>% 
        dplyr::select(CHR, win_start, win_end, prop_roh, UNAFF_n) %>%
        setNames(c("Chromosome", "WinStart", "WinEnd", "% of individuals with ROH", "N (SNPs)")) %>% 
        gt() %>% 
        tab_header(
                title = "Top 0.5% ROH deserts",
                subtitle = "ROH density measured in 500Kb running windows"
        ) #%>% 
       # gtsave(filename = "figs/tables/roh_desert.png")

roh_islands <- running_roh %>% 
        filter(UNAFF_n > 35) %>% 
        mutate(prop_roh = UNAFF_mean/5952) %>%  # 7691 5952
        arrange(desc(prop_roh)) %>% 
        .[1:24, ]

mean(roh_islands$prop_roh)

# make table for supplementary
roh_islands %>% 
        mutate(win_start = round(win_start/1000, 2), win_end = round(win_end/1000, 2), 
               prop_roh = round(prop_roh * 100, 2)) %>% 
        dplyr::select(CHR, win_start, win_end, prop_roh, UNAFF_n) %>%
        setNames(c("Chromosome", "WinStart", "WinEnd", "% of individuals with ROH", "N (SNPs)")) %>% 
        gt() %>% 
        tab_header(
                title = "Top 0.5% ROH islands",
                subtitle = "ROH density measured in 500Kb running windows"
        ) #%>% 
        #gtsave(filename = "figs/tables/roh_islands.png")


quants <- quantile(running_roh$UNAFF_mean, na.rm = TRUE, probs= c(0.01, 0.99))

# check that genome assembly is ok where deserts and islands are ---------------
# load linkage map with interpolated SNP positions
lmap <- read_delim("data/Oar3.1_Interpolated.txt", "\t") %>% 
        setNames(c("chr", "snp_name", "bp", "cM")) %>% 
        mutate(mb_pos = bp/1000000,
               kb_pos = bp/1000)

win_area <- function(win_mid, CHR, ...) {
        diffs <- lmap %>% 
                filter(chr == CHR) %>% 
                mutate(diff = abs(win_mid - kb_pos)) %>% 
                arrange(diff) %>% 
                top_n(-500)
}

# deserts
all_des <- pmap(roh_deserts, win_area) %>% 
                map(as_tibble) %>% 
                bind_rows(.id = "desert_num") %>% 
                mutate(desert = paste0(desert_num, " | Chr. ", chr),
                       desert = fct_inorder(desert))
 
all_des$win_start <- rep(roh_deserts$win_start, each = 500)
all_des$win_end <- rep(roh_deserts$win_end, each = 500)

# ROH deserts and islands tend to be very good at coinciding with genome
# assembly errors. Let's check that this is not the case by 
# by plotting genetic vs. physical SNP positions around deserts and islands.

# deserts
p_des_rec <- ggplot(all_des, aes(mb_pos, cM)) +
        geom_point(shape = 21, fill = "eceff4", stroke = 0.05, size = 2) +
        facet_wrap(~desert, scales = "free") +
        scale_x_continuous( breaks = scales::pretty_breaks(3)) +
        scale_y_continuous( breaks = scales::pretty_breaks(3)) +
        geom_vline(aes(xintercept = win_start/1000)) +
        geom_vline(aes(xintercept = win_end/1000)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
        xlab("Mb") +
        ggtitle("Genetic vs. physical SNP positions in regions with ROH deserts")
p_des_rec 
# ggsave("figs/pot_sup_ROH_des_rec.jpg", width = 8, height = 7)

# islands
all_isl <- pmap(roh_islands, win_area) %>% 
        map(as_tibble) %>% 
        bind_rows(.id = "island_num") %>% 
        mutate(island = paste0(island_num, " | Chr. ", chr),
               island = fct_inorder(island))

all_isl$win_start <- rep(roh_islands$win_start, each = 500)
all_isl$win_end <- rep(roh_islands$win_end, each = 500)

p_isl <- ggplot(all_isl, aes(mb_pos, cM)) +
        geom_point(shape = 21, fill = "eceff4", stroke = 0.05, size = 2) +
        facet_wrap(~island, scales = "free") +
        scale_x_continuous( breaks = scales::pretty_breaks(3)) +
        scale_y_continuous( breaks = scales::pretty_breaks(3)) +
        geom_vline(aes(xintercept = win_start/1000)) +
        geom_vline(aes(xintercept = win_end/1000)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
        xlab("Mb") +
        ggtitle("Genetic vs. physical SNP positions in regions with ROH islands")
p_isl
# ggsave("figs/pot_sup_ROH_isl_rec.jpg", width = 8, height = 7)


# ROH and recombination rate variation analysis --------------------------------
# roh prevalence across the genome
snp_roh <- fread("output/ROH/roh.hom.summary")
# linkage map from Johnston et al. (2016)
lmap <- read_delim("data/7_20200504_Full_Linkage_Map.txt", "\t") %>% 
        rename(SNP = SNP.Name)
# lmap <- read_delim("data/linkage_map_2020.txt", "\t") %>% 
#         rename(SNP = SNP.Name)
snp_df <- snp_roh %>% 
        inner_join(lmap, by = "SNP") %>% 
        mutate(KB = BP/1000)

# recombination fraction vs cM position
ggplot(snp_df, aes(r, cMdiff)) +
        geom_point() +
        geom_smooth(method = "lm") +
        facet_wrap(~CHR)

# mean recombination fraction in 500 Kb non-overlapping running windows
running_snp <- winScan(x = snp_df,
                       groups = "CHR",
                       position = "KB",
                       values = c("UNAFF", "r"),
                       win_size = 500,
                       win_step = 500,
                       funs = c("median", "sum"),
                       cores = 8)

# supplementary plot
running_snp %>% 
        mutate(prop_ROH = (UNAFF_median/5952) * 100) %>% 
        ggplot(aes(r_sum, prop_ROH)) +
        geom_point(shape = 21, fill = "#eceff4", stroke = 0.1, alpha = 1) +
        ylab("% of sheep with ROH") +
        xlab("Recombination rate") +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        geom_smooth(method = "lm", se = FALSE, color = "#5e81ac", size = 1) +
        scale_x_continuous( breaks = scales::pretty_breaks(2)) +
        scale_y_continuous( breaks = scales::pretty_breaks(4)) +
        theme(axis.title = element_text(size = 16)) +
        facet_wrap(~CHR, scales = "free") -> p_rec
p_rec
#ggsave("figs/Sup_ROH_Rec.jpg", width = 8, height = 7)

# model df
running_snp <- running_snp %>% 
                left_join(chr_sizes) %>% 
                mutate(size_MB = size_KB/1000) %>% 
                mutate(r_sum_std = scale(r_sum),
                       size_MB_std = scale(size_MB),
                       ROH_prop = UNAFF_median/5952)

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


#all in one plot
p_roh_rec <- running_snp %>% 
      #  filter(UNAFF_n > 35) %>% 
        mutate(prop_ROH = (UNAFF_median/5952) * 100) %>% 
        ggplot(aes(r_sum, prop_ROH)) +
        geom_point(shape = 21, fill = "#eceff4", stroke = 0.1, alpha = 0.5, size = 2.5) +
        ylab("% of sheep with ROH") +
        xlab("Recombination rate") +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        geom_smooth(method = "lm", se = FALSE, color = "#2e3440", size = 1) 
ggsave("figs/roh_rec.jpg", plot = p_roh_rec, width = 5, height = 3.7)

#geom_boxplot(outlier.shape = NA)


# recombination rate, heterozygosity and roh -----------------------------------
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
        mutate(KB = position / 1000) %>% 
        as_tibble(rownames = "snp_name") %>% 
        left_join(snps_stats, by = "snp_name") %>% 
        rename(het = P.AB)


running_snp_het <- winScan(x = snps_stats_full,
                       groups = "chromosome",
                       position = "KB",
                       values = c("het"),
                       win_size = 500,
                       win_step = 500,
                       funs = c("median", "sum"),
                       cores = 8)

snp_full <- running_snp %>% 
        left_join(running_snp_het)

ggplot(snp_full, aes(r_sum, het_median)) +
        geom_point() +
        geom_smooth(method = "lm") +
        facet_wrap(~CHR, scales = "free_x")

p1 <- ggplot(snp_full %>% filter(CHR %in% 1:3), aes(win_mid, het_median)) +
        geom_smooth() +
        facet_wrap(~CHR, nrow = 1, scales = "free_x")
p2 <- ggplot(snp_full %>% filter(CHR %in% 1:3), aes(win_mid, UNAFF_median)) +
        geom_smooth()+
        facet_wrap(~CHR, nrow = 1, scales = "free_x")

p1/p2

snp_sub <- snp_full %>% filter(CHR == 1)
cor(snp_sub$het_median, snp_sub$UNAFF_median, use = "complete.obs")

p1 <- snp_full %>% 
        mutate(win_num = 1:nrow(.)) %>% 
        ggplot( aes(win_num, het_median)) +
       # geom_smooth(n = 4)+
        stat_smooth(n = 10)

p2 <- snp_full %>% 
        mutate(win_num = 1:nrow(.)) %>% 
        ggplot( aes(win_num, UNAFF_median)) +
        stat_smooth(n = 10)

p1/p2

snps_stats_full %>% 
        filter(chromosome == 1) %>% 
        mutate(snp_num = 1:nrow(.)) %>% 
        ggplot(aes(snp_num, het)) + 
                geom_line(size = 0.01)

roh_d <- roh_deserts %>% 
        rename(chromosome = CHR) %>% 
        left_join(running_snp_het)
        
roh_i <- roh_islands %>% 
        rename(chromosome = CHR) %>% 
        left_join(running_snp_het)

roh_extremes <- bind_rows(roh_d, roh_i, .id = "extremes")

ggplot(roh_extremes, aes(extremes, het_median)) +
        geom_boxplot() 


running_snp %>% 
        rename(chromosome = CHR) %>% 
        left_join(roh_extremes, by = c("chromosome", "win_start")) %>% 
        mutate(extreme = case_when(
                (extremes == 1) ~ "island",
                (extremes == 2) ~ "desert"
        )) %>% 
        mutate(prop_ROH = (UNAFF_median/5952) * 100) %>% 
        ggplot(aes(r_sum, prop_ROH, fill = extreme)) +
        geom_point(shape = 21,  stroke = 0.1, alpha = 0.5, size = 2.5) +
        ylab("% of sheep with ROH") +
        xlab("Recombination rate") +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        geom_smooth(method = "lm", se = FALSE, color = "#2e3440", size = 1) 
      
