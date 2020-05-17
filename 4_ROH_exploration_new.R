# Exploring ROH through plots

library(data.table)
library(RColorBrewer)
library(wesanderson)
library(patchwork)
library(tidyverse)
source("theme_simple.R")
library(ggridges)
library(viridis)
library(GGally)
library("naniar")
options(scipen=999)
library(ggchicklet)
library(windowscanr)
library(cowplot)
library(gt)
library(grid)
library(ggplotify)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Full data           #   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Chr lengths
chr_data <- read_delim("../sheep/data/sheep_genome/chromosome_info_ram.txt", delim = "\t") %>% 
        rename(size_BP = Length,
               CHR = Part) %>% 
        mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
        .[2:27, ] %>% 
        summarise(sum_KB = sum(size_KB)) %>% 
        as.numeric()

#~~ ROH for survival data subset
file_path <- "output/ROH/roh_ram.hom"
roh_lengths <- fread(file_path)

# max ROH length
roh_lengths[which.max(roh_lengths$KB), ]

# ROH overview -----------------------------------------------------------------
# descriptive ROH statistics
num_roh_per_ind <- roh_lengths %>% group_by(IID) %>% tally() 
summary(num_roh_per_ind$n)
sd(num_roh_per_ind$n)

# GENOME SIZE NEEDS fixing (only autosomes)
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

plot(froh$KBSUM, num_roh_per_ind$n)

# ROH length and abundance in the least and most inbred individuals 
num_roh_per_ind %>% 
        left_join(froh) %>% 
        top_frac(-0.01, FROH) %>% 
        #top_frac(0.005, desc(FROH)) %>% 
        # arrange(desc(FROH)) %>% 
        summarise(mean(n), mean(KBAVG))

# Supplementary: FROH / ROH across individuals ---------------------------------
p_froh <- ggplot(froh, aes(FROH)) +
        geom_histogram(bins = 100, fill = "#E5E9F0", color = "black", size = 0.1) +
        ylab("individuals") +
        xlab(expression(inbreeding~coefficient~F[ROH])) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12, 
                     base_family = "Lato") 
p_froh

p_roh <- ggplot(num_roh_per_ind, aes(n)) +
        geom_histogram(binwidth = 1,  fill = "#E5E9F0", color = "black",  size = 0.1) +
        ylab("individuals") +
        xlab("ROH per genome") +
        scale_y_continuous(expand = c(0, 0)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12, 
                     base_family = "Lato") 
p_roh

p_roh_dist <- p_froh + p_roh + plot_annotation(tag_levels = 'A')
p_roh_dist
ggsave("figs/Sup_ROH_dist.jpg", p_roh_dist, width = 6, height = 2.2)

# Supplementary: plot HD vs. imputed individuals -------------------------------
# HD inds 
hd_inds <- read_delim("../sheep/data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram.fam", " ", col_names = FALSE)[[2]]
froh_plot <- froh %>% mutate(hd = ifelse(IID %in% hd_inds, "HD", "imputed")) %>% 
        mutate(hd = as.factor(hd)) %>% mutate(hd = relevel(hd, "imputed"))

# medians
froh_plot %>% group_by(hd) %>% summarise(median(FROH))

#
froh_imp_vs_hd <- ggplot(froh_plot, aes(x=FROH, fill = relevel(hd, "imputed"))) +
        geom_histogram(bins = 100,  color = "black", size = 0.1, position = "identity",
                       alpha = 0.8) +
        ylab("individuals") +
        xlab(expression(inbreeding~coefficient~F[ROH])) +
        scale_fill_manual(name = "Genotypes", values = c("#E5E9F0", "black"),
                          labels = c("partially imputed", "high-density")) +
        scale_y_log10(expand = c(0, 0)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12,  # "#E5E9F0"
                     base_family = "Lato") +
        theme(legend.position="right")
froh_imp_vs_hd
ggsave("figs/Sup_imp_vs_hd.jpg", froh_imp_vs_hd, width = 6, height = 2.7)


# ROH length classes -----------------------------------------------------------
# Specifiy length distributon:
# a separation of m meisies (2m generations) results in segment lengths that are
# exponentially distributied with mean 100/m cM

# 1.451221 cM/Mb for Males
# 1.037693 cM/Mb for Females
# 1.279305 cM/Mb sex-averaged
# 
# 1cM is 0.6890747 in males
# 1cM is 0.9636761 in females
# 1cM is 0.7816743 sex-averaged

length_dist <- data.frame(g = c(1, 2,2^2, 2^3, 2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13)) %>%
        mutate(ROH_length_cM = 100 / (2*g)) %>% 
        mutate(ROH_length_Mb = ROH_length_cM * 0.7816743)

prop_IBD_df <- roh_lengths %>%
        mutate(length_Mb = KB/1000) %>%
        mutate(class = case_when(length_Mb >= 39.083715000 ~ 1,
                                 length_Mb < 39.083715000 & length_Mb >= 19.541857500 ~ 2,
                                 length_Mb < 19.541857500 & length_Mb >= 9.770928750 ~ 4,
                                 #length_Mb < 9.770928750 & length_Mb >= 6.513952500 ~ 6,
                                 length_Mb < 9.770928750& length_Mb >= 4.885464375 ~ 8,
                                 # length_Mb < 4.885464375 & length_Mb >= 3.908371500 ~ 10,
                                 length_Mb < 4.885464375 & length_Mb >= 2.442732188 ~ 16,
                                 length_Mb < 2.442732188 & length_Mb >= 1.221366094 ~ 32,
                                 length_Mb < 1.221366094 & length_Mb >= 0.6 ~ 64)) %>% 
        #length_Mb < 0.610683047 & length_Mb >= 0.30 ~ 128)) %>% # 0.610683047
        mutate(length_class = case_when(
                class == 1 ~ ">39 (1G)",
                class == 2 ~ ">19.5-39 (>1-2G)",
                class == 4 ~ ">9.8-19.5 (>2-4G)",
                # class == 6 ~ "6.5-9.7 (6G)",
                class == 8 ~ ">4.9-6.5 (>4-8G)",
                # class == 10 ~ "3.9-4.9 (10G",
                class == 16 ~ ">2.4-4.9 (>8-16G)",
                class == 32 ~ ">1.2-2.4 (>16-32G)",
                class == 64 ~ ">0.6-1.2 (>32-64G)"
                # class == 128 ~ "0.6-0.3 (128G)"
        )) %>% 
        mutate(length_class = fct_reorder(length_class, class)) %>% 
        mutate(IID = as.character(IID)) %>% 
        group_by(IID, class, length_class) %>%
        dplyr::summarise(prop_IBD = sum(length_Mb / (autosomal_genome_size/1000))) #%>% 
# add IBD of non-ROH snps if wanted
#  bind_rows(homs) 

prop_IBD_df_with_0 <- prop_IBD_df %>% 
        # add missing length classes as 0
        ungroup() %>% 
        tidyr::complete(length_class, nesting(IID)) %>% 
        mutate(class = ifelse(is.na(class), length_class, class)) %>% 
        mutate(prop_IBD = ifelse(is.na(prop_IBD), 0, prop_IBD))

# ROH class distribution across individuals
col_pal <- rev(c(brewer.pal(7, "YlGnBu"))) ##A42820
col_pal <- viridis(7)
prop_IBD_df
set.seed(17)
library(gghalves)

# ROH classes summary
prop_IBD_df_with_0 %>% group_by(length_class) %>% summarise(mean(prop_IBD), sd(prop_IBD))
prop_IBD_df_with_0 %>% group_by(length_class) %>% summarise(prop = sum(prop_IBD > 0) / n() )

# plot
p_roh_classes <- prop_IBD_df_with_0 %>% 
        mutate(prop_IBD = prop_IBD * 100) %>% 
        ggplot(aes(length_class, prop_IBD, fill = length_class)) +
        geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =2,
                        transformation_params = list(height = 0, width = 1.3, seed = 1)) +
        geom_half_boxplot(side = "r", outlier.color = NA,
                          width = 0.6, lwd = 0.3, color = "black",
                          alpha = 0.8) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13) +
        ylab("% genome") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
        theme(legend.position = "none",
              #axis.ticks.x = element_blank(),
              axis.title=element_text(size = rel(1.1)), 
              axis.text = element_text(color = "black")) + 
        xlab("ROH classes in Mb") 

ggsave("figs/fig1a_roh_classes_boxplots.jpg", p_roh_classes, width = 7, height = 3)


#ggsave("figs/roh_classes_subset_inds.jpg", p_roh2, width = 15, height = 5)

#~~ ROH for some indiividuals --------------------------------------------------
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

library(ghibli)
# col <- wes_palette("Darjeeling2")[c(1,2)]
# col <- viridis(2, begin = 0.1, end = 0.9)
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

#png(file="figs/ROH_map.png", units = "in", res = 300, height=4, width=7)

col <- c("#3690c0", "#a6bddb")
col <- c("#4c566a", "#d8dee9")
#col <- c("#eceff4", "#d8dee9")
#col <- viridis(2)
#col <- viridis(6, option = "D")[c(2,6)]
chr_names <- as.character(1:26)
names(chr_names) <- as.character(1:26)
chr_names[c(11, 13, 15, 17, 19, 21, 23, 25)] <- ""

df %>% 
        filter(MB > 5) %>% 
        filter(CHR %in% 1:26) %>% 
        #mutate(CHR = factor(CHR, levels = as.character(1:26))) %>% 
        ggplot() +
        geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
                  alpha=0.5, fill = "#eceff4") + # "#f7f7f7" "#eceff4"
        geom_hline(data = yax, aes(yintercept = yax), color = "#d8dee9", size = 0.4) +
        geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.9, 
                      fill = as.factor(CHR)),  col = "black", size = 0.1, alpha = 1) + 
        scale_fill_manual(values = rep(col, 18)) + 
        scale_color_manual(values = rep(col, 18)) +
        # scale_x_discrete(labels = as.character(c(1:20, 22, 24, 26))) + 
        scale_y_reverse(expand = c(0, 0)) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13) +
        facet_grid(~CHR,scales = 'free_x', space = 'free_x', switch = 'x',
                   labeller = as_labeller(chr_names)) +
        #ggtitle("ROH in the 7 most and least inbred individuals") +
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
                #strip.text.x = element_text(size = 12),
                #plot.margin = unit(c(1, 1, 1, 1), "cm"), 
                axis.text.y = element_text(colour = "white"),
                axis.line.y = element_blank()) +
               # axis.title=element_text(size=17))+
        coord_cartesian(clip = 'off') +
        xlab("Chromosome") +
        ylab("Individuals") -> ROH_per_ind
#ggtitle("ROH > 1Mb in the 10 most and 10 least inbred individuals") 
ROH_per_ind
#ggsave("figs/fig1b_roh_per_ind_5Mb.jpg", ROH_per_ind, width = 7, height = 3)
pg <- ggplotGrob(ROH_per_ind)

for(i in which(grepl("strip-b", pg$layout$name))){
        pg$grobs[[i]]$layout$clip <- "off"
}

ROH_per_ind_grob <- as.ggplot(pg)
#ggsave("figs/fig1b_roh_per_ind_5Mb.jpg", ROH_per_ind_grob, width = 6, height = 3.5)

#p1 / p_roh_classes + plot_layout(heights = c(1, 0.8))

#~~~ ROH density ---------------------------------------------------------------
# do not use pruned data here as this will bias ROH densities
#hom_sum <- fread("output/ROH/roh_nofilt_ram_pruned.hom.summary")
hom_sum <- fread("output/ROH/roh_ram.hom.summary")
#hom_sum <- fread("output/ROH/ROH_surv_subset/roh_nofilt_ram_pruned.hom.summary")
# snps that can have roh
# snps_pos <- read_delim("output/snps_that_can_have_roh_190k", " ")
# hom_sum <- hom_sum %>% filter(snps_pos$snps_ok == "yes")

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

cor(running_roh$UNAFF_mean, running_roh$UNAFF_n, use = "complete.obs")
#filter(CHR == 1) %>% 
running_roh %>% 
        mutate(UNAFF_mean = UNAFF_mean/5952) %>% 
        filter(UNAFF_n > 0) -> running_roh_p

# check distribution of snps
#ggplot(running_roh, aes(win_start, UNAFF_n)) + 
#  geom_point() + geom_smooth(span = 0.1) +
#  facet_wrap(~CHR, scales = "free_x")

# ROH across the genome plot
library(scales)
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
#axis.ticks.x = element_line(size = 1)) 

p2
layout <- c(
        patchwork::area(t = 1, l = 1, b = 5, r = 5),
        patchwork::area(t = 4, l = 3, b = 4, r = 4)
)
p_final <- p1 + p2 + plot_layout(design = layout)
p_final
ggsave("figs/roh_genome_398K.jpg", p1, width = 8, height = 6)
ggsave("figs/roh_genome_legend_398K.jpg", p2, width = 5, height = 3)



# try simple combined plot
p_roh_comb_simple <- plot_grid(ROH_per_ind_grob, p1, nrow = 2, 
                               rel_heights = c(0.658, 1), label_size = 15, 
                               labels = c("A", "B"), align = "v")
p_roh_comb_simple
ggsave("figs/roh_patterns_simple3.jpg", p_roh_comb_simple, width = 7, height = 6.5)





#ggdraw(p_roh_comb_simple) + draw_plot(p2, x = 0.1, y = -0.2, scale = 0.47)


p_roh_comb <- plot_grid(p_roh_classes, ROH_per_ind, nrow = 2, rel_heights = c(2.2,3, 5.2))
p_roh_comb
#ggsave("figs/roh_comb.jpg", p_roh_comb, width = 8, height = 6)


fitness_data

# ROH ISLANDS AND DESERTS

#~~~ ROH density
hom_sum <- fread("output/ROH/roh_ram.hom.summary") # ROH_surv_subset/
head(hom_sum)

hom_sum <- hom_sum %>%
        mutate(MB = BP / 1000000,
               KB = BP / 1000,
               index = 1:nrow(.))# %>% 
# filter(UNAFF > 0) 

# non of the SNPs where UNAFF = 0 can actually have an roh (see script snps_that_can_have_roh.R)
# hom_sum_filt
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

#running_roh <- winScan(x = hom_sum,
#                       groups = "CHR",
#                       #position = "KB",
#                       values = "UNAFF",
#                       win_size = 5,
#                       win_step = 3,
#                       funs = c("mean", "var"))

# check dist
hist(running_roh$UNAFF_mean, xlim = c(0,6000), breaks = 1000)

library(gt)
# roh desert
roh_deserts <- running_roh %>% 
        filter(UNAFF_n > 10) %>% 
        mutate(prop_roh = UNAFF_mean/5952) %>%  # 7691 5952
        arrange(prop_roh) %>% 
        # top 0.5% of windows
        .[1:26, ]
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
        ) %>% 
        gtsave(filename = "figs/tables/roh_desert.png")

roh_islands <- running_roh %>% 
        filter(UNAFF_n > 10) %>% 
        mutate(prop_roh = UNAFF_mean/5952) %>%  # 7691 5952
        arrange(desc(prop_roh)) %>% 
        .[1:26, ]

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
        ) %>% 
        gtsave(filename = "figs/tables/roh_islands.png")


quants <- quantile(running_roh$UNAFF_mean, na.rm = TRUE, probs= c(0.01, 0.99))


#










#




# manhattan plot ROH distribution
hom_sum <- fread("output/ROH/roh_nofilt_ram.hom.summary") %>% 
        rename(snp.name = SNP, roh_count = UNAFF, chromosome = CHR,
               position = BP) %>% 
        dplyr::select(snp.name, roh_count, position, chromosome) 


chr_info <- read_delim("../sheep/data/sheep_genome/chromosome_info_ram.txt", "\t") %>% 
        .[-1, ] %>% 
        rename(chromosome = Part) %>% 
        mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
        mutate(chromosome = as.integer(chromosome)) %>% 
        filter(!is.na(chromosome))

chr_labels <- c(c(1:18),"","20","",  "22","", "24","", "26")
cols <- c("#336B87", "#2A3132")

# get cumsums
chr_info %<>% 
        mutate(tot=cumsum(Length)-Length) %>% 
        dplyr::select(chromosome, tot)

gwas_p <- hom_sum  %>% 
        left_join(chr_info) %>% 
        # Add a cumulative position of each SNP
        arrange(chromosome, position) %>%
        mutate(pos_cum = position + tot) 

axisdf <- gwas_p %>% group_by(chromosome) %>% 
        summarize(center = (max(pos_cum) + min(pos_cum)) / 2 )

cols <- c("#4393c3", "#2A3132")
gwas_p <- gwas_p %>% mutate(roh_prop = roh_count / 7691)
gwas_p %>% 
        sample_frac(0.1) %>% 
        ggplot(aes(x=pos_cum, y=roh_prop)) +
        # Show all points
        geom_point(aes(fill = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
                   size = 2, shape = 21, stroke = 0.05) + # alpha = 1, 
        scale_x_continuous(labels = chr_labels, breaks= axisdf$center ) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
        xlab("Chromosome") + 
        ylab("Proportion ROH") + ## y label from qqman::qq
        scale_fill_manual(values = alpha(cols, 0.2)) +##values = cols
        theme_simple() +
        theme(axis.text.x = element_text(size = 10),
              axis.ticks = element_line(size = 0.1)) +
        geom_hline(yintercept = 0.9, size = 0.1, linetype = "dashed") +
        geom_hline(yintercept = 0.1, size = 0.1, linetype = "dashed") +
        geom_hline(yintercept = 0.32, size = 0.1, linetype = "dashed") +
        #geom_hline(yintercept =1, size = 0.1, linetype = "dashed") +
        guides(fill=FALSE) 

pgwas
ggsave("figs/roh_prop_genome_manhattan.jpg", width = 10, height = 3)


# per chr

# Show all points
hom_sum %>% 
        filter(chromosome %in% c(1,2,3,4,5,6)) %>% 
        ggplot(aes(position, roh_count)) +
        geom_point(aes(fill = chromosome %%2 == 0),#shape = roh_prevalence  #fill = chromosome %%2 == 0
                   size = 2, shape = 21, stroke = 0.05) + # alpha = 1, 
        #scale_x_continuous(labels = chr_labels) +
        #scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
        xlab("Chromosome") + 
        ylab("Proportion ROH") + ## y label from qqman::qq
        scale_fill_manual(values = alpha(cols, 0.2)) +##values = cols
        theme_simple(base_size = 5) +
        theme(axis.text.x = element_text(size = 10),
              axis.ticks = element_line(size = 0.1)) +
        geom_hline(yintercept = 0.9, size = 0.1, linetype = "dashed") +
        geom_hline(yintercept = 0.1, size = 0.1, linetype = "dashed") +
        geom_hline(yintercept = 0.32, size = 0.1, linetype = "dashed") +
        #geom_hline(yintercept =1, size = 0.1, linetype = "dashed") +
        guides(fill=FALSE) +
        facet_wrap(~chromosome, nrow = 26)

# correlation of ROH across classes plot
ROH_classes_pairs <- prop_IBD_df_with_0 %>% 
        split(prop_IBD_df_with_0$length_class) %>% 
        reduce(left_join, by = "IID") %>% 
        dplyr::select(contains("prop_IBD")) %>% 
        rename_all(~ levels(prop_IBD_df_with_0$length_class)) %>% 
        # sample_frac(0.01) %>% 
        # mutate_all(funs(ifelse(. == 0, NA, .)))
        replace_with_na_all(condition = ~.x == 0) %>% 
        ggscatmat(alpha = 0.3) +
        geom_smooth(method = "lm") +
        theme_simple() +
        theme(axis.text = element_blank()) +
        xlab("FROH") +
        ylab("FROH")

ROH_classes_pairs 
#ggsave("figs/SUP1_roh_classes_pairs.jpg", ROH_classes_pairs, width = 8, height = 7)

# ROH per individual plot
prop_IBD_df %>% 
        ungroup() %>% 
        mutate(IID = as.character(IID)) %>% 
        filter(IID %chin% sample(as.character(unique(prop_IBD_df$IID)), 
                                 100, replace = FALSE)) -> plot_ibd_df

col_pal <- rev(c(brewer.pal(7, "YlGnBu"),  "firebrick1")) ##A42820
#col_pal <- rev(c(brewer.pal(6, "YlGnBu"), "goldenrod", "firebrick1")) 
set.seed(153)
prop_IBD_df %>% 
        ungroup() %>% 
        # filter(class %in% c(1,2,4)) %>% 
        filter(IID %chin% sample(as.character(unique(prop_IBD_df$IID)), 500, replace = FALSE)) %>% 
        mutate(prop_IBD = prop_IBD * 100) %>% 
        ggplot(aes(x=IID, y=prop_IBD, fill=length_class)) +
        geom_bar(stat="identity") +
        # scale_fill_brewer(palette = "YlGnBu", name = "ROH class (Mb)", direction = -1) +
        scale_fill_manual(values = col_pal, name = "ROH length class \nin Mb (Generations)") +
        #facet_grid(.~ pop, space = 'free_x', scales = 'free_x', switch = 'x') +
        theme_simple(axis_lines = TRUE) +
        theme(axis.text.x = element_blank()) +
        #ylab("Proportion of genome in ROH") + 
        ylab("% genome") +
        xlab("Individuals") + 
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) +
        theme(legend.position = "right",
              axis.ticks.x = element_blank()) -> p_roh2
p_roh2
