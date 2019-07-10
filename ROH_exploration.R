# Exploring ROH through plots

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(forcats)
library(readxl)
library(tidyverse)
library(zoo)
source("theme_clean.R")
options(scipen=999)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Full data           #   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Chr lengths
chr_data <- read_delim("../sheep/data/chromosome_data.txt", delim = "\t") %>% 
                        mutate(size_KB = `Size (Mb)` * 1000,
                               size_BP = size_KB * 1000) %>% 
                        rename(CHR = Name) %>% 
                        select(CHR, size_KB, size_BP)
                

#2570000
#2439314
# sheep gen    me size 2869490 KB 
#~~ Length distribution
file_path <- "output/ROH/roh_nofilt.hom"
file <- "roh_nofilt"
roh_lengths <- fread(file_path)

froh <- roh_lengths %>%
        dplyr::group_by(FID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/2869490)

# FROH across individuals
p_dist <- ggplot(froh, aes(FROH)) +
    geom_histogram(bins = 1000) +
    theme_clean() +
    ylab("individual count") +
    geom_vline(xintercept = min(froh$FROH), alpha = 0.5) +
    annotate("text",  min(froh$FROH) + 0.04, 80, label = paste0("min = ", round(min(froh$FROH), 2))) +
    geom_vline(xintercept = median(froh$FROH), alpha = 0.8) +
    annotate("text",  median(froh$FROH) + 0.05, 80, label = paste0("median = ", round(median(froh$FROH), 2))) +
    geom_vline(xintercept = max(froh$FROH), alpha = 0.5) +
    annotate("text",  max(froh$FROH) - 0.04, 80, label = paste0("max = ", round(max(froh$FROH), 2))) 
p_dist

ggsave("figs/FROH_dist.jpg", p_dist, width = 5, height = 3)


# Specifiy length distributon:
# a separation of m meisies (2m generations) results in segment lengths that are
# exponentially distributied with mean 100/m cM

length_dist <- data.frame(g = c(1, 2,2^2,2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13)) %>%
        mutate(ROH_length_Mb = 100 / (2*g))

prop_IBD_df <- roh_lengths %>%
        mutate(length_Mb = KB/1000) %>%
        mutate(class = case_when(length_Mb >= 50.00000000 ~ 1,
                                 length_Mb < 50.000000000 & length_Mb >= 25.000000 ~ 2,
                                 length_Mb < 25.000000000 & length_Mb >= 12.500000000 ~ 4,
                                 length_Mb < 12.500000000 & length_Mb >= 6.250000000 ~ 8,
                                 length_Mb < 6.250000000 & length_Mb >= 3.125000000 ~ 16,
                                 length_Mb < 3.125000000 & length_Mb >= 1.562500000 ~ 32,
                                 length_Mb < 1.562500000 & length_Mb >= 0.781250000 ~ 64)) %>%
        mutate(length_class = case_when(
          class == 1 ~ "> 50 (1G)",
          class == 2 ~ "25-50 (2G)",
          class == 4 ~ "12.5-25 (4G)",
          class == 8 ~ "6.3-12.5 (8G)",
          class == 16 ~ "3.13-6.3 (16G)",
          class == 32 ~ "1.6-3.13 (32G)",
          class == 64 ~ "0.78-1.6 (64G)"
        )) %>% 
        mutate(length_class = fct_reorder(length_class, class)) %>% 
        mutate(FID = as.character(FID)) %>% 
        group_by(FID, class, length_class) %>%
        dplyr::summarise(prop_IBD = sum(length_Mb / 2869))

# pop wide
col <- wes_palette("Moonrise2")
col <- col[c(1,2,3)]

p_roh <- ggplot(prop_IBD_df, aes(length_class, prop_IBD)) +
        geom_jitter(alpha = 0.3, width = 0.2) +
        geom_boxplot(outlier.shape = NA, color = "lightgrey", alpha = 0.7) +
        theme_clean() +
        #theme(axis.ticks.x = element_line(colour = "#cccccc", size = 0.3)) +
        xlab("ROH class (Mb)") + ylab("Proportion of genome in ROH")
p_roh
ggsave("figs/roh_classes_boxplots.jpg", p_roh, width = 8, height = 5)

# ind basis
prop_IBD_df %>% 
  ungroup() %>% 
  mutate(FID = as.character(FID)) %>% 
  filter(FID %chin% sample(as.character(unique(prop_IBD_df$FID)), 
                           100, replace = FALSE)) -> plot_ibd_df

col_pal <- rev(c(brewer.pal(6, "YlGnBu"), "#A42820"))

prop_IBD_df %>% 
  ungroup() %>% 
  # filter(class %in% c(1,2,4)) %>% 
  filter(FID %chin% sample(as.character(unique(prop_IBD_df$FID)), 400, replace = FALSE)) %>% 
  ggplot(aes(x=FID, y=prop_IBD, fill=length_class)) +
          geom_bar(stat="identity") +
         # scale_fill_brewer(palette = "YlGnBu", name = "ROH class (Mb)", direction = -1) +
          scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
          #facet_grid(.~ pop, space = 'free_x', scales = 'free_x', switch = 'x') +
          theme_clean() +
          theme(axis.text.x = element_blank()) +
          ylab("Proportion of genome in ROH") + 
          xlab("Individuals") + 
          scale_x_discrete(expand = c(0, 0)) + 
          scale_y_continuous(expand = c(0, 0)) -> p_roh2
p_roh1
ggsave("figs/roh_classes_subset_inds.jpg", p_roh2, width = 20, height = 5)




#~~~ ROH density

hom_sum <- fread("output/ROH/roh_nofilt.hom.summary")

hom_sum <- hom_sum %>%
  mutate(MB = BP / 100000,
         index = 1:nrow(.))


# Here's how you can do this with the dplyr functions group_by and do:
window_width <- 200
jumps <- 200
running_roh <- hom_sum %>% 
  group_by(CHR) %>% 
  do(
    data.frame(
      window.start = rollapply(.$BP, width= window_width, by=jumps, FUN=min, align="left"),
      window.end = rollapply(.$BP, width=window_width, by=jumps, FUN=max, align="left"),
      ninds = rollapply(.$UNAFF, width=window_width, by=jumps, FUN=mean, align="left")
    )
  )

running_roh$index = 1:nrow(x)
p_running_roh <- running_roh %>% 
  filter(CHR %in% 1:26) %>% 
  ggplot(aes(window.start, ninds)) +
  geom_line() +
  scale_y_continuous(breaks = c(1844, 5531), labels = c("25%", "75%")) +
  theme_clean() + 
  ggtitle("Average ROH across the genome of 7374 Soay sheep") +
  xlab("Window start (in Base Pairs) of 200 BP window") +
  ylab("Individuals with ROH") + 
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=1844, ymax=5531, alpha=0.1) +
  facet_grid(CHR~.) 

ggsave("figs/roh_across_genome.jpg", p_running_roh, width = 8, height = 8)


#~~ Runs of ROH
df <- roh_lengths %>%
        mutate(POS1 = POS1 / 1e+6,
               POS2 = POS2 / 1e+6,
               MB = KB / 1000)
num_ind <- 30
df <- df %>% filter(FID %in% sample(unique(df$FID), num_ind))

yax <- data.frame(FID = unique(df$FID)) %>%
        mutate(yax = seq(from = 2,
                         to = 2*length(unique(df$FID)),
                         by = 2)) 

df <- left_join(df, yax, by = "FID")

col <- wes_palette("Moonrise2")
col <- col[c(3,1,2)]

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
df %>% 
  filter(MB > 5) %>% 
  filter(CHR %in% 1:26) %>% 
  ggplot() +
        geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
                  fill = "grey90", alpha=0.3) +
        geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, 
                               ymax = yax + 0.5, col = as.factor(CHR))) +
        #
       # scale_fill_manual(values = col) +
       # scale_color_manual(values = col) +
        geom_hline(aes(yintercept = yax), color = "grey") +
        theme_clean() +
        facet_grid(~CHR, scales = 'free_x', space = 'free_x', switch = 'x') +
        theme(strip.placement = 'outside',
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.spacing = unit(0, "lines"),
              legend.position="none",
              plot.margin = unit(c(1, 1, 1, 1), "cm"), 
              axis.text.y = element_text(colour = "white")) +
        xlab("CHR") +
        ggtitle("ROH > 5Mb in 30 individuals") -> ROH_per_ind

ROH_per_ind
#dev.off()
ggsave("figs/roh_per_ind.jpg", ROH_per_ind, width = 12, height = 8)

