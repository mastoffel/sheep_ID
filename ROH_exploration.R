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
                        dplyr::select(CHR, size_KB, size_BP)
                

# add homozygosity not in roh
homs <- read_delim("output/ROH/roh_nofilt_hom_not_in_roh.txt", delim = " ") %>% 
            rename(prop_IBD = hom, 
                   FID = ID) %>% 
            mutate(FID = as.character(FID)) %>% 
            mutate(length_class = as.factor(rep("not_in_ROH", nrow(.))))
                

#2439314
# sheep gen    me size 2869490 KB 
#~~ Length distribution
file_path <- "output/ROH/roh_nofilt.hom"
file <- "roh_nofilt"
roh_lengths <- fread(file_path)

froh <- roh_lengths %>%
        dplyr::group_by(FID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/2619054)

# FROH across individuals
p_dist <- ggplot(froh, aes(FROH)) +
    geom_histogram(bins = 500) +
    theme_clean() +
    ylab("individual count") 
p_dist

ggsave("figs/FROH_dist.jpg", p_dist, width = 5, height = 3)


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
                                 length_Mb < 1.221366094 & length_Mb >= 0.610683047 ~ 64,
                                 length_Mb < 0.610683047 & length_Mb >= 0.30 ~ 128)) %>% # 0.610683047
        mutate(length_class = case_when(
          class == 1 ~ "> 39 (1G)",
          class == 2 ~ "19.5-39 (2G)",
          class == 4 ~ "9.8-19.5 (4G)",
         # class == 6 ~ "6.5-9.7 (6G)",
          class == 8 ~ "4.9-6.5 (8G)",
         # class == 10 ~ "3.9-4.9 (10G",
          class == 16 ~ "2.4-4.9 (16G)",
          class == 32 ~ "1.2-2.4 (32G)",
          class == 64 ~ "0.6-1.2 (64G)",
          class == 128 ~ "0.6-0.3 (128G)"
        )) %>% 
        mutate(length_class = fct_reorder(length_class, class)) %>% 
        mutate(FID = as.character(FID)) %>% 
        group_by(FID, class, length_class) %>%
        dplyr::summarise(prop_IBD = sum(length_Mb / 2619)) #%>% 
        # add IBD of non-ROH snps if wanted
      #  bind_rows(homs) 

prop_IBD_df_with_0 <- prop_IBD_df %>% 
        # add missing length classes as 0
        ungroup() %>% 
        tidyr::complete(length_class, nesting(FID)) %>% 
        mutate(class = ifelse(is.na(class), length_class, class)) %>% 
        mutate(prop_IBD = ifelse(is.na(prop_IBD), 0, prop_IBD))
 

# ROH class distributions as prop. 

prop_IBD_df
p_roh <- ggplot(prop_IBD_df_with_0, aes(length_class, prop_IBD)) +
        geom_jitter(alpha = 0.3, width = 0.2) +
        geom_boxplot(outlier.shape = NA, color = "lightgrey", alpha = 0.7) +
        theme_clean() +
        #theme(axis.ticks.x = element_line(colour = "#cccccc", size = 0.3)) +
        xlab("ROH class (Mb)") + ylab("Proportion of genome in ROH")
p_roh
ggsave("figs/roh_classes_boxplots.jpg", p_roh, width = 9, height = 5)

library(ggridges)
library(viridis)

IBD_across_classes <- prop_IBD_df %>% 
  #sample_frac(0.05) %>% 
  ggplot(aes(prop_IBD, length_class)) +
  #geom_density_ridges(scale = 2, jittered_points=TRUE, point_shape = "|",
  #                    position = position_points_jitter(height = 0),) +
  stat_density_ridges(quantile_lines = TRUE,bandwidth = 0.003) + # , bandwidth = 0.002
  scale_y_discrete(limits = rev(levels(prop_IBD_df$length_class))) +
  ylab("ROH length class in Mb (generations)") +
  xlab("Proportion of the genome") + 
  #scale_fill_viridis(name = "Temp. [F]", option = "C")  +
  theme_clean()
IBD_across_classes
ggsave("figs/roh_classes_ridgeplot.jpg", IBD_across_classes, width = 8, height = 8)


library(GGally)
library("naniar")

?ggpairs

ROH_classes_pairs <- prop_IBD_df_with_0 %>% 
  split(prop_IBD_df_with_0$length_class) %>% 
  reduce(left_join, by = "FID") %>% 
  dplyr::select(contains("prop_IBD")) %>% 
  rename_all(~ levels(prop_IBD_df_with_0$length_class)) %>% 
 # sample_frac(0.01) %>% 
 # mutate_all(funs(ifelse(. == 0, NA, .)))
  replace_with_na_all(condition = ~.x == 0) %>% 
  ggscatmat(alpha = 0.3) +
  geom_smooth(method = "lm") +
  theme_clean() +
  xlab("FROH") +
  ylab("FROH")

ggsave("figs/roh_classes_pairs.jpg", ROH_classes_pairs, width = 13, height = 12)


# ind basis
prop_IBD_df %>% 
  ungroup() %>% 
  mutate(FID = as.character(FID)) %>% 
  filter(FID %chin% sample(as.character(unique(prop_IBD_df$FID)), 
                           100, replace = FALSE)) -> plot_ibd_df

col_pal <- rev(c(brewer.pal(7, "YlGnBu"),  "firebrick1")) ##A42820
#col_pal <- rev(c(brewer.pal(6, "YlGnBu"), "goldenrod", "firebrick1")) 
prop_IBD_df %>% 
  ungroup() %>% 
  # filter(class %in% c(1,2,4)) %>% 
  filter(FID %chin% sample(as.character(unique(prop_IBD_df$FID)), 500, replace = FALSE)) %>% 
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
          scale_y_continuous(expand = c(0, 0)) +
          theme(legend.position = "bottom") -> p_roh2
p_roh2
ggsave("figs/roh_classes_subset_inds.jpg", p_roh2, width = 15, height = 5)


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

# running_roh$index = 1:nrow(x)
p_running_roh <- running_roh %>% 
  filter(CHR %in% 1:26) %>% 
  ggplot(aes(window.start, ninds)) +
  geom_line() +
  scale_y_continuous(breaks = c(369, 1844, 5531, 7005), labels = c("5%", "25%", "75%", "95%"),
                     limits = c(0,7374)) +
  theme_clean() + 
  ggtitle("Average ROH across the genome of 7374 Soay sheep") +
  xlab("Window start (in Base Pairs) of 200 BP window") +
  ylab("Individuals with ROH") + 
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=1844, ymax=5531, alpha=0.1) +
  geom_hline(yintercept = 369, size = 0.3) +
  geom_hline(yintercept = 7005, size = 0.3) +
  facet_grid(CHR~.) 
p_running_roh 

ggsave("figs/roh_across_genome.jpg", p_running_roh, width = 8, height = 12)



# variance per IBD class
library(mclust)
?mclust
roh_lengths %>% 
  group_by(FID) %>% 
  summarise(KB = mean(KB))
  
out <- Mclust(roh_lengths$KB)
summary(out)
plot(out)


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

