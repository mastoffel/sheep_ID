# Exploring ROH through plots

library(data.table)
library(RColorBrewer)
library(wesanderson)
library(readxl)
library(tidyverse)
library(zoo)
source("theme_simple.R")
library(ggridges)
library(viridis)
library(GGally)
library("naniar")
options(scipen=999)
library(ggchicklet)

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

#~~ Length distribution
file_path <- "output/ROH/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path)

# GENOME SIZE NEEDS fixing (only autosomes)
froh <- roh_lengths %>%
        dplyr::group_by(IID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/autosomal_genome_size) %>% 
        mutate(FROH_cent = FROH - mean(FROH))

# FROH across individuals
p_dist <- ggplot(froh, aes(FROH)) +
    geom_histogram(bins = 100, fill = "#E5E9F0", color = "black", size = 0.2) +
    ylab("individuals") +
    xlab(expression(inbreeding~coefficient~F[ROH])) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12, 
                base_family = "Lato")
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
                                 length_Mb < 1.221366094 & length_Mb >= 0.6 ~ 64)) %>% 
                                 #length_Mb < 0.610683047 & length_Mb >= 0.30 ~ 128)) %>% # 0.610683047
        mutate(length_class = case_when(
          class == 1 ~ "> 39 (1G)",
          class == 2 ~ "19.5-39 (2G)",
          class == 4 ~ "9.8-19.5 (4G)",
         # class == 6 ~ "6.5-9.7 (6G)",
          class == 8 ~ "4.9-6.5 (8G)",
         # class == 10 ~ "3.9-4.9 (10G",
          class == 16 ~ "2.4-4.9 (16G)",
          class == 32 ~ "1.2-2.4 (32G)",
          class == 64 ~ "0.6-1.2 (64G)"
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
 
# ROH class distributions as prop. 
prop_IBD_df
p_roh <- ggplot(prop_IBD_df_with_0, aes(length_class, prop_IBD, fill = length_class)) +
        geom_jitter(alpha = 0.6, width = 0.2, shape = 21, color = "black", stroke = 0.1) +
        geom_boxplot(outlier.shape = NA, color = "#393e46", alpha = 0.7, width = 0.6,
                     lwd = 0.4) +
        #theme_classic() + 
        theme_simple() +
        scale_fill_brewer(palette = "GnBu", direction = -1) +
        #scale_fill_brewer(palette = "YlGnBu", direction = -1) + 
        theme(legend.position = "none") + 
        #theme(axis.ticks.x = element_line(colour = "#cccccc", size = 0.3)) +
        xlab("ROH class in Mb (Generations)") + ylab("Proportion of genome")
p_roh

ggsave("figs/roh_classes_boxplots.jpg", p_roh, width = 9, height = 3.5)



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
  theme_simple()
IBD_across_classes
ggsave("figs/roh_classes_ridgeplot.jpg", IBD_across_classes, width = 8, height = 8)


?ggpairs

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

ggsave("figs/roh_classes_pairs.jpg", ROH_classes_pairs, width = 8, height = 7)


# ind basis
prop_IBD_df %>% 
  ungroup() %>% 
  mutate(IID = as.character(IID)) %>% 
  filter(IID %chin% sample(as.character(unique(prop_IBD_df$IID)), 
                           100, replace = FALSE)) -> plot_ibd_df

col_pal <- rev(c(brewer.pal(7, "YlGnBu"),  "firebrick1")) ##A42820
#col_pal <- rev(c(brewer.pal(6, "YlGnBu"), "goldenrod", "firebrick1")) 
prop_IBD_df %>% 
  ungroup() %>% 
  # filter(class %in% c(1,2,4)) %>% 
  filter(IID %chin% sample(as.character(unique(prop_IBD_df$IID)), 500, replace = FALSE)) %>% 
  ggplot(aes(x=IID, y=prop_IBD, fill=length_class)) +
          geom_bar(stat="identity") +
         # scale_fill_brewer(palette = "YlGnBu", name = "ROH class (Mb)", direction = -1) +
          scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
          #facet_grid(.~ pop, space = 'free_x', scales = 'free_x', switch = 'x') +
          theme_simple() +
          theme(axis.text.x = element_blank()) +
          ylab("Proportion of genome in ROH") + 
          xlab("Individuals") + 
          scale_x_discrete(expand = c(0, 0)) + 
          scale_y_continuous(expand = c(0, 0)) +
          theme(legend.position = "bottom") -> p_roh2
p_roh
ggsave("figs/roh_classes_subset_inds.jpg", p_roh2, width = 15, height = 5)


#~~~ ROH density
hom_sum <- fread("output/ROH/roh_nofilt_ram.hom.summary")

hom_sum <- hom_sum %>%
  mutate(MB = BP / 1000000,
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


#~~ ROH for some indiividuals
all_roh <- roh_lengths %>% 
  group_by(IID) %>% 
  summarise(sum_roh = sum(KB)) %>% 
  ungroup() %>% 
  arrange(desc(sum_roh))

longest_roh <- all_roh %>% 
  top_n(10)
shortest_roh <- all_roh %>% 
  top_n(-10)

num_ind <- 20

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
col <- wes_palette("Darjeeling2")[c(1,2)]
col <- viridis(2, begin = 0.1, end = 0.9)
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
col <- c("#82687D", "#2C5475" )
col <- c("#d8b365", "#5ab4ac")
df %>% 
  filter(MB > 1) %>% 
  filter(CHR %in% 1:26) %>% 
  #mutate(CHR = factor(CHR, levels = as.character(1:26))) %>% 
  ggplot() +
        geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
                  alpha=0.8, fill = "#f7f7f7") +
        geom_hline(data = yax, aes(yintercept = yax), color = "grey", size = 0.2) +
        geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.5, 
                     fill = as.factor(CHR)),  col = "black", size = 0.1, alpha = 1) +
        scale_fill_manual(values = rep(col, 18)) +
        scale_color_manual(values = rep(col, 18)) +
        #scale_x_discrete(labels = as.character(c(1:20, 22, 24, 26))) + 
        scale_y_reverse(expand = c(0, 0)) +
        theme_clean() +
        facet_grid(~CHR, scales = 'free_x', space = 'free_x', switch = 'x') +
        theme(strip.placement = 'outside',
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.spacing = unit(0, "lines"),
              legend.position="none",
              plot.margin = unit(c(1, 1, 1, 1), "cm"), 
              axis.text.y = element_text(colour = "white"),
              axis.line.y = element_blank(),
              axis.title=element_text(size=15))+
        xlab("chromosome") +
        ylab("individuals") +
        ggtitle("ROH > 1Mb in the 10 most and 10 least inbred individuals") -> ROH_per_ind

ROH_per_ind
#dev.off()
ggsave("figs/roh_per_ind_1Mb.jpg", ROH_per_ind, width = 11, height = 5)


# test with ggchicklets
ggplot(df, aes(IID, KB, group = POS1)) +
  geom_chicklet() + 
  coord_flip()
