# Exploring ROH through plots

library(data.table)
library(RColorBrewer)
library(wesanderson)
library(readxl)
library(patchwork)
library(tidyverse)
library(zoo)
source("theme_simple.R")
library(ggridges)
library(viridis)
library(GGally)
library("naniar")
options(scipen=999)
library(ggchicklet)
library(windowscanr)
library(cowplot)
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

#~~ ROH
file_path <- "output/ROH/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path)

# survival data 
load("data/survival_mods_data.RData")
annual_survival <- fitness_data %>% 
  # filter na rows
  filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
  as.data.frame() 
ids_surv <- unique(as.character(annual_survival$id))

# subset roh to do all further plots for ids from survival analysis
roh_lengths <- roh_lengths %>% filter(IID %in% ids_surv)

# descriptive ROH statistics
num_roh_per_ind <- roh_lengths %>% group_by(IID) %>% tally() 
summary(num_roh_per_ind$n)

# GENOME SIZE NEEDS fixing (only autosomes)
froh <- roh_lengths %>%
        dplyr::group_by(IID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/autosomal_genome_size) %>% 
        mutate(FROH_cent = FROH - mean(FROH))
mean(froh$FROH)
range(froh$FROH)

# FROH / ROH across individuals
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
ggsave("figs/Sup_ROH_dist.jpg", p_roh_dist, width = 6, height = 2.7)


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
 
# ROH class distribution across individuals
col_pal <- rev(c(brewer.pal(7, "YlGnBu"))) ##A42820
col_pal <- viridis(7)
prop_IBD_df
set.seed(17)
library(gghalves)
prop_IBD_df_with_0 %>% 
  mutate(prop_IBD = prop_IBD * 100) %>% 
  ggplot(aes(length_class, prop_IBD, fill = length_class)) +
  geom_half_point(side = "l", shape = 21, alpha = 0.5, stroke = 0.1, size =2,
                  transformation_params = list(height = 0, width = 1.3, seed = 1)) +
  geom_half_boxplot(side = "r", outlier.color = NA,
                    width = 0.6, lwd = 0.3, color = "black",
                    alpha = 0.8) +
  theme_simple(axis_lines = TRUE, base_size = 13) +
  ylab("% genome") +
  scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black")) + 
  xlab("ROH classes") -> p_roh_classes

# prop_IBD_df_with_0 %>% 
#   mutate(prop_IBD = prop_IBD * 100) %>% 
#   ggplot(aes(length_class, prop_IBD, fill = length_class)) +
#         geom_jitter(alpha = 0.6, width = 0.2, shape = 21, color = "black", stroke = 0.1) +
#         geom_boxplot(outlier.shape = NA, color = "#4c566a", alpha = 0.01, width = 0.6, ##4c566a
#                      lwd = 0.4) +
#         theme_simple(axis_lines = TRUE) +
#         ylab("% genome") +
#          scale_fill_manual(values = col_pal, name = "ROH class (Mb)") +
#         theme(legend.position = "none",
#               axis.ticks.x = element_blank()) + 
#         xlab("ROH classes") -> p_roh_classes
# p_roh_classes

# ggsave("figs/roh_classes_boxplots.jpg", p_roh, width = 9, height = 3.5)


#ggsave("figs/roh_classes_subset_inds.jpg", p_roh2, width = 15, height = 5)
p_roh_final <- p_roh2 / p_roh + plot_layout(guides = "collect") +
              plot_annotation(tag_levels = 'A') #& 
p_roh_final
#ggsave("figs/roh_classes.jpg", p_roh_final, width = 8, height = 5)


#~~ ROH for some indiividuals
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
#col <- viridis(6, option = "D")[c(2,6)]
chr_names <- as.character(1:26)
names(chr_names) <- as.character(1:26)
chr_names[c(11, 13, 15, 17, 19, 20, 21, 23,24, 25)] <- ""

df %>% 
  filter(MB > 5) %>% 
  filter(CHR %in% 1:26) %>% 
  #mutate(CHR = factor(CHR, levels = as.character(1:26))) %>% 
  ggplot() +
  geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
            alpha=0.4, fill = "#eceff4") + # "#f7f7f7"
  geom_hline(data = yax, aes(yintercept = yax), color = "#d8dee9", size = 0.3) +
  geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.9, 
                fill = as.factor(CHR)),  col = "black", size = 0.1, alpha = 1) +
  scale_fill_manual(values = rep(col, 18)) +
  scale_color_manual(values = rep(col, 18)) +
 # scale_x_discrete(labels = as.character(c(1:20, 22, 24, 26))) + 
  scale_y_reverse(expand = c(0, 0)) +
  theme_simple(axis_lines = TRUE, base_size = 12) +
  facet_grid(~CHR,scales = 'free_x', space = 'free_x', switch = 'x',
             labeller = as_labeller(chr_names)) +
  ggtitle("ROH in the 7 most and least inbred individuals") +
  theme(#strip.placement = 'outside',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.line.x = element_blank(),
        legend.position="none",
        #plot.margin = unit(c(1, 1, 1, 1), "cm"), 
        axis.text.y = element_text(colour = "white"),
        axis.line.y = element_blank(),
        axis.title=element_text(size=15))+
  xlab("chromosome") +
  ylab("individuals") -> ROH_per_ind
#ggtitle("ROH > 1Mb in the 10 most and 10 least inbred individuals") 
ROH_per_ind
#dev.off()
#ggsave("figs/roh_per_ind_5Mb.jpg", ROH_per_ind, width = 12, height = 5)

#~~~ ROH density
hom_sum <- fread("output/ROH/roh_nofilt_ram.hom.summary")

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
                       funs = c("mean", "var"))

# ROH sharing 1
running_roh %>% 
  mutate(UNAFF_mean = UNAFF_mean/7691) %>% 
  filter(!is.na(UNAFF_mean))-> running_roh_p
  #filter(CHR == 1) %>% 
library(scales)
fill_cols <- viridis(20, option = "D")
qn <- scales::rescale(quantile(running_roh_p$UNAFF_mean,
                       probs=seq(0, 1, length.out=length(fill_cols))))

p1 <- ggplot(running_roh_p, aes(x = win_start, y = 0.5, fill = UNAFF_mean)) + 
  geom_tile(color = "black", size = 0.01) +
  theme_simple(base_size = 15) + 
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), labels = c(0, 100, 200, 300))+
  ylab("Chromosome") +
  scale_fill_gradientn("Proportion of Sheep with ROH",
                       colors = rev(fill_cols), values = qn,
                       breaks = c(0.2, 0.5, 0.8)) +
  facet_grid(CHR~., switch="both") +
  xlab("Position in mb") +
  theme(panel.spacing.y=unit(0.1, "lines"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.ticks.x = element_line(),
        #axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.605,0.15),
        legend.direction = "horizontal",
        strip.text.y = element_text(size = 10, angle = 180),
        axis.line.x = element_blank()) +
  guides(fill = guide_colourbar(title.position = "bottom" ,
                                barwidth = 10.95, barheight = 0.5))

p2 <- ggplot(running_roh_p, aes(UNAFF_mean, "test", fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3) +
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8)) +
  theme_void() +
  scale_fill_gradientn("Proportion of Sheep with ROH",
                       colors = rev(fill_cols), values = qn,
                       breaks = c(0.2, 0.5, 0.8)) +
  theme(legend.position = "none")
        #axis.ticks.x = element_line(size = 1)) 

p2
layout <- c(
  area(t = 1, l = 1, b = 5, r = 5),
  area(t = 4, l = 3, b = 4, r = 4)
)
p_final <- p1 + p2 + plot_layout(design = layout)
p_final
ggsave("figs/roh_genome.jpg", p_final, width = 8, height = 6)

p_roh_comb <- plot_grid(p_roh_classes, ROH_per_ind, nrow = 2, rel_heights = c(2.2,3, 5.2))
p_roh_comb
ggsave("figs/roh_comb.jpg", p_roh_comb, width = 8, height = 6)

plot_grid(p_roh_classes, NULL, ROH_per_ind, p1, ncol = 2, rel_heights = 1, 1.5, 2.5)


fitness_data




# ggplot(aes(x = win_start, y = 0.5, fill = UNAFF_mean)) + 
#   geom_tile(color = "black", size = 0.01) +
#   theme_simple(base_size = 12) + 
#   scale_y_continuous(expand = c(0,0))+
#   ylab("Chromosome") +
#   scale_fill_gradientn("Proportion of\nSheep with ROH",
#                        colors = c("#053061",rep("#f7f7f7", 3), "#b2182b"),
#                        values = c(0, 0.1,0.3, 0.5, 0.7, 0.8, 1),
#                        breaks = c(0.2, 0.5, 0.8)) +
#   facet_grid(CHR~., switch="both") +
#   xlab("Position in KB") +
#   theme(panel.spacing.y=unit(0.1, "lines"),
#         axis.text.y = element_blank(),
#         #axis.title.y = element_blank(),
#         axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         legend.position = c(0.6,0.2),
#         legend.direction = "horizontal",
#         strip.text.y = element_text(size = 8, angle = 180),
#         axis.line.x = element_blank()) +
#   guides(fill = guide_colourbar(title.position = "top" ,
#                                 barwidth = 10, barheight = 0.5)) -> p_roh_across_genome
# p_roh_across_genome
# 
# ggsave("figs/roh_across_genome.jpg", p_roh_across_genome, width = 8, height = 6)


# ROH sharing 2
# p_running_roh <- running_roh %>% 
#   filter(CHR %in% 1:26) %>% 
#   ggplot(aes(win_start, UNAFF_mean)) +
#   geom_line() +
#   scale_y_continuous(breaks = c(369, 1844, 5531, 7005), labels = c("5%", "25%", "75%", "95%"),
#                      limits = c(0,7374)) +
#   theme_simple() + 
#   ggtitle("Average ROH across the genome of 7374 Soay sheep") +
#   xlab("Window start (in Base Pairs) of 200 BP window") +
#   ylab("Individuals with ROH") + 
#   annotate("rect", xmin=-Inf, xmax=Inf, ymin=1844, ymax=5531, alpha=0.1) +
#   geom_hline(yintercept = 369, size = 0.3) +
#   geom_hline(yintercept = 7005, size = 0.3) +
#   facet_grid(CHR~.) 
# p_running_roh 
# 
# ggsave("figs/roh_across_genome.jpg", p_running_roh, width = 8, height = 12)



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
