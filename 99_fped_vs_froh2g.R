library(tidyverse)
library(data.table)
library(optiSel)
library(viridis)
library(gghalves)
source("theme_simple.R")

load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")

# roh data
file_path <- "output/ROH/roh.hom"
roh_lengths <- fread(file_path) 

# Chr lengths
chr_data <- read_delim("data/chromosome_info_oar31.txt", delim = "\t") %>% 
        rename(size_BP = Length,
               CHR = Part) %>% 
        mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
        .[2:27, ] %>% 
        summarise(sum_KB = sum(size_KB)) %>% 
        as.numeric()

# fped
ped <- sheep_ped[, c(1,3,2)]
ped_fin <- prePed(ped)
ID <- pedInbreeding(ped_fin) %>% as_tibble() %>% rename(id = Indiv, fped = Inbr)
ggplot(ID, aes(fped)) +
        geom_histogram(bins = 500) +
        scale_y_log10() 

#froh classes
length_dist <- data.frame(g = c(1, 2,2^2, 2^3, 2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13)) %>%
        mutate(ROH_length_cM = 100 / (2*g)) %>% 
        mutate(ROH_length_Mb = ROH_length_cM * 0.7816743)

prop_IBD_df <- roh_lengths %>%
        mutate(length_Mb = KB/1000) %>%
        mutate(class = case_when(#length_Mb >= 39.083715000 ~ 1,
                # length_Mb < 39.083715000 & length_Mb >= 19.541857500 ~ 2,
                length_Mb >= 19.541857500 ~ 2,
                length_Mb < 19.541857500 & length_Mb >= 9.770928750 ~ 4,
                #length_Mb < 9.770928750 & length_Mb >= 6.513952500 ~ 6,
                length_Mb < 9.770928750& length_Mb >= 4.885464375 ~ 8,
                # length_Mb < 4.885464375 & length_Mb >= 3.908371500 ~ 10,
                length_Mb < 4.885464375 & length_Mb >= 2.442732188 ~ 16,
                length_Mb < 2.442732188 & length_Mb >= 1.2 ~ 32)) %>% # 0.610683047 1.221366094
        mutate(length_class = case_when(
                #class == 1 ~ ">39 (1G)",
                #class == 2 ~ ">19.5-39 (2g)",
                class == 2 ~ ">19.5 (2g)",
                class == 4 ~ "9.8-19.5 (2-4g)",
                #class == 6 ~ "6.5-9.7 (6G)",
                class == 8 ~ "4.9-9.8 (4-8g)",
                # class == 10 ~ "3.9-4.9 (10G",
                class == 16 ~ "2.4-4.9 (8-16g)",
                class == 32 ~ "1.2-2.4 (16-32g)"
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


f_df <- prop_IBD_df_with_0 %>% 
       # filter(length_class == ">19.5 (2g)") %>% 
        group_by(IID) %>% 
        #summarise(mean_froh = mean(prop_IBD)) %>% 
        rename(id = IID) %>% 
        left_join(ID)


f_df <- f_df %>% 
        mutate(fped_class = case_when(
                #length_class == ">19.5 (2g)" & (fped >= 0 & fped < 0.075) ~ "<0.075",
                #length_class == ">19.5 (2g)" & (fped >= 0.075 & fped < 0.175) ~ "0.075-0.175",
                #length_class == ">19.5 (2g)" & (fped >= 0.175) ~ ">0.175",
                
                length_class == ">19.5 (2g)" & (fped >= 0 & fped < 0.1) ~ "<0.1",
                length_class == ">19.5 (2g)" & (fped >= 0.1 & fped < 0.2) ~ "0.1-0.2",
                length_class == ">19.5 (2g)" & (fped >= 0.2) ~ ">0.2",
                
                # length_class == ">19.5 (2g)" & (fped > 0.075 & fped < 0.175) ~ "0.125 +- 0.05",
                # length_class == ">19.5 (2g)" & (fped >= 0.175 & fped < 0.325) ~ "0.25 +- 0.075",
               # length_class == ">19.5 (2g)" & (fped > 0.20) ~ ">0.2",
                TRUE ~ "<0.1"
        )) %>% 
        mutate(fped_class = as.factor(fped_class))

# f_df <- f_df %>% 
#         mutate(fped_class = case_when(
#                 (fped > 0.075 & fped < 0.175) ~ 0.125,
#                 (fped > 0.20 & fped < 0.3) ~ 0.25,
#                 TRUE ~ 0
#         )) %>% 
#         mutate(fped_class = as.factor(fped_class))

col_pal <- plasma(6)
#col_pal <- paste0("#", (c("21295c","204683","1763a1","0a96d6","65bee2","BAE2F2")))
col_pal <- paste0("#", (c("D8DEE9", "9F6976","432371")))
#col_pal <- paste0("#", c("432371","714674","9f6976","cc8b79","e39d7a","faae7b"))
p_roh_length <- f_df %>% 
        mutate(prop_IBD = prop_IBD * 100) %>% 
        ggplot(aes(length_class, prop_IBD, fill = fped_class)) +
        geom_half_point(side = "l", shape = 21, alpha = 0.6, stroke = 0.1, size =2, color = "#4c566a",
                        transformation_params = list(height = 0, width = 0.5, seed = 1)) +
        geom_half_boxplot(side = "r", outlier.color = NA,
                          width = 0.6, lwd = 0.3, color = "black",
                          alpha = 0.8) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13,
                     base_family = "Helvetica") +
        ylab("% genome") +
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
        scale_fill_manual(values = col_pal, name = "Fped") +
        theme(#legend.position = "none",
              plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
              #axis.ticks.x = element_blank(),
              axis.title=element_text(size = rel(1.1)), 
              axis.text = element_text(color = "black")) + 
        xlab("ROH length class in Mb (~ generations to MRCA)") 
ggsave("figs/Sup_Fig1B_by_fped.jpg", p_roh_length, width = 8, height = 3)

f_df
