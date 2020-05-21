
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

