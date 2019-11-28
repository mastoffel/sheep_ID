library(gganimate)
library(gapminder)
library(ggplot2)
library(tidyverse)
library(data.table)
library(magrittr)
library(rlang)
source("theme_clean.R")
# annual measures of traits and fitness
load("data/fitness_roh_df.RData")
library(windowscanr)
library(wesanderson)

#fitness_data %<>% 
#        filter(ID %in% sample(unique(ID), 3500))

# number of individuals in each age class
max_age_df <- fitness_data %>% 
        group_by(ID) %>% 
        mutate(max_age = max(Age)) %>% 
        filter(Age == max_age) %>% 
        # filter(BIRTHYEAR %in% c(2017,2018)) %>% 
        # filter(is.na(Survival)) %>% 
        mutate(age_class_0 = ifelse(Age >= 0, 1, 0),
               age_class_1 = ifelse(Age >= 1, 1, 0),
               age_class_2 = ifelse(Age >= 2, 1, 0),
               age_class_3 = ifelse(Age >= 3, 1, 0),
               age_class_4 = ifelse(Age >= 4, 1, 0),
               age_class_5 = ifelse(Age >= 5, 1, 0),
               age_class_6 = ifelse(Age >= 6, 1, 0),
               age_class_7 = ifelse(Age >= 7, 1, 0),
               age_class_8 = ifelse(Age >= 8, 1, 0),
               age_class_9 = ifelse(Age >= 9, 1, 0),
               age_class_10 = ifelse(Age >= 10, 1, 0)
        ) %>% 
        ungroup() %>% 
        filter(!is.na(FROH_all))# %>% 
# filter(FROH_all < 0.375)


# number of individuals with genotype data per age class
# num_ind_per_age <- max_age_df %>% 
#         summarise_at(vars(age_class_0:age_class_10), sum)

ids_per_age <- map(0:10, function(x) {
        to_filt <- paste0("age_class_", x)
        max_age_df %>% 
                filter(get(to_filt) == 1) %>% 
                dplyr::select(ID) %>% 
                unlist()
})
names(ids_per_age) <- paste0("age_", c(0:10))

num_ind_per_age <- unlist(map(ids_per_age, length))

# load ROH info
file_path <- "output/ROH/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path) 

# calculate length classes
calc_froh_classes <- function(roh_crit, roh_lengths) {
        
        roh_filt <- dplyr::case_when(
                roh_crit == "short"  ~ expr(KB < 1221),
                roh_crit == "medium" ~ expr((KB > 1221)&(KB < 4885)),
                roh_crit == "long"   ~ expr(KB > 4885),
                roh_crit == "all" ~ expr(KB > 0)
        )
        
        roh_lengths %>%
                dplyr::group_by(IID) %>%
                #filter({{ roh_filt }}) %>% 
                filter(!!roh_filt) %>% 
                dplyr::summarise(KBSUM = sum(KB)) %>% 
                mutate(FROH = KBSUM / 2869914) %>% 
                dplyr::select(IID, FROH) %>% 
                rename(ID = IID, 
                       !! paste0("FROH_", roh_crit) := FROH)
                      # !! paste0("ROH_", roh_crit) := ROH)
}

# proportion of ROH length classes in each genome. Individuals which
# do not have long ROH have 0 for this class.
froh <- purrr::map(c("short", "medium", "long", "all"), calc_froh_classes,  roh_lengths) %>% 
        purrr::reduce(left_join, by = "ID") %>% 
        replace_na(list(FROH_long = 0)) %>% 
        pivot_longer(cols = starts_with("FROH"), names_to = "roh_class", values_to = "prop")


roh_per_age <- map(c(1:length(ids_per_age)), function(x) {
        out <- froh %>% 
          filter(ID %in% ids_per_age[[x]]) %>% 
          mutate(age = paste0("age_", x))
        out
})

roh_plot <- roh_per_age %>% 
        bind_rows() %>% 
        mutate(age = fct_inorder(age))
        
p1 <- ggplot(roh_plot, aes(prop, fill = roh_class)) +
        #facet_wrap(age~roh_class, scales = "free_y" ,nrow = 4) + # scales = "free",
        geom_histogram(bins = 100, color = "grey", lwd = 0.2) +
        scale_y_sqrt() + 
        scale_fill_brewer(palette = "YlGnBu", direction = -1) + 
        theme_clean() +
        theme(axis.text.y = element_blank(),
              axis.line.y = element_blank()) +
        transition_states(age,
                          transition_length = 2,
                          state_length = c(7,6,5,4,3,2,1,1,1,1,1)) + 
        ease_aes('linear') + 
        labs(title = "{closest_state}") +
        enter_fade() +
        exit_fade() +
        view_follow() 
p1      

animate(p1, fps = 5, duration = 15)
p1
anim_save("figs/roh_dist_across_life.gif", animation = p1)



# another animated plot

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
                                 length_Mb < 1.221366094 & length_Mb >= 0.6 ~ 64 )) %>%
                                # length_Mb < 0.610683047 & length_Mb >= 0.30 ~ 128)  ) %>% # 0.610683047
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
                #class == 128 ~ "0.6-0.3 (128G)"
        )) %>% 
        mutate(length_class = fct_reorder(length_class, class)) %>% 
        mutate(IID = as.character(IID)) %>% 
        group_by(IID, class, length_class) %>%
        dplyr::summarise(prop_IBD = sum(length_Mb / 2869),
                         sum_IBD = sum(length_Mb)) 
# add IBD of non-ROH snps if wanted
#  bind_rows(homs) 

prop_IBD_df_with_0 <- prop_IBD_df %>% 
        # add missing length classes as 0
        ungroup() %>% 
        tidyr::complete(length_class, nesting(IID)) %>% 
        mutate(class = ifelse(is.na(class), length_class, class)) %>% 
        mutate(prop_IBD = ifelse(is.na(prop_IBD), 0, prop_IBD),
               sum_IBD = ifelse(is.na(sum_IBD), 0, sum_IBD))


# ROH class distributions as prop. 
prop_IBD_across_time <- prop_IBD_df_with_0 %>% 
        rename(ID = IID) %>% 
        mutate(ID = as.numeric(ID)) %>% 
        left_join(dplyr::select(max_age_df, ID, starts_with("age_class")), by = "ID") %>% 
        pivot_longer(cols = starts_with("age_class"), names_to = "age_class") %>% 
        filter(value == 1) %>% 
        mutate(age = as.numeric(str_replace(age_class, "age_class_", "")))

# 
prop_IBD_across_time %>% 
  filter(length_class == "> 39 (1G)") %>% 
ggplot(aes(sum_IBD, fill= as.factor(class))) + 
      geom_histogram(bins = 1000) +
      scale_y_log10() + 
      facet_grid(class~age)

library(ggbeeswarm)
p1 <- prop_IBD_across_time %>% 
        filter(!is.na(length_class)) %>% 
        group_by(ID) %>% 
        sample_frac(1) %>% 
        ungroup() %>% 
ggplot(aes(length_class, prop_IBD)) +
       #geom_jitter(alpha = 0.3, width = 0.2, size = 2) +
        geom_quasirandom(alpha = 0.2, width = 0.3, size = 4) +
        #geom_beeswarm(alpha = 0.3, cex = 0.7, size = 2) +
        geom_boxplot(outlier.shape = NA, color = "lightgrey", alpha = 0.7) +
        theme_clean() +
        theme(axis.text = element_text(size = 15),
              axis.title = element_text(size = 18),
              plot.title=element_text(size=20,face="bold")) +
        #theme(axis.ticks.x = element_line(colour = "#cccccc", size = 0.3)) +
      #  xlab("ROH class (Mb)") + ylab("Proportion of genome in ROH") +
        labs(title = 'Age: {closest_state}', x = "ROH class (Mb)", y = "Proportion of genome in ROH") +
        transition_states(factor(age)) +
        ease_aes('cubic-in-out') 

p_rend <- animate(p1, fps = 20, duration = 5, start_pause = 20, height = 560, width =1120)
anim_save("figs/roh_dist_across_life.gif", animation = p_rend)




# FROH across age cohorts 

all_roh_per_age <- purrr::map(1:length(ids_per_age), function(x) {
  out <- roh_lengths %>% 
    rename(ID = IID) %>% 
    filter(ID %in% ids_per_age[[x]]) %>% 
    mutate(age = paste0("age_", x-1))
  out
})

roh_plot <- all_roh_per_age %>% 
  bind_rows() %>% 
  mutate(age = fct_inorder(age)) %>% 
  mutate(MB = KB/1000) %>% 
  mutate(class = dplyr::case_when(
    KB < 1221 ~ "short",
    (KB > 1221)&(KB < 4885) ~ "medium",
    KB > 4885 ~ "long"
  ))

# roh_plot %>% 
#   sample_frac(0.1) %>% 
#   ggplot(aes(MB)) + geom_histogram(bins = 100) +
#   facet_wrap(age~class, scales = "free_y", ncol = 3)

library(ggridges)
library(viridis)

roh_plot %>% 
  #filter(MB != 0) %>% 
  #sample_frac(1) %>% 
  group_by(ID, class, age) %>% 
  summarise(MB_mean = mean(MB),
            MB_sum = sum(MB)) -> roh_plot1

library(ghibli)
library(wesanderson)
library(viridis)
library(ggthemes)
library(nord)
?nord
#cols <- fivethirtyeight_pal()(3)
cols <- nord("victory_bonds", 3)
p <- ggplot(roh_plot1, aes(x = MB_sum, y = age, fill =class, point_color = class)) +
  geom_density_ridges(scale = 0.9,
                               jittered_points=TRUE,
                               point_shape = "|", point_size = 3, size = 0.4,
                               position = position_points_jitter(height = 0),
                               point_alpha = 1,
                               alpha = 0.5)+
 # scale_fill_viridis() +
  scale_fill_manual(values =cols ) + 
  scale_discrete_manual("point_colour", values = cols ) + 
  theme_ridges(font_family = "Avenir") + 
  facet_wrap(~class, scales = "free_x") +
  ylab("") +
  xlab("Sum of ROH per individual in Mb") +
  theme(strip.background = element_blank(),
        strip.text = element_text(vjust=2),
        panel.spacing.x = unit(1, "lines"),
        strip.text.x = element_text(margin = margin(1,0,0,0, "cm")),
        legend.position = "none") +
  labs(title = "ROH across life in Soay sheep",
       subtitle = "- individuals with lots of long ROHs rarely survive their first year")
  

ggsave("figs/roh_across_life_ridges_sum.jpg", plot = p, height = 7, width = 8)        
  


# FROH across life
roh_plot %>% 
  #filter(MB != 0) %>% 
  #sample_frac(1) %>% 
  group_by(ID, age) %>% 
  summarise(MB_mean = mean(MB),
            MB_sum = sum(MB)) %>% 
  mutate(FROH = MB_sum/2869.898) %>% 
  mutate(age = str_replace(age, "_", " ")) -> roh_plot2

roh_plot2 %>% 
  #ungroup() %>% 
  #sample_frac(0.1) %>% 
  ggplot(aes(x = FROH, y = fct_inorder(age))) +
  geom_density_ridges(scale = 0.85,
                      jittered_points=TRUE, color = "#465881", 
                      fill = "#c9d1d3",
                      point_shape = "|", point_size = 2.5, size = 0.4,
                      position = position_points_jitter(height = 0),
                      point_alpha = 1, point_color = "#1b2a49", 
                      alpha = 0.9, bandwidth = 0.006) +
  # scale_fill_viridis() +
  theme_ridges(font_family = "Avenir") + 
  ylab("") +
  xlab("FROH across age cohorts") +
  theme(strip.background = element_blank(),
        strip.text = element_text(vjust=2),
        panel.spacing.x = unit(1, "lines"),
        strip.text.x = element_text(margin = margin(1,0,0,0, "cm")),
        legend.position = "none") -> p2 #+
  #labs(title = "Inbreeding across life in Soay sheep",
   #    subtitle = "- individuals with lots of long ROHs rarely survive their first year")
p2

ggsave("figs/FROH_across_ages_col.jpg", plot = p2, width = 5, height = 6.5)
