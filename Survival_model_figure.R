library(ghibli)
library(wesanderson)
library(viridis)
library(ggthemes)
library(nord)
library(ggridges)
library(viridis)
library(tidyverse)
library(data.table)
library(magrittr)
library(rlang)
source("theme_simple.R")
library(patchwork)
# annual measures of traits and fitness
load("data/fitness_roh_df.RData")


# PLOT 1: FROH ACROSS AGE CLASSES
# number of individuals in each age class
max_age_df <- fitness_data %>% 
        group_by(ID) %>% 
        mutate(max_age = max(Age)) %>% 
        filter(Age == max_age) %>% 
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

# FROH across life
roh_plot %>% 
        group_by(ID, age) %>% 
        summarise(MB_mean = mean(MB),
                  MB_sum = sum(MB)) %>% 
        mutate(FROH = MB_sum/2655) %>% 
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
        xlab(expression(F[ROH]~across~age~cohorts)) +
        theme(strip.background = element_blank(),
              strip.text = element_text(vjust=2),
              panel.spacing.x = unit(1, "lines"),
              strip.text.x = element_text(margin = margin(1,0,0,0, "cm")),
              legend.position = "none") -> p2 #+
#labs(title = "Inbreeding across life in Soay sheep",
#    subtitle = "- individuals with lots of long ROHs rarely survive their first year")
p2
ggsave("figs/FROH_across_ages_col.jpg", plot = p2, width = 4.5, height = 5.5)

## plot 2

mod_inla <- readRDS("output/inla_survival_model_full.rds")

trans_link_to_dat <- function(pred, mod_inla) {
        inla.rmarginal(n = 10000, marginal = mod_inla$marginals.fixed[[pred]]) %>% 
                exp() %>% 
                quantile(probs = c(0.025,0.5, 0.975))  %>% 
                t() %>% 
                as.data.frame() %>% 
                bind_cols() %>% 
                add_column(pred = pred, .before = 1)
}

fix_eff <- map_df(names(mod_inla$marginals.fixed), trans_link_to_dat, mod_inla) %>% 
        .[c(2,4,5,6,8), ] %>% 
        filter(!(pred == "(Intercept)")) %>% 
        mutate(pred = fct_inorder(as.factor(pred)))
names(fix_eff) <- c("Predictor", "lower_CI", "median", "upper_CI")

p_surv_mod <- ggplot(fix_eff, aes(median, Predictor, xmax = upper_CI, xmin = lower_CI)) +
        geom_vline(xintercept = 1, linetype='dashed', colour =  "#4c566a", size = 0.7) +
        geom_errorbarh(alpha = 1, height = 0.5,
                       size = 0.5) +
        geom_point(size = 3, shape = 21, col = "#4c566a", fill = "#eceff4", # "grey69"
                   alpha = 1, stroke = 0.7) + 
        scale_y_discrete(limits = rev(levels(fix_eff$Predictor)),
                         #labels = rev(c("FROH", "Age", "Sex(Male)", "Twin", "FROH:Age"))) +
                         labels = rev(c(expression(F[ROH]), "Age", expression(Sex[Male]), "Twin", expression(F[ROH]:Age)))) + 
        #labels = rev(c("Froh", "Lamb", "Age", "Sex(Male)", "Twin", "Froh:lamb", "Froh:age"))) +
        scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1, 1.25)) +
        theme_simple(axis_lines = TRUE) +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.line.x = element_line(colour = "black",size = 0.8),
                axis.ticks.x = element_line(size = 0.8),
                axis.text = element_text(),
                axis.title.x = element_text(margin=margin(t=8))
        ) +
        xlab(expression(beta~and~95*"%"*~CI~(odds~of~survival)))
p_surv_mod
#ggsave("figs/surv_mod_fixs.jpg", p_surv_mod, width = 3.5, height = 3)


p2 + p_surv_mod

