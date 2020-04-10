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
library(lme4)
library(broom.mixed)
library(INLA)
library(snpStats)
library(furrr)
library(brinla)
library(performance)
library(sjPlot)
library(ggeffects)
library(patchwork)
# data
# annual measures of traits and fitness
load("data/survival_mods_data.RData")

# Plot A: froh across age classes ----------------------------------------------
# number of individuals in each age class
max_age_df <- fitness_data %>% 
        group_by(id) %>% 
        mutate(max_age = max(age)) %>% 
        filter(age == max_age) %>% 
        mutate(age_class_0 = ifelse(age >= 0, 1, 0),
               age_class_1 = ifelse(age >= 1, 1, 0),
               age_class_2 = ifelse(age >= 2, 1, 0),
               age_class_3 = ifelse(age >= 3, 1, 0),
               age_class_4 = ifelse(age >= 4, 1, 0),
               age_class_5 = ifelse(age >= 5, 1, 0),
               age_class_6 = ifelse(age >= 6, 1, 0),
               age_class_7 = ifelse(age >= 7, 1, 0),
               age_class_8 = ifelse(age >= 8, 1, 0),
               age_class_9 = ifelse(age >= 9, 1, 0),
               age_class_10 =ifelse(age >= 10, 1, 0)
        ) %>% 
        ungroup() %>% 
        filter(!is.na(froh_all))# %>% 

ids_per_age <- map(0:10, function(x) {
        to_filt <- paste0("age_class_", x)
        max_age_df %>% 
                filter(get(to_filt) == 1) %>% 
                dplyr::select(id) %>% 
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
        mutate(age = str_replace(age, "_", " ")) %>% 
        mutate(age_num = as.numeric(str_replace(age, "age ", ""))) %>% 
        as.data.frame() %>% 
        mutate(age = factor(age)) %>% 
        mutate(age = fct_reorder(age, age_num))-> roh_plot2

roh_plot2 %>% 
        #ungroup() %>% 
        #sample_frac(0.1) %>% 
        ggplot(aes(x = FROH, y = age)) +
        geom_density_ridges(scale = 0.85,
                            jittered_points=TRUE, color = "#4c566a",  # # "#4c566a"  "#eceff4"
                            fill = "#eceff4",
                            point_shape = "|", point_size = 2.5, size = 0.4,
                            position = position_points_jitter(height = 0),
                            point_alpha = 1, point_color = "#4c566a", 
                            alpha = 0.9, bandwidth = 0.006) +
        # scale_fill_viridis() +
        #theme_ridges(font_family = "Avenir") + 
        theme_simple(axis_lines = FALSE, base_size = 14) +
        ylab("") + 
        xlab(expression(F[ROH]~across~age~cohorts)) +
        theme(strip.background = element_blank(),
              panel.grid.major = element_line(colour = "#d8dee9", size = 0.5),
              strip.text = element_text(vjust=2),
              panel.spacing.x = unit(1, "lines"),
              strip.text.x = element_text(margin = margin(1,0,0,0, "cm")),
              legend.position = "none") -> p_froh_across_ages #+
#labs(title = "Inbreeding across life in Soay sheep",
#    subtitle = "- individuals with lots of long ROHs rarely survive their first year")
p_froh_across_ages
# ggsave("figs/FROH_across_ages_col.jpg", plot = p_froh_across_ages, width = 4.5, height = 5.5)



# Plot B: effect sizes ---------------------------------------------------------
mod_inla <- readRDS("output/AS_mod_INLA.rds")

trans_link_to_dat <- function(pred, mod_inla) {
       trans_pred <-  inla.rmarginal(n = 10000, marginal = mod_inla$marginals.fixed[[pred]]) %>% 
                exp() 
       sum_pred <- c(mean = mean(trans_pred), quantile(trans_pred, probs = c(0.025, 0.975)))
       sum_pred %>%  t() %>% 
                as.data.frame() %>% 
                bind_cols() %>% 
                add_column(pred = pred, .before = 1)
}

fix_eff <- map_df(names(mod_inla$marginals.fixed), trans_link_to_dat, mod_inla) %>% 
        .[c(2,7,8), ] %>% 
        mutate(pred = factor(pred, levels = rev(pred)))
names(fix_eff) <- c("Predictor", "mean", "lower_CI", "upper_CI")

p_surv_mod <- ggplot(fix_eff, aes(mean, Predictor, xmax = upper_CI, xmin = lower_CI)) +
        geom_vline(xintercept = 1, linetype='dashed', colour =  "#4c566a", size = 0.3) +     # "#4c566a"  "#eceff4"
        geom_errorbarh(alpha = 1, height = 0.5,
                       size = 0.5) +
        geom_point(size = 3, shape = 21, col = "#4c566a", fill = "#eceff4", # "grey69"
                   alpha = 1, stroke = 0.7) + 
        # when only showing roh related effs
        scale_x_log10(breaks = c(0.1, 0.3, 1, 2, 5), limits = c(0.1, 5), labels = c( "0.1", "0.3", "1", "2", "5")) +
        scale_y_discrete(labels = rev(c(expression(F[ROH]), expression(F[ROH]:Age), expression(F[ROH]:Lamb)))) +
        theme_simple(axis_lines = TRUE, base_size = 14) +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_text(),
                axis.title.x = element_text(margin=margin(t=8))
        ) +
        xlab(expression(beta~and~95*"%"*~CI~(odds~of~survival))) -> p_forest


# Plot C marginal effects ------------------------------------------------------

mod_inla <- readRDS("output/AS_mod_INLA.rds")

# plot INLA marginal effects
fun <- function(...) {
        one <-  invlogit(Intercept + 
                                 df1$x1 * froh_all10_cent + 
                                 df1$x2 * age_cent + 
                                 df1$x3 * lamb1 + 
                                 df1$x4 * twin1 + 
                                 df1$x5 * sexM + 
                                 df1$x6 * `froh_all10_cent:lamb1` + 
                                 df1$x7 * `froh_all10_cent:age_cent`)
        
        return (list(one))
}

invlogit <- function(x) exp(x)/(1+exp(x))

froh <- seq(from = min(annual_survival$froh_all10_cent), to = (max(annual_survival$froh_all10_cent)), by = 0.1)
age <- c(-2.4, -1.4, 1.6, 4.6)
combined_df <- expand_grid(froh, age) %>% 
        mutate(lamb = ifelse(age == -2.4, 1, 0),
               twin = 0.15,
               sex = 0.4,
               frohxlamb = froh*lamb,
               frohxage = froh*age) 
names(combined_df) <- paste0("x", 1:7)

xx <- inla.posterior.sample(1000, mod_inla)
marg_means <- purrr::map(1:nrow(combined_df), function(x) {
        df1 <<- combined_df[x, ]
        out <- inla.posterior.sample.eval(fun, xx)
}) 


d <- marg_means %>% 
        map(as_tibble) %>% 
        bind_rows() %>% 
        pmap_df(function(...) {
                samp <- as.numeric(unlist(list(...)))
                c(mean = mean(samp), quantile(samp, probs = c(0.025, 0.975)))
        }) %>% 
        bind_cols(combined_df) %>% 
        .[, 1:5] %>% 
        setNames(c("prediction", "ci_lower", "ci_upper", "froh", "age")) %>% 
        mutate(age = as.factor(round(age + mean(annual_survival$age), 0)),
               froh = (froh + mean(annual_survival$froh_all10))/10)


#saveRDS(d, file = "output/AS_mod_INLA_predictions_for_plot.rds")
inla_preds <- readRDS("output/AS_mod_INLA_predictions_for_plot.rds") %>% 
        mutate(prediction = prediction * 100,
               ci_lower = ci_lower * 100,
               ci_upper = ci_upper * 100)

p_marginal_effs <- ggplot(inla_preds, aes(froh, prediction)) +
        geom_line(aes(color = age), size = 1.5) +
        geom_ribbon(aes(x=froh, ymin = ci_lower, ymax = ci_upper, fill = age, color = age),
                    alpha = 0.2, linetype = 2, size = 0.1)+
        scale_color_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        scale_fill_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        scale_y_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 25),
                           limits = c(0,100)) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 14) +
        theme(axis.line.y = element_blank(),
              legend.position = "top") +
        xlab(expression(F[ROH])) +
        ylab("Predicted\nsurvival probability %")


p_froh_across_ages + (p_forest / p_marginal_effs)

library(cowplot)
p_final <- cowplot::plot_grid(p_froh_across_ages, plot_grid(p_forest, p_marginal_effs, nrow = 2, rel_heights = c(1,1.5),
                                                            labels = c('B', 'C')),
                              labels = c('A', ''))
ggsave("figs/Fig2_inla.jpg", height = 6, width = 9)






