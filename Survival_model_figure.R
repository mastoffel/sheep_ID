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
              legend.position = "none") -> p_froh_across_ages #+
#labs(title = "Inbreeding across life in Soay sheep",
#    subtitle = "- individuals with lots of long ROHs rarely survive their first year")
p_froh_across_ages
# ggsave("figs/FROH_across_ages_col.jpg", plot = p_froh_across_ages, width = 4.5, height = 5.5)


# age interaction and lamb
mod <- readRDS("output/survival_mod_full_lme4")

df1 <- get_model_data(mod, type = "eff", 
                      terms = c("froh_all10_cent [all]", "age_cent[-1.4, 2.6, 5.6]", "lamb")) %>% 
        as_tibble() %>% 
        filter(facet == 0)
df2 <- get_model_data(mod, type = "eff", 
                      terms = c("froh_all10_cent [all]",  "age_cent[-2.4, -1.4]", "lamb")) %>% 
        as_tibble() %>% 
        filter(facet == 1) %>% 
        filter(group == -2.4)
df_full <- bind_rows(df1, df2) %>% 
        mutate(group = as.character(as.numeric(group) + 2.4),
               x = ((x + mean(annual_survival$froh_all10))/10),
               predicted = predicted * 100,
               conf.low = conf.low * 100,
               conf.high = conf.high * 100)

p_marginal_effs <- ggplot() +
        #geom_jitter(data = annual_survival, aes(froh_all10_cent, y = -0.1), height = 0.05,
        #            alpha = 0.2, size = 3) +
        #geom_point(data = df_full, aes(x = froh_all10_cent, fit, color = age_cent)) +
        geom_line(data = df_full, aes(x = x, predicted, color = group), size = 1.5) +
        geom_ribbon(data= df_full, aes(x=x, ymin=conf.low, ymax=conf.high, fill = group, color = group), alpha= 0.2,
                    linetype = 2, size = 0.1) +
        scale_color_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        scale_fill_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        scale_y_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 100, by = 25),
                           limits = c(0,100)) +
        theme_simple(axis_lines = TRUE, line_size = 0.3) +
        theme(axis.line.y = element_blank()) +
        xlab(expression(F[ROH])) +
        ylab("Predicted survival probability %")
p_marginal_effs  

get_model_data(mod, type = "est") %>% 
        as_tibble() %>% 
        ggplot(aes(estimate, term, xmax = conf.high, xmin = conf.low)) +
        geom_vline(xintercept = 1, linetype='dashed', colour =  "#4c566a", size = 0.3) +
        geom_errorbarh(alpha = 1, height = 0.5,
                       size = 0.5) +
        geom_point(size = 3, shape = 21, col = "#4c566a", fill = "#eceff4", # "grey69"
                   alpha = 1, stroke = 0.7) + 
        scale_x_log10(breaks = c(0.1, 0.3, 0.5, 1, 2), limits = c(0.028, 4), labels = c( "0.1", "0.3", "0.5", "1", "2")) +
        scale_y_discrete(labels = rev(c(expression(F[ROH]), "Age", "Lamb", expression(Sex[Male]), "Twin", expression(F[ROH]:Age), expression(F[ROH]:Lamb)))) +
        theme_simple(axis_lines = TRUE) +
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



p_froh_across_ages + (p_forest / p_marginal_effs)






## ALTERNATIVE PLOT 2 WITH INLA

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


p_froh_across_ages + p_surv_mod

# model predictions, third plot
# Run survival models 

#load("data/fitness_roh_df.RData") # formerly
load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")

# roh data
file_path <- "data/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Annual survival~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# survival data preprocessing
annual_survival <- fitness_data %>% 
        # filter na rows
        filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
        mutate(age_cent = age - mean(age, na.rm = TRUE),
               age_cent2 = age_cent^2,
               age_std = as.numeric(scale(age)),
               age_std2 = age_std^2,
               # times 10 to estimate a 10% percent increase
               froh_all10 = froh_all * 10,
               froh_all10_cent = froh_all10 - mean(froh_all10, na.rm = TRUE),
               lamb = ifelse(age == 0, 1, 0),
               lamb_cent = lamb - mean(lamb, na.rm = TRUE),
               lamb = as.factor(lamb)) %>% 
        as.data.frame() 


# a few lme4 trials
# time saver function for modeling
nlopt <- function(par, fn, lower, upper, control) {
        .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
                                  opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                              maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
        list(par = res$solution,
             fval = res$objective,
             conv = if (res$status > 0) 0 else res$status,
             message = res$message
        )
}

mod7 <- glmer(survival ~ froh_all10_cent * (age_cent + lamb) + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
plot_model(mod7)
check_model(mod7)
performance(mod7)

mod8 <- glmer(survival ~ froh_all10_cent * age_cent + lamb + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
              family = binomial, data = annual_survival,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
performance(mod8)
check_collinearity(mod8)
tidy(mod8, conf.int = TRUE)

plot_model(mod8, type = "eff", terms = c("froh_all10_cent [all]", "lamb"))

roh_eff_lamb <- ggeffect(mod8, type = "eff", terms = c("froh_all10_cent [all]", "lamb [1]")) %>% 
        as_tibble() #%>% 
#filter(group == 1)
summary(roh_eff_lamb)
roh_eff_adult <- ggeffect(mod8, type = "eff", 
                          terms = c("froh_all10_cent [all]", "age_cent[-1.4, 3.6, 7.6]")) %>%  # -1.4, 1.6, 3.6, 5.6, 7.6
        as_tibble() 
summary(roh_eff_adult)
p_lamb <- ggplot(roh_eff_lamb, aes(x, predicted)) +
        geom_line(size = 1) +
        #geom_point(data = annual_survival %>% filter(age == 0), 
        #           aes(x = froh_all10_cent, y = survival),
        #           size = 4, alpha = 0.1) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                    alpha = 0.1, linetype = 2, color = "grey") +
        theme_simple() +
        scale_y_continuous(labels = paste0(seq(from = 0, to = 100, by = 25), "%"), limits = c(0,1)) +
        scale_x_continuous(labels = c("0.24", "0.34", "0.44"), breaks = c(0, 1, 2)) +
        ylab("Survival probability") +
        #ggtitle("Lambs") +
        xlab(expression(F[ROH])) 
p_lamb
p_adult <-  ggplot(roh_eff_adult, aes(x, predicted)) +
        geom_line(aes(color = group), size = 1) +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group, color = group),
                    alpha = 0.05, size = 0.2, linetype = 2) +
        theme_simple() +
        scale_color_viridis_d("age", labels = c(1, 4, 8)) +
        scale_fill_viridis_d("age", labels = c(1, 4, 8)) +
        scale_y_continuous(labels = paste0(seq(from = 0, to = 100, by = 25), "%"), limits = c(0,1)) +
        scale_x_continuous(labels = c("0.24", "0.34", "0.44"), breaks = c(0, 1, 2)) +
        ylab("Survival probability") +
        xlab(expression(F[ROH])) 
p_adult  

p_final <- p_froh_across_ages + p_surv_mod + p_lamb + p_adult + plot_annotation(tag_levels = "A") +
        plot_layout(nrow = 2, heights = c(1.4,1))
ggsave("figs/model_plot.jpg", p_final, width = 10, height = 8)
