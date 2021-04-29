# Plotting survival model and raw data (Figure 3)

library(ghibli)
library(tidyverse)
library(data.table)
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
library(inlafuns)
library(ggridges)
library(gt)
library(Hmisc)
# data
# annual measures of traits and fitness
load("data/survival_mods_data.RData")
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
        as.data.frame()  %>% 
        mutate(life_stage = case_when(
                age == 0 ~ "lamb",
                age > 0 & age <= 2 ~ "early_life",
                age > 2 & age <= 4 ~ "mid_life",
                age >= 5 ~ "late_life",
        )) %>% 
        mutate(life_stage = factor(life_stage, levels = c( "lamb","early_life", "mid_life", "late_life")))

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

#saveRDS(ids_per_age, file = "output/ids_per_age.rds")
# load ROH info
file_path <- "output/ROH/roh.hom"
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
        mutate(FROH = MB_sum/2452.07) %>% ### check this / checked
        mutate(age = str_replace(age, "_", " ")) %>% 
        mutate(age_num = as.numeric(str_replace(age, "age ", ""))) %>% 
        mutate(age_num_fct = factor(age_num, levels = as.character(0:10))) %>% 
        as.data.frame() %>% 
        mutate(age = factor(age)) %>% 
        mutate(age = fct_reorder(age, age_num))-> roh_plot2

roh_plot2 %>% 
        filter(age_num %in% c(0:9)) %>% 
        #ungroup() %>% 
        # sample_frac(0.1) %>% 
        ggplot(aes(x = FROH, y = age_num_fct)) +
        geom_density_ridges(scale = 0.83,
                            jittered_points=TRUE, color = "#4c566a",  # # "#4c566a"  "#eceff4"
                            fill = "#eceff4",
                            point_shape = "|", point_size = 2.5, size = 0.4,
                            position = position_points_jitter(height = 0),
                            point_alpha = 1, point_color = "#4c566a", 
                            alpha = 0.9, bandwidth = 0.006) +
        theme_simple(axis_lines = FALSE, base_size = 14, base_family = "Helvetica") +
        #scale_x_continuous(limits = c(0.199, 0.55)) +
        ylab("Age") + 
        scale_y_discrete(expand = expand_scale(add = c(0.2, 1.01))) +
        xlab(expression(italic(F)[ROH])) + # ~across~age~cohorts
        theme(strip.background = element_blank(),
              panel.grid.major = element_line(colour = "#d8dee9", size = 0.5),
              strip.text = element_text(vjust=2),
              panel.spacing.x = unit(1, "lines"),
              #plot.margin=unit(c(1,0,0,0),"cm"),
              axis.text = element_text(size = 13, color = "black"),
              strip.text.x = element_text(margin = margin(1,0,0,0, "cm")),
              legend.position = "none") -> p_froh_across_ages #+
#labs(title = "Inbreeding across life in Soay sheep",
#    subtitle = "- individuals with lots of long ROHs rarely survive their first year")
#p_froh_across_ages
# ggsave("figs/FROH_across_ages_col.jpg", plot = p_froh_across_ages, width = 4.5, height = 5.5)


# Plot B: effect sizes ---------------------------------------------------------
mod_inla <- readRDS("output/AS_mod_oar2_life_stage.rds")

set.seed(145)
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
        .[c(2,8:10), ] %>% 
        mutate(pred = factor(pred, levels = rev(pred)))
names(fix_eff) <- c("Predictor", "mean", "lower_CI", "upper_CI")

# decrease in odds of survival with 10% increase in FROH
1-fix_eff

p_surv_mod <- ggplot(fix_eff, aes(mean, Predictor, xmax = upper_CI, xmin = lower_CI)) +
        geom_vline(xintercept = 1, linetype='dashed', colour =  "#4c566a", size = 0.3) +     # "#4c566a"  "#eceff4"
        geom_errorbarh(alpha = 1, height = 0,
                       size = 0.5) +
        geom_point(size = 3, shape = 21, col = "#4c566a", fill = "#eceff4", # "grey69"
                   alpha = 1, stroke = 0.7) + 
        # when only showing roh related effs
        scale_x_log10(breaks = c(0.1, 0.3, 1, 2, 5), limits = c(0.1, 5), labels = c( "0.1", "0.3", "1", "2", "5")) +
        scale_y_discrete(labels = rev(c(expression(F[ROH]), 
                                        expression(F[ROH]~'*'~Early~life), 
                                        expression(F[ROH]~'*'~Mid~life),
                                        expression(F[ROH]~'*'~Late~life))),
                         ) +
        theme_simple(axis_lines = TRUE, base_size = 14, base_family = "Helvetica") +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.text = element_text(size = 12),
                axis.title.x = element_text(margin=margin(t=8))
        ) +
        xlab("Odds-ratio and 95% CI")-> p_forest
#xlab(expression(beta~and~95*"%"*~CI~(odds~of~survival))) -> p_forest
p_forest

# Plot C marginal effects ------------------------------------------------------
#mod_inla <- readRDS("output/AS_mod_INLA_400k.rds")

# plotting marginal effects (predictions) from INLA models is not 
# entirely straightforward. This solution is based on an answer by H. Rue from
# the INLA google group.

# plot INLA marginal effects
fun <- function(...) {
        # plogis = inverse logit
        one <-  plogis(Intercept + 
                               df1$x1 * froh_all10_cent + 
                               df1$x2 * life_stageearly_life + 
                               df1$x3 * life_stagemid_life + 
                               df1$x4 * life_stagelate_life +
                               df1$x5 * sexM +  
                               df1$x6 * twin1 +
                               df1$x7 * `froh_all10_cent:life_stageearly_life` + 
                               df1$x8 * `froh_all10_cent:life_stagemid_life` +
                               df1$x9 * `froh_all10_cent:life_stagelate_life`)
        
        return (list(one))
}

# values here look odd, because variables have been centered around the mean
# for easier modeling
froh <- seq(from = min(annual_survival$froh_all10_cent), to = (max(annual_survival$froh_all10_cent)), by = 0.1)
age <- c("life_stagelamb","life_stageearly_life", "life_stagemid_life", "life_stagelate_life")
#age <- c(-2.4, -1.4, 0.6, 3.6, 6.6)
#age <- c(-2.4, -1.4, -0.4, 0.6, 1.6, 2.6, 3.6, 4.6, 5.6, 6.6)
combined_df <- expand_grid(froh, age) %>% 
        mutate(
               # lamb = ifelse(age == "life_stagelamb", 1, 0),
               early = ifelse(age == "life_stageearly_life", 1, 0),
               mid = ifelse(age == "life_stagemid_life", 1, 0),
               late = ifelse(age == "life_stagelate_life", 1, 0),
               sex = 1,
               twin = 0,
               #twin =  0.5,
               #sex = 0.5,
               frohxearly = froh*early,
               frohxmid = froh*mid,
               frohxlate = froh*late) %>% 
        select(-age)

names(combined_df) <- paste0("x", 1:9)

set.seed(144)
xx <- inla.posterior.sample(1000, mod_inla)

marg_means <- purrr::map(1:nrow(combined_df), function(x) {
        df1 <<- combined_df[x, ]
        out <- inla.posterior.sample.eval(fun, xx)
}) 
rm(df1)

d <- marg_means %>% 
        map(as_tibble) %>% 
        bind_rows() %>% 
        pmap_df(function(...) {
                samp <- as.numeric(unlist(list(...)))
                c(mean = mean(samp), quantile(samp, probs = c(0.025, 0.975)))
        }) %>% 
        bind_cols(combined_df) %>% 
        mutate(life_stage = expand_grid(froh, age)$age) %>% 
        mutate(life_stage = case_when(
                life_stage == "life_stagelamb" ~ "Lamb",
                life_stage =="life_stageearly_life" ~ "Early life",
                life_stage == "life_stagemid_life" ~ "Mid life",
                life_stage == "life_stagelate_life" ~ "Late life"
        )) %>% 
        .[, c(1:4, 13)] %>% 
        setNames(c("prediction", "ci_lower", "ci_upper", "froh", "age")) %>% 
        mutate(froh = (froh + mean(annual_survival$froh_all10))/10)


#saveRDS(d, file = "output/AS_mod_INLA_398k_predictions_for_plot2.rds")
#d <- readRDS("output/AS_mod_INLA_398k_predictions_for_plot2.rds")
inla_preds <- d %>% 
        mutate(prediction = prediction * 100,
               ci_lower = ci_lower * 100,
               ci_upper = ci_upper * 100)

inla_preds <- inla_preds %>% 
        mutate(age = factor(age, levels = c("Lamb", "Early life", "Mid life", "Late life")))

p_marginal_effs <- ggplot(inla_preds, aes(froh, prediction)) +
        geom_ribbon(aes(x=froh, ymin = ci_lower, ymax = ci_upper, fill = age, color = age),
                    alpha = 0.03, linetype = 2, size = 0.2) +
        geom_line(aes(color = age), size = 0.9, alpha = 1) +
        scale_color_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        scale_fill_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 14,
                     base_family = "Helvetica") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
        theme(axis.line.y = element_blank(),
              legend.position = "none",
              #axis.text.x = element_text(size = 11),
              axis.text = element_text(color = "black", size = 12)) +
        xlab(expression(italic(F)[ROH])) +
        ylab("Predicted\nsurvival probability\n(% per year)")
p_marginal_effs


# plot raw data
surv_per_F <- function(iter) {
        out <- annual_survival %>%
                dplyr::mutate(binned_froh = cut2(froh_all10, cuts = c(1.8, 2, 2.2, 2.4, 2.6, 2.8, 5.1))) %>% 
                dplyr::group_by(life_stage, binned_froh) %>% 
               # filter(n() >= 20) %>% 
                dplyr::sample_n(size = nrow(annual_survival), replace = TRUE) %>% 
                dplyr::select(froh_all10, binned_froh, survival) %>% 
                dplyr::summarise(binned_survival = mean(survival)) 
        out
}

all_boot <- map_df(1:10, surv_per_F)

emp <- annual_survival %>% 
        dplyr::mutate(binned_froh = cut2(froh_all10, cuts = c(1.8, 2, 2.2, 2.4, 2.6, 2.8, 5.1))) %>% 
        dplyr::group_by(life_stage, binned_froh) %>% 
        dplyr::select(froh_all10, binned_froh, survival) %>% 
        dplyr::summarise(binned_survival = mean(survival)) 
emp
# make source data for nat comms
write_xlsx(emp, path = "output/source_data_figure3b.xlsx")

p_raw <- ggplot(emp , aes(binned_froh, binned_survival, fill = life_stage)) +
        #geom_jitter(width = 0.05, size = 2, alpha = 0.3, shape = 21, stroke=0.2) +
        geom_point(shape = 21, size = 3, stroke = 0.4, color = "black")  +
        scale_fill_viridis_d("Life Stage (age)", 
                             labels = c("Lamb (0)", "Early life (1,2)",
                                        "Mid life (3,4)",
                                        "Late life (5+)"), option = "D") +
        scale_y_continuous(breaks = seq(0.25, 1, 0.25), labels = seq(25, 100, 25),
                           limits = c(0.25, 1)) +
        scale_x_discrete(labels = c("[0.18,0.20)", "[0.20,0.22)", "[0.22,0.24)", 
                                    "[0.24,0.26)", "[0.26,0.28)", "[0.28,0.51)"),
                         guide = guide_axis(n.dodge = 2)) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE,
                     base_family = "Helvetica") +
        theme(axis.line.y = element_blank(),
              axis.text.y = element_text(color = "black"),
              axis.text.x = element_text(color = "black"))+
        xlab(expression(italic(F)[ROH]~class)) +
        ylab("Survivors (% per year)")
#geom_line(mapping = aes(group =age), size = 0.2, alpha = 1) 
p_raw


p_final <- p_froh_across_ages + (p_raw / p_marginal_effs) +
        plot_layout(guides = "collect", widths = c(0.57, 0.43)) &
        plot_annotation(tag_levels = "a") &
        theme(axis.text = element_text(size = 11),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 11),
              legend.title = element_text(size = 12),
              plot.tag = element_text(face = "bold"))

p_final <- p_froh_across_ages + (p_raw / p_marginal_effs) +
  plot_layout(guides = "collect", widths = c(0.57, 0.43)) &
  plot_annotation(tag_levels = "a") &
  theme(axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        plot.tag = element_text(face = "bold"))
p_final

ggsave("figs/Fig3_inla2.jpg", plot = p_final, height = 6.5, width = 9.5, dpi = 900)
ggsave("figs/Fig3_inla2.pdf", plot = p_final, height = 6.5, width = 9.5)
#ggsave("figs/Fig2_inla2_2.jpg", plot = p_final, height = 8.5, width = 12.5)

# library(cowplot)
# p_final <- cowplot::plot_grid(p_froh_across_ages, 
#                               cowplot::plot_grid(p_raw, p_marginal_effs, label_size = 16, nrow = 2, rel_heights = c(1,1.8), labels = c('B', 'C')),
#                               labels = c('A', ''), rel_widths = c(1, 0.9), label_size = 16)
# ggsave("figs/Fig2_inla2.jpg", height = 6, width = 8.5)

# predictions on original scale
inla_preds <- d

p_marginal_effs_logit <- ggplot(inla_preds, aes(froh, prediction)) +
        geom_line(aes(color = age), size = 1) +
        geom_ribbon(aes(x=froh, ymin = ci_lower, ymax = ci_upper, fill = age, color = age),
                    alpha = 0.1, linetype = 2, size = 0.2) +
        scale_color_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        scale_fill_viridis_d("Age", labels = c(0, 1, 4, 7)) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 14) +
        theme(axis.line.y = element_blank(),
              legend.position = "top",
              axis.text = element_text(size = 13)) +
        xlab(expression(F[ROH])) +
        ylab("log-odds(survival)")
ggsave("figs/Fig2_survival_logodd.jpg", height = 3.5, width = 3.5)

#mod_inla <- readRDS("output/AS_mod_INLA_397k.rds")



# make Supplementary table / figure for model

fix_eff <- mod_inla$summary.fixed %>% 
        map_df(round, 2) %>% 
        select(-kld) %>% 
        setNames(c("mean", "std_err", "ci_lower", "median", "ci_upper", "mode")) %>% 
        mutate(term = rownames(mod_inla$summary.fixed)) %>% 
        select(term, everything()) %>% 
        mutate(mean = paste0(mean,' (', round(exp(mean),2), ')'),
               std_err = paste0(std_err, ' (', round(exp(std_err),2), ')'),
               ci_lower = paste0(ci_lower, ' (', round(exp(ci_lower),2), ')'),
               median = paste0(median, ' (', round(exp(median),2), ')'),
               ci_upper = paste0(ci_upper, ' (', round(exp(ci_upper),2), ')'),
               mode = paste0(mode, ' (', round(exp(mode),2), ')'))
        

fix_eff2 <- fix_eff %>% 
        mutate(add_info = c("", "continuous", "categorical (0=no, 1=yes)",
                            "categorical (0=no, 1=yes)",
                            "categorical (0=no, 1=yes)",
                            "categorical (0=female, 1=male)", 
                            "categorical (0=no, 1=yes)",
                            "","","")) %>% 
        mutate(standardisation = c("", "(x * 10)-mean(x * 10)", rep("", 8))) %>% 
        mutate(term = c("Intercept", "F<sub>ROH</sub>", 
                        "LifeStage: EarlyLife",
                        "LifeStage: MidLife",
                        "LifeStage: LateLife",
                        "Sex", "Twin",
                        "F<sub>ROH</sub> * (LifeStage: EarlyLife)", 
                        "F<sub>ROH</sub> * (LifeStage: MidLife)",
                        "F<sub>ROH</sub> * (LifeStage: LateLife)"
                        ))
library(inlafuns)
raneff <- get_raneff(mod_inla, scales = "sd") %>% 
        mutate(across(is.numeric, round, 2)) %>% 
        mutate(term = c("Birth year", "Capture year", "Individual", "Add. genetic")) %>% 
        mutate(add_info = c("n = 40", "n = 40", "n = 5952", "Pedigree-based")) %>% 
        mutate(standardisation = "") %>% 
        mutate_if(is.numeric, as.character)
#mutate(groups = c(44, 44, 5952, ""))

mod_tab <- bind_rows(fix_eff2, raneff) %>% 
        select(term, mean, std_err, ci_lower, ci_upper, add_info, standardisation) %>% 
        mutate(effect = c(rep("Fixed effects", 10), 
                          rep("Random effects (variances)", 4))) %>% 
        select(8, 1:5, 7, 6) %>% 
        setNames(c("effect", "Term", "Post.Mean", "Std.Error", "CI (2.5%)", "CI (97.5%)", "Standardisation", "Info"))

mod_tab %>% gt(
        rowname_col = "term",
        groupname_col = "effect"
) %>% 
        tab_style(
                style = cell_text( weight = "bold"),
                locations = cells_column_labels(columns = TRUE)
        ) %>% 
        fmt_markdown(columns = TRUE) %>% 
        gtsave("AS_model_table.png", path = "figs/tables/")

# differences in inbreeding depression across life
mod_inla2 <- readRDS("output/AS_mod_oar2_life_stage_early_ref.rds")
fix_eff2 <- mod_inla2$summary.fixed %>% 
  map_df(round, 2) %>% 
  select(-kld) %>% 
  setNames(c("mean", "std_err", "ci_lower", "median", "ci_upper", "mode")) %>% 
  mutate(term = rownames(mod_inla2$summary.fixed)) %>% 
  select(term, everything()) %>% 
  mutate_if(is.numeric, exp)


# supplementary figure like 2C but with all age classes and log-odds
age <- c(-2.4, -1.4, -0.4, 0.6, 1.6, 2.6, 3.6, 4.6, 5.6, 6.6)

combined_df <- expand_grid(froh, age) %>% 
        mutate(lamb = ifelse(age == -2.4, 1, 0),
               twin = 0,
               sex = 1,
               #twin =  0.5,
               #sex = 0.5,
               frohxlamb = froh*lamb,
               frohxage = froh*age) 

names(combined_df) <- paste0("x", 1:7)

set.seed(144)
xx <- inla.posterior.sample(1000, mod_inla)

marg_means <- purrr::map(1:nrow(combined_df), function(x) {
        df1 <<- combined_df[x, ]
        out <- inla.posterior.sample.eval(fun, xx)
}) 
marg_means_link <- purrr::map(1:nrow(combined_df), function(x) {
        df1 <<- combined_df[x, ]
        out <- inla.posterior.sample.eval(fun_nolink, xx)
})
rm(df1)

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

d_link <- marg_means_link %>% 
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

inla_preds <- d %>% 
        mutate(prediction = prediction * 100,
               ci_lower = ci_lower * 100,
               ci_upper = ci_upper * 100)
inla_preds_link <- d_link

p_marginal_effs <- ggplot(inla_preds, aes(froh, prediction)) +
        geom_line(aes(color = age), size = 0.5) +
        geom_ribbon(aes(x=froh, ymin = ci_lower, ymax = ci_upper, fill = age, color = age),
                    alpha = 0.01, linetype = 2, size = 0.2) +
        scale_color_viridis_d("Age", labels = 0:9) +
        scale_fill_viridis_d("Age", labels = 0:9) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 14) +
        theme(axis.line.y = element_blank(),
              legend.position = "top",
              axis.text = element_text(size = 13)) +
        xlab(expression(F[ROH])) +
        ylab("Predicted\nsurvival probability %") +
        theme(legend.position = "top")

p_marginal_effs_link <- ggplot(inla_preds_link, aes(froh, prediction)) +
        geom_line(aes(color = age), size = 0.5) +
        geom_ribbon(aes(x=froh, ymin = ci_lower, ymax = ci_upper, fill = age, color = age),
                    alpha = 0.01, linetype = 2, size = 0.2) +
        scale_color_viridis_d("Age", labels = 0:9) +
        scale_fill_viridis_d("Age", labels = 0:9) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 14) +
        theme(axis.line.y = element_blank(),
              legend.position = "top",
              axis.text = element_text(size = 13)) +
        xlab(expression(F[ROH])) +
        ylab("Predicted\n survival (log-odds scale)") +
        theme(legend.position = "top")

p_effs <- p_marginal_effs + p_marginal_effs_link + 
        plot_annotation(tag_levels = 'A') +
        plot_layout(guides = 'collect') & theme(legend.position = 'right')

ggsave("figs/Sup_model_pred_surv.jpg", p_effs, height = 4, width = 9)


# lme4
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

# age interaction and lamb
library(lme4)
library(broom.mixed)
mod_lme4 <- glmer(survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
                  family = binomial, data = annual_survival,
                  control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
tidy(mod_lme4, effects = "ran_pars", scales = "vcov")

mod_lme42 <- lmer(survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id),
                  data = annual_survival)

tidy(mod_lme42, effects = "ran_pars", scales = "vcov")
