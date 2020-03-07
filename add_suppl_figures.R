# additional figures of interest for supplementary material

library(tidyverse)
library(optiSel)
library(ggExtra)
source("theme_simple.R")
#load("data/fitness_roh_df.RData") # formerly
load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")

# roh data
file_path <- "data/roh_nofilt_ram_pruned.hom"
roh_lengths <- fread(file_path) 

ped <- sheep_ped[, c(1,3,2)]
ped_fin <- prePed(ped)
ID <- pedInbreeding(ped_fin) %>% as_tibble() %>% rename(id = Indiv, fped = Inbr)
ggplot(ID, aes(Fped)) +
        geom_histogram(bins = 500) +
        scale_y_log10() 
        

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

# correlation FROH, FPED
F_df <- annual_survival %>% 
        left_join(ID) %>% 
        dplyr::select(id, sex, fped, froh_all, froh_long, froh_medium, froh_short)

lm(froh_all ~ fped, data = F_df) %>% summary()
df_rel <- data.frame(rel = fct_inorder(c("First cousins once removed", "First cousins", "Half-siblings", "Parents-Children")), 
                     fped = c(0.03125, 0.0625, 0.125, 0.25))
p_F <- ggplot(F_df, aes(fped, froh_all)) +
        geom_point(size = 1.3, alpha = 1, shape = 21, fill = "#d8dee9", color = "#2e3440", stroke = 0.1) +        
        geom_vline(data = df_rel, aes(xintercept = fped, color = rel), 
                   linetype = 2, size = 0.3) +
        labs(color = "Parent relatedness") +
        scale_color_viridis_d() +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 14) +
        ylab(expression(F[ROH])) +
        xlab(expression(F[PED]))

ggsave("figs/sup_Froh_vs_Fped.jpg", height = 3, width = 7)
# num individuals
annual_survival %>% group_by(id) %>% tally()
annual_survival %>% group_by(id) %>% summarise(age = max(age)) %>% summarise(mean_age = mean(age))
# num observations
nrow(annual_survival)


# annual survival plot
p_age_df <- annual_survival %>% group_by(age, sex) %>% tally()
p_age <- ggplot(p_age_df , aes(age,n,fill=sex)) +
        geom_col(position = "dodge",color = "#2e3440", size = 0.2) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12, 
                     base_family = "Lato") +
        scale_fill_manual("Sex", values = c("#4c566a", "#e5e9F0")) +
        ylab("Individuals") +
        xlab("Age") +
        scale_y_continuous(expand = c(0,0))

ggsave("figs/sup_age_dist.jpg", height = 2.5, width = 4)
