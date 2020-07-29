# additional figures of interest for supplementary material

library(tidyverse)
library(optiSel)
library(ggExtra)
library(data.table)
source("theme_simple.R")
library(snpStats)
#load("data/fitness_roh_df.RData") # formerly
load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")
library(broom)

# genotype data
sheep_plink_name <- "data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# roh data
file_path <- "output/ROH/roh.hom"
roh_lengths <- fread(file_path) 

ped <- sheep_ped[, c(1,3,2)]
ped_fin <- prePed(ped)
ID <- pedInbreeding(ped_fin) %>% as_tibble() %>% rename(id = Indiv, fped = Inbr)
ggplot(ID, aes(fped)) +
        geom_histogram(bins = 500) +
        scale_y_log10() 
        
# survival data 
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

# correlation FROH, FPED--------------------------------------------------------
F_df <- annual_survival %>% 
        left_join(ID) %>% 
        dplyr::select(id, sex, fped, froh_all, froh_long, froh_medium, froh_short)

lm(froh_all ~ fped, data = F_df) %>% summary()
df_rel <- data.frame(rel = fct_inorder(c("First cousins once removed", "First cousins", "Half-siblings", "Parents-offspring")), 
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

# annual survival plot----------------------------------------------------------
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



# additional figures from review -----------------------------------------------

# MAF distribution
snp_stats <- col.summary(full_sample$genotypes)

p_maf <- ggplot(snp_stats, aes(MAF)) + 
        geom_histogram(binwidth = 0.005,  color = "white", size = 0.1, position = "identity",
                       alpha = 1, fill = "#4c566a") +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
        xlab("minor allele frequency")
        
ggsave("figs/Sup_MAF.jpg", width = 4, height = 1.8)

library(gghalves)
# FROH change over time.
annual_survival$birth_year <- as.numeric(as.character(annual_survival$birth_year))
annual_survival$birth_year <- as.factor((annual_survival$birth_year))

mod <- lm(froh_all ~ as.numeric(birth_year), data = annual_survival)
summary(mod)
tidy(mod, conf.int = TRUE)
coefs <- coef(lm(froh_all ~ as.numeric(birth_year), data = annual_survival))
coefs

p_froh_time <- ggplot(annual_survival, aes(birth_year, froh_all)) +
       # geom_point(shape = 21, size =0.5, stroke = 0.05) +
        geom_half_point(side = "l", shape = 21, alpha = 1, stroke = 0.1, size = 1,
                        transformation_params = list(height = 0, width = 0.5, seed = 1),
                        fill = "#d8dee9") +
        geom_half_boxplot(side = "r", outlier.color = NA,
                          width = 0.6, lwd = 0.3, color = "black",
                          alpha = 0.8, fill = "#d8dee9") +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE) +
       # scale_fill_manual(values = c("#4c566a", "#d8dee9")) +
        #scale_color_manual(values = c( "#d8dee9", "#4c566a")) +
        geom_abline(intercept = coefs[1], slope = coefs[2], color = "black", size = 0.3) +
        xlab("Birth Year") +
        ylab(expression(F[ROH])) +
        theme(axis.text.x = element_text(angle = 45)) 

ggsave("figs/sup_froh_over_time.jpg", p_froh_time, width = 8, height = 2.5)



