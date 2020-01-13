library(tidyverse)
source("theme_clean.R")
library(boot)
library(wesanderson)
# Part 1: overall ==============================================================
# results
out <- readRDS("output/inla_survival_models.rds")

inla_mod1 <- out[[1]]
inla_mod2 <- out[[2]]

summary(inla_mod1)
summary(inla_mod2)

sampvars <- 1/inla.hyperpar.sample(1000, inla_mod1)
sampicc <- sampvars[,3]/ (rowSums(sampvars))
quantile(sampicc, c(0.025, 0.5, 0.975))

# inv logit exp(x)/(1+exp(x))
inv_logit <- function(x) exp(x)/(1+exp(x))

# plot
pred_names <- c("FROH", "FROH long", "FROH medium", "FROH short", "sex (male)", "twin", "age", "age2")
surv_mod_df <- out %>% 
                map("summary.fixed") %>% 
                map(as_tibble, rownames = "predictor") %>% 
                bind_rows(.id = "mod") %>%
                rename(lower_ci = `0.025quant`,
                        upper_ci = `0.975quant`)

# formatC(inv_logit(seq(0, -50, by = -10)), format = "e", digits = 1)
# mean FROH_all
annual_survival %>% group_by(id) %>% summarise(froh = mean(froh_all)) %>% summarise(mean(froh))

# calculate FROH as odds of surviving when offspring of first cousins vs average
surv_mod_cousin <- surv_mod_df %>% 
        filter(str_detect(predictor, "froh")) %>% 
        mutate_at(c("mean", "lower_ci", "upper_ci"), function(x) exp((0.226+0.0625) * x) / exp(0.226 * x)) 

surv_mod_plot <- surv_mod_df %>% 
        filter(!str_detect(predictor, "froh")) %>% 
        mutate_at(c("mean", "lower_ci", "upper_ci"), exp) %>% 
        #filter(!str_detect(predictor, "age")) %>% 
        rbind(surv_mod_cousin, .) 

surv_mod_plot %>% 
        filter(mod == 1 & predictor %in% c("froh_all", "sexM", "twin")) %>% 
        ggplot(aes(mean, predictor)) +
        geom_errorbarh(aes(xmax = upper_ci, xmin = lower_ci), height = 0, size = 0.5, color = "black") +
        geom_point(size = 3, shape=21,fill = "grey", color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
        scale_y_discrete(labels = rev(c("FROH\n(Offspring of cousins)", "Sex\n(Male)", "Twin"))) +
        #scale_y_discrete(labels = rev(c(expression(atop(F[ROH],(Offspring~of~cousins))), "Sex\n(Male)", "Twin"))) +
        scale_x_continuous(limits=c(0,1)) +
        xlab("odds of survival\n(mean & credible interval)") +
        theme_clean() +
        theme(axis.line.y = element_blank(),
              axis.ticks.x = element_line()) +
        ylab("") -> p_froh

p_froh
ggsave(filename = "figs/surv_mod_froh.jpg", p_froh, width = 4, height = 2.3)        


# FROH in different length classes
surv_mod_plot2 <- surv_mod_plot%>% 
        filter(str_detect(predictor, "froh")) %>% 
        mutate(predictor = as.factor(predictor))
        
ggplot(surv_mod_plot2, aes(mean, predictor)) +
        geom_errorbarh(aes(xmax = upper_ci, xmin = lower_ci), height = 0, size = 0.5, color = "black") +
        geom_point(size = 3, shape=21,fill = "grey", color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
        scale_y_discrete(limits =  rev(levels(surv_mod_plot2$predictor)),
                         labels = rev(c("FROH", "FROH long\n(> 5MB)", 
                                        "FROH medium\n(1MB - 4MB)", 
                                        "FROH short\n(< 1MB)"))) +
        #scale_x_continuous(limits=c(0,1)) +
        xlab("odds of survival\n(mean & credible interval)") +
        theme_clean() +
        theme(axis.line.y = element_blank(),
              axis.ticks.x = element_line()) +
        ylab("") -> p_froh2
ggsave(filename = "figs/surv_mod_froh_by_length.jpg", p_froh2, width = 4, height = 2.6)     

# rpt and herit
sampvars <- 1/inla.hyperpar.sample(1000, inla_mod1)
sampicc <- sampvars[, 4]/ (rowSums(sampvars) + (pi^2)/3)
quantile(sampicc, c(0.025, 0.5, 0.975))



# Part 2: FROH in different age classes ========================================

all_inla_mods <- readRDS("output/old_models/inla_survival_models_diff_ages_froh_all_new.rds")
all_inla_mods_froh <- readRDS("output/old_models/inla_survival_models_diff_ages_new.rds")
inla_mods <- list(all_inla_mods, all_inla_mods_froh)
# extract all fixed effects
fix_effs <- inla_mods %>% 
        flatten() %>% 
        map("summary.fixed") %>% 
        #compact() %>% 
        map(rownames_to_column, var = "var") %>% 
        bind_rows(.id = "age") %>% 
        mutate(age = as.numeric(age) - 1) %>% 
        mutate(age = ifelse(age > 8, age - 9, age)) %>% 
        #filter(age != 9) %>% 
        #filter(var %in% c("froh_all_std", "froh_long_std", "froh_medium_std", "froh_short_std"))
        filter(var %in% c("froh_all", "froh_long", "froh_medium", "froh_short")) %>% 
        mutate(var = case_when(
                var =="froh_all" ~ "FROH",
                var == "froh_long" ~ "FROH long",
                var == "froh_medium" ~ "FROH medium",
                var == "froh_short" ~ "FROH short"
        ))

#filter(str_detect(var, "froh"))
names(fix_effs) <- c("age", "var", "mean", "sd", "lower", "median", "upper", "mode", "kld")
# transform
fix_effs_trans <- fix_effs %>% 
        mutate_at(c("mean", "lower", "upper"), function(x) exp((0.226+0.0625) * x) / exp(0.226 * x)) 
fix_effs_trans <- fix_effs %>% 
        mutate_at(c("mean", "lower", "upper"), function(x) 0.0625 * x )

# only FROH
ggplot(fix_effs_trans %>% filter(var == "FROH"), aes(age, mean)) + 
        geom_hline(aes(yintercept = 0), color = "lightgrey") +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, color = "#3B4252") +
        geom_point(size = 2.5, shape = 21, fill = "#046C9A") +
        #geom_smooth(method="lm", se = FALSE) +
        # facet_wrap(~var, scales = "free_y", ncol = 2) +
        theme_clean() +
        ylab("FROH estimate\n(log odds)") + 
        # scale_color_manual(values = wes_palette("Darjeeling2")) +
        theme(legend.position = "none") -> p_roh_froh
ggsave("figs/FROH_across_ages_full.jpg", width = 3.3, height = 2.4)

ggplot(fix_effs, aes(age, mean, color = var)) + 
        geom_hline(aes(yintercept = 0), color = "lightgrey") +
        geom_point() +
        #geom_smooth(method="lm", se = FALSE) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
        facet_wrap(~var, scales = "free_y", ncol = 2) +
        theme_clean() +
        ylab("inbreeding depression estimate") + 
        scale_color_manual(values = wes_palette("Darjeeling2")) +
        theme(legend.position = "none") -> p_roh
p_roh
ggsave("figs/ID_across_life.jpg", p_roh, width = 4, height = 3)

annual_survival %>% group_by(age, sex) %>% tally() %>% 
        mutate(sex = as.character(sex)) %>% 
        filter(!is.na(sex)) %>% 
        filter(age %in% c(1:8)) %>% 
        ggplot(aes(age, n, color = sex)) + 
        geom_point() +
        geom_line() +
        theme_clean()

# inbreeding depression across life ============================================
inla_mods <- readRDS("output/inla_survival_models_interaction_excl_age0.rds")
surv_mod_df <- inla_mods %>% 
        map("summary.fixed") %>% 
        map(as_tibble, rownames = "predictor") %>% 
        bind_rows(.id = "mod") %>%
        rename(lower_ci = `0.025quant`,
               upper_ci = `0.975quant`)

# calculate FROH as odds of surviving when offspring of first cousins vs average
# surv_mod_cousin <- surv_mod_df %>% 
#         filter(str_detect(predictor, "froh")) %>% 
#         mutate_at(c("mean", "lower_ci", "upper_ci"), function(x) exp(0.0625 * x) / exp(0)) 

surv_mod_plot <- surv_mod_df %>% 
        filter(mod == 1 & predictor %in% c("froh_all_cent", "froh_all_cent:age_cent", "sexM", "twin")) %>% 
        mutate(mean = ifelse(predictor == "froh_all_cent", 0.0625 * mean, mean),
               lower_ci = ifelse(predictor == "froh_all_cent", 0.0625 * lower_ci, lower_ci),
               upper_ci = ifelse(predictor == "froh_all_cent", 0.0625 * upper_ci, upper_ci)) %>% 
        mutate(predictor = as.factor(predictor)) 


ggplot(surv_mod_plot, aes(mean, predictor)) +
        geom_errorbarh(aes(xmax = upper_ci, xmin = lower_ci), height = 0, size = 0.5, color = "black") +
        geom_point(size = 2.5, shape=21,fill = "grey", color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
        scale_y_discrete(limits = rev(levels(surv_mod_plot$predictor)),
                         labels = rev(c("FROH\n(Offspring \nof cousins)", "FROH * age", "Sex\n(Male)", "Twin"))) +
        xlab("log(odds) of survival") +
        theme_clean() +
        theme(axis.line.y = element_blank(),
              axis.ticks.x = element_line()) +
        ylab("") -> p_froh

p_froh
ggsave(filename = "figs/surv_mod_froh_interaction_incl_age0.jpg", p_froh, width = 3.3, height = 2.4)        


# plot per froh length class
surv_mod_plot <- surv_mod_df %>% 
        filter(mod == 2 & predictor %in% c("froh_long_cent",
                                           "froh_medium_cent",
                                           "froh_short_cent",
                                           "froh_long_cent:age_cent",
                                           "age_cent:froh_medium_cent",
                                           "age_cent:froh_short_cent")) %>% 
        mutate(mean = ifelse(predictor %in% c("froh_long_cent", "froh_medium_cent", "froh_short_cent"), 0.0625 * mean, mean),
               lower_ci = ifelse(predictor %in% c("froh_long_cent", "froh_medium_cent", "froh_short_cent"), 0.0625 * lower_ci, lower_ci),
               upper_ci = ifelse(predictor %in% c("froh_long_cent", "froh_medium_cent", "froh_short_cent"), 0.0625 * upper_ci, upper_ci)) %>% 
        mutate(predictor = as.factor(predictor)) %>% 
        mutate(predictor = fct_relevel(predictor, "froh_long_cent", "froh_medium_cent",
                                       "froh_short_cent", "froh_long_cent:age_cent",
                                       "age_cent:froh_medium_cent", "age_cent:froh_short_cent"))


ggplot(surv_mod_plot, aes(mean, predictor)) + 
        geom_errorbarh(aes(xmax = upper_ci, xmin = lower_ci), height = 0, size = 0.5, color = "black") +
        geom_point(size = 2, shape=21,fill = "grey", color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
        scale_y_discrete(limits = rev(levels(surv_mod_plot$predictor)),
                         labels = rev(c("FROH long", "FROH medium", "FROH short",
                                        "FROH long * age", "FROH medium * age", "FROH short * age"))) +
        xlab("log(odds) of survival") +
        theme_clean() +
        theme(axis.line.y = element_blank(),
              axis.ticks.x = element_line()) +
        ylab("") -> p_froh

ggsave(filename = "figs/surv_mod_froh_interaction_excl_age0_roh_classes_wide.jpg", p_froh, width = 4.7, height = 2.4) 
