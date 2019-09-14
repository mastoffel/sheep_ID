# plot 
library(tidyverse)
library(brinla)
source("theme_clean.R")
mod_lambs <- readRDS("output/models/mods_lambs.RDS")

all_fixed <- map(mod_lambs$fitted, function(x) {
                        out <- tibble::rownames_to_column(x$summary.fixed, var = "var")
                        out %>% 
                                filter(str_detect(var, "froh"))
            }) %>% 
            .[c(1,2,4,5)] %>% 
            bind_rows() %>% 
            mutate(response = rep(c("hindleg", "weight"), each = 4))

all_heritabilities <- map(mod_lambs$fitted, function(x) {
        
                        rand_vars <- 1/inla.hyperpar.sample(1000, x)
                        add_gen <- which(str_detect(colnames(rand_vars), "IndexA"))
                        herit <- rand_vars[ ,add_gen]/ (rowSums(rand_vars))
                        out <- quantile(herit, prob = c(0.025, 0.5, 0.975))
                        tibble(var = "heritability", mean = mean(herit), 
                                      `0.025quant` = out[[1]], `0.975quant` = out[[3]])
                        }) %>% 
        .[c(2,5)] %>% 
        bind_rows() %>% 
        mutate(response = c("hindleg", "weight"))
        

lambs_plot_df <- bind_rows(all_fixed, all_heritabilities) 

names(all_fixed) <- c("var", "mean", "sd", "lower_ci", "median", "upper_ci", "mode", "kld", "response")

all_fixed %>% 
        mutate(var = as.factor(var)) %>% 
        mutate(var = forcats::fct_rev(factor(var))) %>% 
ggplot(aes(mean, var)) + 
        geom_point() +
        #scale_y_discrete(limits = rev(levels(var))) +
        geom_errorbarh(aes(xmax = lower_ci, xmin= upper_ci), height = 0) +
        facet_wrap(~response, scales = "free_x") +
        theme_bw() +
        geom_vline(xintercept=0) +
        ylab("") +
        xlab("Estimate +- 95%CI / units in mm (hindleg) /kg (weight)") +
        theme(panel.spacing = unit(4, "lines")) -> plot_lambs
        
ggsave("figs/roh_lambs.jpg", plot_lambs, height = 3, width = 8)


rownames(test$summary.hyperpar)

herit <- (inla.rmarginal(1000, 1/(test$marginals.hyperpar[[4]])) / rowSums(map_df(test$marginals.hyperpar[1:4], function(x) 1/(inla.rmarginal(1000, x)))))
quantile(herit, prob = c(0.025, 0.5, 0.975))



# adults
mod_adults <- readRDS("output/models/mods_adults.rds")


all_fixed <- map(mod_adults$fitted, function(x) {
        out <- tibble::rownames_to_column(x$summary.fixed, var = "var")
        out %>% 
                filter(str_detect(var, "froh"))
}) %>% 
        .[c(1,6,7,12)] %>% 
        bind_rows() %>% 
        mutate(response = rep(c("hindleg", "weight"), each = 4))


names(all_fixed) <- c("var", "mean", "sd", "lower_ci", "median", "upper_ci", "mode", "kld", "response")

all_fixed %>% 
        mutate(var = as.factor(var)) %>% 
        mutate(var = forcats::fct_rev(factor(var))) %>% 
        ggplot(aes(mean, var)) + 
        geom_point() +
        #scale_y_discrete(limits = rev(levels(var))) +
        geom_errorbarh(aes(xmax = lower_ci, xmin= upper_ci), height = 0) +
        facet_wrap(~response, scales = "free_x") +
        theme_bw() +
        geom_vline(xintercept=0) +
        ylab("") +
        xlab("Estimate +- 95%CI / units in mm (hindleg) /kg (weight)") +
        theme(panel.spacing = unit(4, "lines")) -> plot_adults

ggsave("figs/roh_adults.jpg", plot_adults, height = 3, width = 8)


# adults
mod_lbs <- readRDS("output/models/mods_lbs.rds")


all_fixed <- map(mod_lbs$fitted, function(x) {
        out <- tibble::rownames_to_column(x$summary.fixed, var = "var")
        out %>% 
                filter(str_detect(var, "froh"))
}) %>% 
        .[c(1,2,4,5)] %>% 
        bind_rows() %>% 
        mutate(response = rep(c("F", "M"), each = 4))


names(all_fixed) <- c("var", "mean", "sd", "lower_ci", "median", "upper_ci", "mode", "kld", "sex")

all_fixed %>% 
        mutate(var = as.factor(var)) %>% 
        mutate(var = forcats::fct_rev(factor(var))) %>% 
        ggplot(aes(mean, var)) + 
        geom_point() +
        #scale_y_discrete(limits = rev(levels(var))) +
        geom_errorbarh(aes(xmax = lower_ci, xmin= upper_ci), height = 0) +
        facet_wrap(~sex, scales = "free_x") +
        theme_bw() +
        geom_vline(xintercept=0) +
        ylab("") +
        xlab("Estimate on logscale. \n A completely inbred female has exp(-2.54) ~ 8% of the lbs
of a completely outbred female. For males, its 0.02%") +
        theme(panel.spacing = unit(4, "lines")) -> plot_lbs
ggsave("figs/roh_lbs.jpg", plot_lbs, height = 3, width = 8)

