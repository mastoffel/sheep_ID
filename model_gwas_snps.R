library(lme4)
library(performance)
library(tidyverse)
library(data.table)
library(snpStats)
library(performance)
early_survival_top_snps <- read_delim("output/early_survival_top_snps.txt", delim = " ") %>% 
        mutate(roh_all = rowSums(select_at(., vars(contains("roh"))), na.rm = TRUE))

# model without roh
glme_form <- reformulate(c("sex", "twin", "age_std", "age2_std", "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
mod1 <- glmer(glme_form, data = early_survival_top_snps, family = "binomial")
summary(mod1)
plot(mod1)
r2_nakagawa(mod1)

# model with roh
glme_form2 <- reformulate(c("sex", "twin", "age_std", "age2_std", "roh_all", "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
mod2 <- glmer(glme_form2, data = early_survival_top_snps, family = "binomial")
summary(mod2)
r2_nakagawa(mod2)
