library(boot)
# Modelling inbreeding depression in survival
# Using binomial mixed (animal) models with logit link, annual survival as response
# main modeling package is INLA
library(lme4)
library(tidyverse)
library(broom.mixed)
source("theme_simple.R")
library(INLA) # Downloaded from http://www.r-inla.org/download
library(AnimalINLA) # Downloaded from http://www.r-inla.org/related-projects/animalinla
library(MCMCglmm)
library(sjPlot)
library(brinla)
# data
load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")
ped <- sheep_ped

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Annual survival~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# survival data preprocessing
annual_survival <- fitness_data %>% 
        # filter na rows
        filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
        mutate(age_cent = age - mean(age, na.rm = TRUE),
               age_cent2 = age_cent^2,
               age_std = as.numeric(scale(age)),
               age_std2 = age_std^2,
               froh_all_cent = froh_all - mean(froh_all, na.rm = TRUE),
               # times 10 to estimate a 10% percent increase
               froh_all10 = froh_all * 10,
               froh_all10_cent = froh_all10 - mean(froh_all10, na.rm = TRUE),
               lamb = ifelse(age == 0, 1, 0),
               lamb_cent = lamb - mean(lamb, na.rm = TRUE),
               lamb = as.factor(lamb)) %>% 
        as.data.frame() 


ped <- sheep_ped[, c(1,3,2)]
ped_fin <- prePed(ped)
ID <- pedInbreeding(ped_fin) %>% as_tibble() %>% rename(id = Indiv, fped = Inbr)

froh <- annual_survival %>% select(id, froh_all)
f_vals <- froh %>% left_join(ID) %>% distinct()

ggplot(f_vals, aes(fped, froh_all)) +
        geom_point()


annual_survival %>% filter(froh_all > 0.34) %>% nrow(.)
annual_survival %>% filter(froh_all > 0.24) %>% nrow(.)

# average age according to Suppl. Figure 9
av_age <- (2134*1 + 1589*2 + 1282*3 + 1100*4 + 910*5 + 720*6 + 554*7 + 422*8 + 254*9)/(2134 + 1589 + 1282 + 1100 + 910 + 720 + 554 + 422 + 254)

# parameter estimates in linear predictor according to Suppl. Table 4 for males (sex=1) and not twins (twin=0)
int <- 3.34
Froh <- -1.14
Age <- -0.21
Lamb <- -3.41
Sex <- -0.63
fa_int <- 0.17
fl_int <- 0.62




# calculate linear predictor based on the parameter estimates above

# linear predictor for lamb (age=0) at average F_ROH (=0)
linpred_L_0 <- int + Froh*0 + Age*(0-av_age) + Lamb + Sex*1 + fa_int*(0-av_age)*0 + fl_int*0*1
surv_L_0 <- inv.logit(linpred_L_0)
surv_L_0

# linear predictor for lamb (age=0) at F_ROH of 0.34 (=1)
linpred_L_1 <- int + Froh*1 + Age*(0-av_age) + Lamb + Sex*1 + fa_int*(0-av_age)*1 + fl_int*1*1
surv_L_1 <- inv.logit(linpred_L_1)
surv_L_1

linpred_L_2 <- int + Froh*2 + Age*(0-av_age) + Lamb + Sex*1 + fa_int*(0-av_age)*1 + fl_int*1*1
surv_L_2 <- inv.logit(linpred_L_2)
surv_L_2


# inbreeding depression lambs
ID_Lamb <- (surv_L_0 - surv_L_1)/surv_L_0
# odds ratio for inbreeding effect in lambs
OR_Lamb <- (surv_L_1/(1-surv_L_1)/(surv_L_0/(1-surv_L_0)))

# linear predictor for age=1 at average F_ROH (=0)
linpred_1_0 <- int + Froh*0 + Age*(1-av_age) + Sex*1+ fa_int*(1-av_age)*0
surv_1_0 <- inv.logit(linpred_1_0)
surv_1_0

# linear predictor for age=1 at F_ROH of 0.34 (=1)
linpred_1_1 <- int + Froh*1 + Age*(1-av_age) + Sex*1 + fa_int*(1-av_age)*1
surv_1_1 <- inv.logit(linpred_1_1)
surv_1_1

# inbreeding depression age=1
ID_1 <- (surv_1_0 - surv_1_1)/surv_1_0
# odds ratio for inbreeding effect in sheep age=1
OR_1 <- (surv_1_1/(1-surv_1_1)/(surv_1_0/(1-surv_1_0)))

# linear predictor for age=4 at average F_ROH (=0)
linpred_4_0 <- int + Froh*0 + Age*(4-av_age) + Sex*1 + fa_int*(4-av_age)*0
surv_4_0 <- inv.logit(linpred_4_0)
surv_4_0

# linear predictor for age=4 at F_ROH of 0.34 (=1)
linpred_4_1 <- int + Froh*1 + Age*(4-av_age) + Sex*1 + fa_int*(4-av_age)*1
surv_4_1 <- inv.logit(linpred_4_1)
surv_4_1

# inbreeding depression age=4
ID_4 <- (surv_4_0 - surv_4_1)/surv_4_0
# odds ratio for inbreeding effect in sheep age=4
OR_4 <- (surv_4_1/(1-surv_4_1)/(surv_4_0/(1-surv_4_0)))

# linear predictor for age=7 at average F_ROH (=0)
linpred_7_0 <- int + Froh*0 + Age*(7-av_age) + Sex*1 + fa_int*(7-av_age)*0
linpred_7_0
surv_7_0 <- inv.logit(linpred_7_0)
surv_7_0

# linear predictor for age=7 at F_ROH of 0.34 (=1)
linpred_7_1 <- int + Froh*1 + Age*(7-av_age) + Sex*1 + fa_int*(7-av_age)*1
linpred_7_1
surv_7_1 <- inv.logit(linpred_7_1)
surv_7_1

# inbreeding depression age=7
ID_7 <- (surv_7_0 - surv_7_1)/surv_7_0
# odds ratio for inbreeding effect in sheep age=2
OR_7 <- (surv_7_1/(1-surv_7_1)/(surv_7_0/(1-surv_7_0)))

# compare inbreeding depression and odds ratios
ID_Lamb
OR_Lamb
ID_1
OR_1
ID_4
OR_4
ID_7
OR_7


# Estimate age-specific beta_FROH (as in Suppl. Fig. 9) and plot them

beta_FROH <- rep(NA,10)

beta_FROH[1] <- Froh + fl_int + fa_int*(0-av_age)
beta_FROH[2] <- Froh + fa_int*(1-av_age)
beta_FROH[3] <- Froh + fa_int*(2-av_age)
beta_FROH[4] <- Froh + fa_int*(3-av_age)
beta_FROH[5] <- Froh + fa_int*(4-av_age)
beta_FROH[6] <- Froh + fa_int*(5-av_age)
beta_FROH[7] <- Froh + fa_int*(6-av_age)
beta_FROH[8] <- Froh + fa_int*(7-av_age)
beta_FROH[9] <- Froh + fa_int*(8-av_age)
beta_FROH[10] <- Froh + fa_int*(9-av_age)

plot(0:9,beta_FROH,xlab='age')
