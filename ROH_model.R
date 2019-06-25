library(sheepDB)
library(RJDBC)
library(dplyr)
library(dbplyr)
library(lubridate)
library(asreml)
library(tidyverse)
library(data.table)
library("MasterBayes")
library(pedigreemm)
#### database connection ######

dbname <- "../sheep/data/db/StKilda_Data.accdb"
driver <- "net.ucanaccess.jdbc.UcanloadDriver"
driverpath <- "../sheep/data/db/UCanAccess/loader/ucanload.jar"
options <- paste0("jdbc:ucanaccess://", dbname, ";memory=false")

con <- DBI::dbConnect(JDBC(driver, driverpath), options)
# src <- src_dbi(con)

tbls <- dbGetTables(con)
flds <- dbGetFields(con, "Sheep")
Sheep <- dbGetQuery(con, "Select * from Sheep")
Census <- dbGetQuery(con, "Select * from CensusData")
Capture <- dbGetQuery(con, "Select * from CaptureData")
tblPreg <- dbGetQuery(con, "Select * from tblPregnancies")
dbDisconnect(con)


#### filter and assemble data #####

hindleg_df <- Capture %>% 
        #filter(!is.na(Hindleg)) %>% 
        # only august
        filter(CapMonth == 8) %>% 
        left_join(Sheep, by = "ID") %>% 
        left_join(tblPreg, by = "BirthRef") %>% 
        dplyr::select(ID, CapYear,CapDay, Weight, Hindleg, Sex, BirthWt, Coat, MumID, BirthYear) %>% 
        # only traits in august, birth year (~4 month old)
        filter(CapYear == BirthYear) %>% 
        
       # filter(!is.na(Hindleg)) %>% 
        filter(!is.na(Weight)) %>% 
        group_by(ID) %>% 
        arrange(CapDay) %>% 
        # if multiple catches, take the first
        top_n(n = 1, CapDay)


##### ROH data #####
file_path <- "output/ROH/roh_nofilt.hom"
file <- "roh_nofilt"
roh_lengths <- fread(file_path)



ROH_mod <- function(ROH_criterion) {
        
froh <- roh_lengths %>%
        dplyr::group_by(FID) %>%
        # all ROH greater than 5Mb
        filter(!!ROH_criterion) %>% 
        #filter((KB < 3130)) %>% 
        #filter((KB > 3130)&(KB < 12500)) %>% 
        #filter((KB > 12500)) %>% 
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(MBSUM = KBSUM / 1000) %>% 
        mutate(FROH = KBSUM/2869490) %>% 
        rename(ID = FID)

# create data.frame for modeling
df_ROH_hindleg <- left_join(hindleg_df , froh, by = "ID")

##### Prepare pedigree for asreml #####
# pedigree format: Needs to be sorted, and individual, father, mother
sheep_ped <- read_delim("../sheep/data/SNP_chip/20190208_Full_Pedigree.txt", delim = "\t")[c(1,3,2)] %>% 
                        as.data.frame() %>% 
                        filter(ID %in% df_ROH_hindleg$ID) %>% 
                        orderPed() 
                        # rename(Calf = ID,
                        #        Sire = FATHER,
                        #        Dam = MOTHER)
# inverse relationship mat
sheep_ainv <- asreml::ainverse(sheep_ped)

df_ROH_hindleg <- df_ROH_hindleg %>% 
        ungroup() %>% 
        mutate(ID = as.factor(ID),
               MumID = as.factor(MumID),
               BirthYear = as.factor(BirthYear))

##### modeling ######
# random effects: birth year, maternal ID, additive genetic
# fixed: sex, FROH
# ggplot(data = df_ROH_hindleg, aes(FROH, Hindleg, color = Sex)) + geom_point() + geom_smooth(method="lm")

mod <- asreml(fixed = Weight ~ 1 + FROH + Sex, random = ~vm(ID, sheep_ainv) + idv(MumID) + idv(BirthYear), 
               data = df_ROH_hindleg, na.action = na.method(x=c("omit"))) # "omit" "include"
# Hindleg
mod_sum <- summary(mod, coef = TRUE)
mod_sum$varcomp
mod_sum$coef.fixed

out <- mod_sum$coef.fixed %>% 
        as.data.frame() %>% 
        rownames_to_column("var") %>% 
        filter(var == "FROH") 

}

crits <- list(expr(KB < 3130), expr((KB > 3130)&(KB < 12500)), expr(KB > 12500))
        
library(purrr)
all_mods <- purrr::map_df(crits, ROH_mod)
# all_mods$run <- c("ROH < 3Mb", "3Mb < ROH < 12.5Mb", "ROH > 12.5Mb ")
all_mods$run <- c("short", "medium", "long")
mod_plot <- all_mods %>% 
        rename(std_error = "std error") %>% 
        mutate(err_low = solution - std_error/2,
               err_high = solution + std_error/2) %>% 
        mutate(CI_low = solution - 1.96*std_error,
               CI_high = solution + 1.96*std_error)

source("theme_clean.R")
p1 <- ggplot(data = mod_plot, aes(solution, run, xmax = err_low, xmin = err_high)) +
        geom_point(size = 4, shape = 21, col = "black", fill = "grey69",
                   alpha = 0.7) +
        geom_errorbarh(alpha = 1, color = "black", height = 0) +
        theme_clean() +
        ylab("ROH length class") +
        xlab("Estimate (SE)") +
        xlim(-15, 5) +
        ggtitle("ROH and 4-month weight") +
        geom_vline(xintercept = 0, color = "grey80")

p2 <- ggplot(data = mod_plot, aes(solution, run, xmax = CI_low, xmin = CI_high)) +
        geom_point(size = 4, shape = 21, col = "black", fill = "grey69",
                   alpha = 0.7) +
        geom_errorbarh(alpha = 1, color = "black", height = 0) +
        theme_clean() +
        ylab("ROH length class") +
        xlab("Estimate (CI)") +
        xlim(-20, 10) +
        geom_vline(xintercept = 0, color = "grey80")

library(cowplot)
p_hindleg <- plot_grid(p1, p2)
p_hindleg
ggsave("figs/ROH_weight.jpg", plot = p_hindleg, width = 7, height = 3)


library(asremlPlus)
simulate(mod, nsim = 10)
