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
library(brms)
library(MCMCglmm)
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

# get LRS
sheep_ped <- read_csv("../sheep_pedigree/Jisca_sequoia/sheeppedigree/Pedigree_SoaySheep_2019-07-03_toDB.csv") 

sheep_ped <- read_delim("../sheep/data/SNP_chip/20190208_Full_Pedigree.txt", delim = "\t")[c(1,3,2)] %>%
        as.data.frame() %>%
        orderPed() %>%
        rename(Mother = MOTHER, Father = FATHER)
#save(sheep_ped,  file = "model_in/sheep_ped.RData")

# females: offspring
sheep_ped %>% 
        group_by(Mother) %>% 
        tally() %>% 
        filter(!is.na(Mother)) %>% 
        rename(ID = Mother) -> female_LRT

# males: offspring
sheep_ped %>% 
        group_by(Father) %>% 
        tally() %>% 
        filter(!is.na(Father)) %>% 
        rename(ID = Father) -> male_LRT
# Inds with no offspring
no_offspring <- sheep_ped %>% 
      #  rename(ID = Code) %>% 
        mutate(n = ifelse( !((ID %in% .$Mother) | (ID %in% .$Father)), 0, NA)) %>% 
        filter(n == 0) %>% 
        dplyr::select(ID, n)

# put together
LRT_df <- rbind(female_LRT, male_LRT) %>% 
                rbind(no_offspring) %>% 
                rename(n_offspring = n)

Sheep$ID <- as.numeric(Sheep$ID) # for joining

# add other variables and filter out living individuals
LRT <- LRT_df %>% 
        left_join(Sheep, by="ID") %>% 
        left_join(tblPreg, by = "BirthRef") %>% 
        dplyr::select(ID, n_offspring, BirthRef, Status, Sex, DeathYear, BirthWt,
                      BirthYear, MumID) %>% 
        filter(Status %in% c("Dead", "EstDead")) %>% 
        mutate_at(c("ID", "MumID", "BirthYear", "DeathYear", "Sex"), as.character)


##### ROH data #####
file_path <- "output/ROH/roh_nofilt.hom"
file <- "roh_nofilt"
roh_lengths <- fread(file_path)

froh <- roh_lengths %>%
        dplyr::group_by(FID) %>%
        #filter((KB < 3130)) %>% 
        #filter((KB > 3130)&(KB < 12500)) %>% 
        #filter((KB > 12500)) %>% 
        #filter(KB < 3130) %>% 
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(MBSUM = KBSUM / 1000) %>%
        mutate(FROH = KBSUM / 2534344) %>%
        rename(ID = FID) %>% 
        mutate(ID = as.character(ID))

# dataset for modeling

# individuals with both lrt and roh info
lrt_roh_df <- LRT %>% 
               # mutate(ID = as.character(ID)) %>% 
                inner_join(froh, by = "ID") %>% 
                filter(!is.na(Sex))

# check distribution
lrt_roh_df %>% 
        filter(!is.na(Sex)) %>% 
        mutate(Sex = as.factor(Sex)) %>% 
        ggplot(aes(n_offspring)) + 
                geom_histogram(bins = 100) +
                scale_y_sqrt() + 
                facet_wrap(Sex~., scales = "free_x")

save(lrt_roh_df, "model_in/lrt_roh_df.RData")
save(sheep_ped, "model_in/sheep_ped.RData")
