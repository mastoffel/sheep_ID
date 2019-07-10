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

# get LRS
sheep_ped <- read_csv("../sheep_pedigree/Jisca_sequoia/sheeppedigree/Pedigree_SoaySheep_2019-07-03_toDB.csv") 

sheep_ped %>% 
        group_by(Mother) %>% 
        tally() %>% 
        filter(!is.na(Mother)) %>% 
        rename(ID = Mother) -> female_LRT
sheep_ped %>% 
        group_by(Father) %>% 
        tally() %>% 
        filter(!is.na(Father)) %>% 
        rename(ID = Father) -> male_LRT

no_offspring <- sheep_ped %>% 
        rename(ID = Code) %>% 
        mutate(n = ifelse( !((ID %in% .$Mother) | (ID %in% .$Father)), 0, NA)) %>% 
        filter(n == 0) %>% 
        dplyr::select(ID, n)

LRT_df <- rbind(female_LRT, male_LRT) %>% 
                rbind(no_offspring) %>% 
                rename(n_offspring = n)

Sheep %>% 
        left_join(LRT_df, by )


