library(sheepDB)
library(RJDBC)
library(dplyr)
library(dbplyr)
library(lubridate)
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
dbDisconnect(con)


Sheep[5000:5010, ]
Capture[5000:5010, ]

capt_sub <- Capture %>% 
        filter(!is.na(Hindleg)) %>% 
        filter(CapMonth == 8) %>% 
        group_by(ID) %>% 
        #arrange(ID, CapYear) %>% 
        top_n(-1, CapYear) %>% 
        filter(Weight < 21)
Census

hist(capt_sub$Weight, breaks = 500)

library(data.table)
file_path <- "output/ROH/roh_nofilt.hom"
file <- "roh_nofilt"
roh_lengths <- fread(file_path)

froh <- roh_lengths %>%
        dplyr::group_by(FID) %>%
        filter(KB < 5000) %>% 
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/2869490) %>% 
        rename(ID = FID)


df_mod <- left_join(capt_sub, froh, by = "ID")

library(ggplot2)
ggplot(df_mod, aes(Hindleg, FROH)) +
        geom_point() +
        geom_smooth(method = "lm")


       


