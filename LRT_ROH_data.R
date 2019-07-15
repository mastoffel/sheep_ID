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
library(rlang)

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
#sheep_ped <- read_csv("../sheep_pedigree/Jisca_sequoia/sheeppedigree/Pedigree_SoaySheep_2019-07-03_toDB.csv") 

sheep_ped <- read_delim("../sheep/data/SNP_chip/20190208_Full_Pedigree.txt", delim = "\t")[c(1,3,2)] %>%
        as.data.frame() %>%
        orderPed() %>%
        rename(Mother = MOTHER, Father = FATHER)
#save(sheep_ped,  file = "model_in/sheep_ped.RData")

# annual fitness
annual_fitness <- read_delim("data/1_Annual_Fitness_Measures_April_20190501.txt", delim = "\t")
annual_fitness %<>% 
  dplyr::select(-contains("EASTING")) %>% 
  dplyr::select(-contains("NORTHING")) %>% 
  dplyr::select(-contains("Cartesian"))

##### ROH data #####
file_path <- "output/ROH/roh_nofilt.hom"
file <- "roh_nofilt"
roh_lengths <- fread(file_path)


# roh_crit = c("short", "medium", "long", "all")
calc_froh_classes <- function(roh_crit, roh_lengths) {
  
  roh_filt <- dplyr::case_when(
    roh_crit == "short"  ~ expr(KB < 3130),
    roh_crit == "medium" ~ expr((KB > 3130)&(KB < 12500)),
    roh_crit == "long"   ~ expr(KB > 12500),
    roh_crit == "all" ~ expr(KB > 0)
  )
  
  roh_lengths %>%
    dplyr::group_by(FID) %>%
    #filter({{ roh_filt }}) %>% 
    filter(!!roh_filt) %>% 
    dplyr::summarise(KBSUM = sum(KB)) %>% 
    mutate(FROH = KBSUM / 2534344) %>% 
    dplyr::select(FID, FROH) %>% 
    rename(ID = FID, !! paste0("FROH_", roh_crit) := FROH)
  
}

# proportion of ROH length classes in each genome. Individuals which
# do not have long ROH have 0 for this class.
froh <- purrr::map(c("short", "medium", "long", "all"), calc_froh_classes,  roh_lengths) %>% 
          purrr::reduce(left_join, by = "ID") %>% 
          replace_na(list(FROH_long = 0))


# dataset for modeling
fitness_data <- annual_fitness %>% 
  dplyr::select(ID, OffspringBorn, OffspringSurvived, SheepYear, Age, Survival,
         BIRTHYEAR, SEX, MOTHER, BIRTHWT, CapMonth, Weight, Hindleg, ID.CapYear) %>% 
  left_join(froh, by = "ID")

save(fitness_data, file = "model_in/fitness_roh_df.RData")
save(sheep_ped, file = "model_in/sheep_ped.RData")



# # females: offspring
# sheep_ped %>% 
#         group_by(Mother) %>% 
#         tally() %>% 
#         filter(!is.na(Mother)) %>% 
#         rename(ID = Mother) -> female_LRT
# 
# # males: offspring
# sheep_ped %>% 
#         group_by(Father) %>% 
#         tally() %>% 
#         filter(!is.na(Father)) %>% 
#         rename(ID = Father) -> male_LRT
# # Inds with no offspring
# no_offspring <- sheep_ped %>% 
#       #  rename(ID = Code) %>% 
#         mutate(n = ifelse( !((ID %in% .$Mother) | (ID %in% .$Father)), 0, NA)) %>% 
#         filter(n == 0) %>% 
#         dplyr::select(ID, n)
# 
# # put together
# LRT_df <- rbind(female_LRT, male_LRT) %>% 
#                 rbind(no_offspring) %>% 
#                 rename(n_offspring = n)
# 
# Sheep$ID <- as.numeric(Sheep$ID) # for joining
# 
# # add other variables and filter out living individuals
# LRT <- LRT_df %>% 
#         left_join(Sheep, by="ID") %>% 
#         left_join(tblPreg, by = "BirthRef") %>% 
#         dplyr::select(ID, n_offspring, BirthRef, Status, Sex, DeathYear, BirthWt,
#                       BirthYear, MumID) %>% 
#         filter(Status %in% c("Dead", "EstDead")) %>% 
#         mutate_at(c("ID", "MumID", "BirthYear", "DeathYear", "Sex"), as.character)
# individuals with both lrt and roh info
# lrt_roh_df <- LRT %>% 
#   # mutate(ID = as.character(ID)) %>% 
#   inner_join(froh, by = "ID") %>% 
#   filter(!is.na(Sex))
# 
# # check distribution
# lrt_roh_df %>% 
#   filter(!is.na(Sex)) %>% 
#   mutate(Sex = as.factor(Sex)) %>% 
#   ggplot(aes(n_offspring)) + 
#   geom_histogram(bins = 100) +
#   scale_y_sqrt() + 
#   facet_wrap(Sex~., scales = "free_x")
# 
# save(lrt_roh_df, "model_in/lrt_roh_df.RData")
# save(sheep_ped, "model_in/sheep_ped.RData")
