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
        filter(Status %in% c("Dead", "EstDead"))


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
                mutate(ID = as.character(ID)) %>% 
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

##### Prepare pedigree for asreml #####
# pedigree format: Needs to be sorted, and individual, father, mother

# lrt_roh_df_m <- lrt_roh_df %>% filter(Sex == 1) %>% setDF() %>% 
#                         mutate(MumID = as.factor(MumID),
#                                BirthYear = as.factor(BirthYear),
#                                DeathYear = as.factor(DeathYear))

# sheep_ped_asreml <- sheep_ped[, c(1,3,2)] %>% 
#                         rename(ID = Code) %>% 
#                         filter(!is.na(ID)) %>% 
#                         #filter(ID %in% lrt_roh_df$ID) %>% 
#                         as.data.frame() %>% 
#                         orderPed() 

ggplot(lrt_roh_df, aes(FROH, n_offspring)) +
        geom_smooth(method = "lm") +
        geom_point(size = 0.1, alpha = 0.3) +
        theme_classic() +
        facet_wrap(Sex~., scales = "free_y") +
        scale_y_sqrt()


# brms
library(brms)

lrt_roh_df <- lrt_roh_df %>% 
                        rename(animal = ID)

inv.phylo <- MCMCglmm::inverseA(sheep_ped)
sheep_ped_inv <- solve(inv.phylo$Ainv)
rownames(sheep_ped_inv) <- rownames(inv.phylo$Ainv)

#save(lrt_roh_df, file = "model_in/lrt_roh_df.RData")
#save(sheep_ped_inv, file = "model_in/sheep_ped_inv.RData")

lrt_roh_df_mod1 <- lrt_roh_df %>% 
        filter(Sex == 1) %>% 
        filter(!is.na(MumID)) %>% 
        filter(!is.na(BirthYear)) %>% 
        filter(!is.na(DeathYear)) %>% 
        sample_frac(0.05) 
        
model1 <- brm(
        n_offspring ~ 1 + FROH + (1|MumID) + (1|BirthYear) + (1|DeathYear) + (1|animal), data = lrt_roh_df_mod1 , 
        family = negbinomial(), cov_ranef = list(animal = sheep_ped_inv),
        chains = 2, cores = 2
        #control = list(adapt_delta = 0.9)
)

model2 <- brm(
        n_offspring ~ 1 + FROH + (1|MumID) + (1|BirthYear) + (1|DeathYear) + (1|animal), data = lrt_roh_df_mod1 , 
        family = poisson(), cov_ranef = list(animal = A),
        chains = 2, cores = 2
        #control = list(adapt_delta = 0.9)
)

model3 <- brm(
        n_offspring ~ 1 + FROH + (1|MumID) + (1|BirthYear) + (1|DeathYear) + (1|animal), data = lrt_roh_df_mod1 , 
        family = hurdle_negbinomial(), cov_ranef = list(animal = A),
        chains = 2, cores = 2
        #control = list(adapt_delta = 0.9)
)






summary(model2)
plot(model2)
coef(model2)

pp_check(model2)


# inverse relationship mat
sheep_ainv <- asreml::ainverse(sheep_ped)

#attr(sheep_ainv, "rowNames")
# 1 + FROH, random = ~vm(sheep_ped) + idv(MumID) + idv(BirthYear) + idv(DeathYear)
mod <- asreml(fixed = n_offspring ~ 1 + FROH, random = ~vm(sheep_ainv) + idv(MumID), 
              data = lrt_roh_df_mod,
              #family = asr_negative.binomial(link = "log", dispersion = 1, phi = NA),
              #   family = list(asr_poisson()),
              na.action = na.method(x=c("omit")))

mod_sum <- summary(mod, coef = TRUE)
mod_sum$varcomp
mod_sum$coef.fixed
