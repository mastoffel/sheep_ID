library(sheepDB)
library(tidyverse)
?get_sheep_data
library(lme4)
library(lme4qtl)
library(data.table)
library(brms)
library(MCMCglmm)
library("MasterBayes")
library(broom)
# list of tables
all_tables <- get_sheep_data(dbpath = "../sheep/data/db/", list_tables=TRUE)
# get all of them
sheep_dat <- get_sheep_data(dbpath = "../sheep/data/db/",
               dataset = all_tables)
# transform all into tibbles
sheep_dat <- map(sheep_dat, as_tibble)

# for hindleg length
table(sheep_dat$CaptureData$CapMonth)

names(sheep_dat)
sheep_dat$Sheep

# birth date
sheep_dat$zzz_SheepObsoleteBirthFieldsBackup

# extract 4-month hindleg length
hindleg_dat <- sheep_dat$CaptureData %>% 
        filter(CapMonth == 8) %>% 
        left_join(sheep_dat$zzz_SheepObsoleteBirthFieldsBackup, by = "ID") %>% 
        filter(CapYear == BirthYear) %>% 
        left_join(sheep_dat$Sheep, by = "ID") %>% 
        select(ID, Hindleg, Sex)

weight_dat <- sheep_dat$CaptureData %>% 
        filter(CapMonth == 8) %>% 
        left_join(sheep_dat$zzz_SheepObsoleteBirthFieldsBackup, by = "ID") %>% 
        filter(CapYear == BirthYear) %>% 
        left_join(sheep_dat$Sheep, by = "ID") %>% 
        select(ID, Weight, Sex, CapYear)

# get roh data
file_path <- "output/ROH/roh_nofilt_new.hom"
file <- "roh_nofilt"
roh_lengths <- fread(file_path)

froh <- roh_lengths %>%
        #filter(KB > 1000) %>% 
        dplyr::group_by(FID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/2869490) %>% 
        rename(ID = FID)

# combine 
sheep <- weight_dat %>% left_join(froh) %>% 
        filter(!is.na(Weight)) %>% 
        filter(!is.na(FROH)) 

hist(sheep$Weight)

sheep %>% 
        ggplot(aes(Weight, FROH)) +
        geom_point() +
        geom_smooth(method = "lm") +
                facet_wrap(~Sex)
length(unique(sheep$ID))


# modeling
sheep_ped <- read_delim("../sheep/data/SNP_chip/20190208_Full_Pedigree.txt", delim = '\t')
names(sheep_ped) <- c("animal", "dam", "sire")
sheep_ped <- orderPed(as.data.frame(sheep_ped))
#any(!(sheep_ped$dam[!is.na(sheep_ped$dam)] %in% sheep_ped$animal))
##Prepare the relationship matrix
inv.A <- MCMCglmm::inverseA(as.data.frame(sheep_ped))
A <- solve(inv.A$Ainv)
rownames(A) <- rownames(inv.A$Ainv)

# subset
sheep_mod <- sheep %>% 
        #sample_n(2000) %>% 
        rename(ANIMAL = ID)

m1 <- relmatLmer(Weight ~ FROH + Sex + (1|CapYear) + (1|ANIMAL), data = sheep_mod, relmat = list(ANIMAL = A))

m1 <- brm(Weight ~ FROH + Sex + (1|CapYear) + (1|ANIMAL), data = sheep_mod, 
          family = "gaussian", cov_ranef = list(ANIMAL = A), chains = 2, cores = 2,
          file = "output/brms_models/mod1")

summary(m1)

plot(m1)
library(lme4)
mod <- lmer()
summary(lm(Weight ~ FROH + Sex, data = sheep))

hist(roh_lengths$NSNP, breaks = 1000, xlim = c(0,2000))
hist(roh_lengths$KB)
roh_lengths[which.max(roh_lengths$KB), ]
hist(roh_lengths$PHET)


sheep <- froh %>% mutate(fetus = ifelse(ID %in% sheep_dat$tblFoetusMeasures$ID, "yes", "no"))

ggplot(sheep, aes(fetus, FROH)) +
        geom_boxplot() +
        geom_point(size = 4)

sheep_fetus <- froh %>% 
        filter 
        mutate(fetus = "yes")
        
non_fetus <- froh %>% 
        filter(!(ID %in% sheep_fetus$ID)) %>% 
        mutate(fetus = "no")

sheep_full <- rbind(sheep_fetus, non_fetus)
