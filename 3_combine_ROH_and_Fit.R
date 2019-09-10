library(tidyverse)
library(MasterBayes)
library(data.table)
# annual measures of traits and fitness
annual_fitness <- read_delim("../sheep/data/1_Annual_Fitness_Measures_April_20190501.txt", delim = "\t")
# sheep_ped <- annual_fitness %>% 
#         dplyr::select(ID, MOTHER, FATHER) %>% 
#         group_by(ID) %>% 
#         summarise(Mother = first(MOTHER),
#                   Father = first(FATHER))

# get LRS
sheep_ped <- read_delim("../sheep/data/SNP_chip/20190711_Full_Pedigree.txt", 
                        delim = "\t",
                        col_types = "ccc") %>%
        as.data.frame() %>%
        MasterBayes::orderPed() 
#save(sheep_ped,  file = "model_in/sheep_ped.RData")

##### ROH data #####
file_path <- "output/ROH/roh_nofilt.hom"
file <- "roh_nofilt"
roh_lengths <- fread(file_path)
hist(roh_lengths$KB, breaks = 1000, xlim = c(0,10000))

# roh_crit = c("short", "medium", "long", "all")
calc_froh_classes <- function(roh_crit, roh_lengths) {
        
        roh_filt <- dplyr::case_when(
                roh_crit == "short"  ~ expr(KB < 2443),
                roh_crit == "medium" ~ expr((KB > 2443)&(KB < 9771)),
                roh_crit == "long"   ~ expr(KB > 9771),
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

# all FROH are negatively correlated
library(GGally)
ggpairs(froh[-1])

# add death year
# load("model_in/lrt_roh_df.RData")
# lrt_roh_df <- lrt_roh_df %>% rename(ID = animal) %>% 
#         dplyr::select(ID, DeathYear) %>% 
#         mutate(ID = as.numeric(as.character(ID)))

# dataset for modeling
fitness_data <- annual_fitness %>% 
       # dplyr::select(ID, OffspringBorn, OffspringSurvived, SheepYear, Age, Survival,
      #                   BIRTHYEAR, SEX, MOTHER, BIRTHWT, CapMonth, Weight, Hindleg, ID.CapYear) %>% 
       # rename() %>% 
        left_join(froh, by = "ID")# %>% 
     #   left_join(lrt_roh_df, by = "ID")

# add homozygosity not in roh
homs <- read_delim("output/ROH/roh_nofilt_hom_not_in_roh.txt", delim = " ")

fitness_data <- fitness_data %>% 
        left_join(homs, by = "ID")

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
