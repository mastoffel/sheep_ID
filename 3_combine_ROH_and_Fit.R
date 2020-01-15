library(tidyverse)
library(MasterBayes)
library(data.table)
library(rlang)

# annual measures of traits and fitness
fitness_path <- "../sheep/data/1_Annual_Fitness_Measures_April_20190501.txt"
annual_fitness <- read_delim(fitness_path, delim = "\t")
names(annual_fitness)

# prepare, order pedigree
sheep_ped <- read_delim("../sheep/data/SNP_chip/20190711_Soay_Pedigree.txt", 
                        delim = "\t",
                        col_types = "ccc") %>%
        as.data.frame() %>%
        MasterBayes::orderPed() 
#save(sheep_ped,  file = "model_in/sheep_ped.RData")

# IDs which weren't well imputed and should be removed
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

##### ROH data #####
file_path <- "output/ROH/roh_nofilt_ram_pruned.hom"
file <- "roh_nofilt_pruned"
roh_lengths <- fread(file_path)
hist(roh_lengths$KB, breaks = 1000, xlim = c(0,10000))

# Chr lengths
chr_data <- read_delim("../sheep/data/sheep_genome/chromosome_info_ram.txt", delim = "\t") %>% 
  rename(size_BP = Length,
         CHR = Part) %>% 
  mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
                          .[2:27, ] %>% 
                          summarise(sum_KB = sum(size_KB)) %>% 
                          as.numeric()

# define ROH length classes
calc_froh_classes <- function(roh_crit, roh_lengths) {
        
        roh_filt <- dplyr::case_when(
                roh_crit == "short"  ~ expr(KB < 1221),
                roh_crit == "medium" ~ expr((KB > 1221)&(KB < 4885)),
                roh_crit == "long"   ~ expr(KB > 4885),
                roh_crit == "all" ~ expr(KB > 0)
        )
        
        roh_lengths %>%
                dplyr::group_by(IID) %>%
                #filter({{ roh_filt }}) %>% 
                filter(!!roh_filt) %>% 
                dplyr::summarise(KBSUM = sum(KB)) %>% 
                mutate(FROH = KBSUM / autosomal_genome_size) %>% 
                dplyr::select(IID, FROH) %>% 
                rename(ID = IID, !! paste0("FROH_", roh_crit) := FROH)
        
}

# proportion of ROH length classes in each genome. Individuals which
# do not have long ROH have 0 for this class.
ROH_classes <- c("short", "medium", "long", "all")
froh <- purrr::map(ROH_classes, calc_froh_classes, roh_lengths) %>% 
        purrr::reduce(left_join, by = "ID") %>% 
        replace_na(list(FROH_long = 0))

# add FROH minus each of the chromosomes as new variabels for gwas
calc_froh_minus_chr <- function(chr, roh_lengths) {
    
    # get chromosome length and substract from autosomal genome length
    chr_size <- chr_data %>% 
                filter(CHR == paste0("Chromosome ", chr)) %>% 
                .[["size_KB"]]
    genome_size_minus_chr <- autosomal_genome_size - chr_size
    
    roh_lengths %>%
      as_tibble() %>% 
      dplyr::group_by(IID) %>%
      #filter({{ roh_filt }}) %>% 
      dplyr::filter(CHR != !!chr) %>% 
      dplyr::summarise(KBSUM = sum(KB)) %>%
      mutate(FROH = KBSUM / genome_size_minus_chr ) %>%
      dplyr::select(IID, FROH) %>%
      rename(ID = IID, !! paste0("froh_no_chr", chr) := FROH) 
  
}

froh_no_chr <- map(1:26, calc_froh_minus_chr, roh_lengths) %>% 
               reduce(left_join, by = "ID")

# add to froh
froh <- froh %>% 
        left_join(froh_no_chr, by = "ID")

# add froh to fitness
fitness_data <- annual_fitness %>% 
        left_join(froh, by = "ID")

# add homozygosity not in roh
homs <- read_delim("output/ROH/roh_nofilt_hom_not_in_roh.txt", delim = " ")
fitness_data <- fitness_data %>% 
        left_join(homs, by = "ID")

# prepare for analysis
fitness_data <- fitness_data %>% 
        dplyr::rename(birth_year = BIRTHYEAR,
                sheep_year = SheepYear,
                age = Age,
                id = ID,
                twin = TWIN,
                sex = SEX,
                mum_id = MOTHER,
                froh_short = FROH_short,
                froh_medium = FROH_medium,
                froh_long = FROH_long,
                froh_all = FROH_all,
                froh_not_roh = hom,
                survival = Survival,
                comment = Comment,
                seen_in_rut = SeenInRut,
                dad_id = FATHER,
                weight = Weight, 
                birth_weight = BIRTHWT,
                cap_month = CapMonth,
                hindleg = Hindleg,
                foreleg = Foreleg,
                rut_measure = RutMeasure,
                horn_len = HornLen,
                horn_circ = HornCirc,
                bol_circ = BolCirc,
                bol_len = BolLen,
                horn = Horn,
                freq = Freq,
                offspr_born = OffspringBorn,
                offspr_surv = OffspringSurvived,
                offspr_surv_aug = AugSurvOffspring,
                offspr_surv_oct = OctSurvOffspring,
                offspr_surv_wint = OverWinterOffspring) %>% 
        dplyr::select(id, survival, comment, sheep_year, age, birth_year,
                      sex, mum_id, twin, froh_all, froh_long, froh_medium, 
                      froh_all, everything()) %>% 
        # for sex, check what Cas is
        filter(sex %in% c("F", "M")) %>% 
        # some individuals arent imputed well and should be discarded 
        filter(!(id %in% IDs_lots_missing$id)) %>% 
        # na rows should be discarded
        mutate_at(c("id", "sheep_year", "birth_year", "sex",
                    "mum_id", "twin"), as.factor)

# save data
save(fitness_data, file = "data/survival_mods_data.RData") # formerly fitness_roh_df
save(sheep_ped, file = "data/sheep_ped.RData") # ordered ped


