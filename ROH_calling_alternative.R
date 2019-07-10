library(tidyverse)
library(snpStats)
library(detectRUNS)

gen_map <- read_delim("../sheep/data/2_Soay_Linkage_Map_b.txt", delim = "\t")
head(gen_map)

genotypeFilePath <- "../sheep/data/SNP_chip/20140214_SheepHD_QC1_Polym.ped"
mapFilePath <- "../sheep/data/SNP_chip/20140214_SheepHD_QC1_Polym.map"

slidingRuns <- slidingRUNS.run(
        genotypeFile = genotypeFilePath, 
        mapFile = mapFilePath, 
        windowSize = 50, 
        threshold = 0.05,
        minSNP = 50, 
        ROHet = FALSE, 
        maxOppWindow = 1, 
        maxMissWindow = 5,
        maxGap = 10^8, 
        minLengthBps = 1500000, 
        minDensity = 1/10^3, # SNP/kbps
        maxOppRun = NULL,
        maxMissRun = NULL
) 
