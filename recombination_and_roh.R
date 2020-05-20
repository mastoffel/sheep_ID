# recombination rate and roh

library(data.table)
library(tidyverse)
library(windowscanr)
snp_roh <- fread("output/ROH/roh_ram_long.hom.summary")
lmap <- read_delim("../sheep_ID_paper/comments/susie/7_20200504_Full_Linkage_Map.txt", "\t") %>% 
        rename(SNP = SNP.Name)


snp_df <- snp_roh %>% 
                inner_join(lmap, by = "SNP") %>% 
                mutate(KB = BP/1000)

ggplot(snp_df, aes(r, cMdiff)) +
        geom_point() +
        geom_smooth(method = "lm") +
        facet_wrap(~CHR)

running_snp <- winScan(x = snp_df,
                       groups = "CHR",
                       position = "KB",
                       values = c("UNAFF", "r"),
                       win_size = 2000,
                       win_step = 500,
                       funs = c("mean"),
                       cores = 8)

ggplot(running_snp , aes(r_mean, UNAFF_mean)) +
        geom_point() +
        xlab("ROH count") +
        ylab("recomb_rate") +
        geom_smooth(method = "lm") +
        facet_wrap(~CHR, scales = "free")
