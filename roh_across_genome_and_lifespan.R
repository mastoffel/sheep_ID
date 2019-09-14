library(tidyverse)
library(data.table)
library(magrittr)
library(rlang)
source("theme_clean.R")
# annual measures of traits and fitness
load("model_in/fitness_roh_df.RData")
library(windowscanr)

# number of individuals in each age class
max_age_df <- fitness_data %>% 
        group_by(ID) %>% 
        mutate(max_age = max(Age)) %>% 
        filter(Age == max_age) %>% 
        filter(BIRTHYEAR != 2018) %>% 
       # filter(is.na(Survival)) %>% 
        mutate(age_class_0 = ifelse(Age >= 0, 1, 0),
               age_class_1 = ifelse(Age >= 1, 1, 0),
               age_class_2 = ifelse(Age >= 2, 1, 0),
               age_class_3 = ifelse(Age >= 3, 1, 0),
               age_class_4 = ifelse(Age >= 4, 1, 0),
               age_class_5 = ifelse(Age >= 5, 1, 0),
               age_class_6 = ifelse(Age >= 6, 1, 0),
               age_class_7 = ifelse(Age >= 7, 1, 0),
               age_class_8 = ifelse(Age >= 8, 1, 0),
               age_class_9 = ifelse(Age >= 9, 1, 0),
             age_class_10 = ifelse(Age >= 10, 1, 0)
        ) %>% 
        ungroup() %>% 
        filter(!is.na(FROH_all)) 
        

# number of individuals with genotype data per age class
# num_ind_per_age <- max_age_df %>% 
#         summarise_at(vars(age_class_0:age_class_10), sum)

ids_per_age <- map(0:10, function(x) {
        to_filt <- paste0("age_class_", x)
        max_age_df %>% 
                filter(get(to_filt) == 1) %>% 
                dplyr::select(ID) %>% 
                unlist()
})
names(ids_per_age) <- paste0("age_", c(0:10))

num_ind_per_age <- unlist(map(ids_per_age, length))

# load ROH info
file_path <- "output/ROH/roh_nofilt_1Mb.hom"
file <- "roh_nofilt"
roh_lengths <- fread(file_path)

# this is needed for sliding window analyses
hom_sum <- fread("output/ROH/roh_nofilt_1Mb.hom.summary") 

# define vectorized seq to work with mutate
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))


# 6908
# work on subsets for now 
hom_sum_sub <- hom_sum %>% 
        filter(CHR %in% c(19,20))

roh_lengths_sub <- roh_lengths %>% 
        filter(CHR %in% c(19,20))


# function to calculate how many individuals have a given SNP in an ROH
roh_per_snps <- function(roh_lengths, hom_sum, ids = NULL, subsamp = NULL) {
        
        if (is.null(ids)) ids <- unique(roh_lengths$FID)
        if (is.null(subsamp)) subsamp <- length(ids)
        
        # create indices for all rohs
        roh_lengths <- roh_lengths %>% 
                as_tibble() %>% 
                filter(FID %in% ids) %>% 
                filter(FID %in% sample(ids, subsamp)) %>% 
                mutate(index1 = match(SNP1, hom_sum$SNP),
                       index2 = index1 + NSNP - 1) %>% 
                mutate(all_ind = seq2(from = index1, to = index2))
        num_roh_per_snp <- tabulate(unlist(roh_lengths$all_ind), nbins = nrow(hom_sum))
}

first_age <- "age_1"
second_age <- "age_2"

roh_first_age <- roh_per_snps(roh_lengths_sub, hom_sum_sub, ids = ids_per_age[[first_age]]) / length(ids_per_age[[first_age]])
roh_second_age <- roh_per_snps(roh_lengths_sub, hom_sum_sub, ids = ids_per_age[[second_age]]) / length(ids_per_age[[second_age]])
roh_subsamples <- 1000 %>% rerun(roh_per_snps(roh_lengths_sub, hom_sum_sub, 
                                                  # subsample only ids which were present at first age class
                                                  ids = ids_per_age[[first_age]], 
                                                  # subsample as many as survived into the second age class
                                                  subsamp = num_ind_per_age[[second_age]])) %>% 
                        bind_cols() %>% 
                        mutate_all(function(x) x/length(ids_per_age[[second_age]]))

all_quants <- apply(roh_subsamples, 1, function(x) as.data.frame(t(quantile(x, probs = c(0.005, 0.995))))) %>% 
                bind_rows() %>% 
                rename(lower_ci = `0.5%`,
                       upper_ci = `99.5%`)
                # rename(lower_ci = `2.5%`,
                #       upper_ci = `97.5%`)


roh_age_0_1 <- hom_sum_sub %>% 
        bind_cols(tibble(age_0 = roh_first_age, age_1 = roh_second_age), all_quants) %>% 
        mutate(KB = BP/1000)

   

running_roh <- winScan(x = roh_age_0_1, 
                   groups = NULL, 
                   position = "KB", 
                   values = c("age_0", "age_1", "lower_ci", "upper_ci"), 
                   win_size = 500,
                   win_step = 100,
                   funs = c("mean"))

# mark points where ROH is lower than CI_lower 
running_roh %<>% 
        mutate(outliers = ifelse(age_1_mean < lower_ci_mean, 1, 0)) %>% 
        mutate(outliers = as.factor(outliers))

running_roh %>% 
        ggplot() +
        geom_line(aes(win_start, age_0_mean)) +
        geom_line(aes(win_start, lower_ci_mean), color = "lightgrey") +
        geom_line(aes(win_start, upper_ci_mean), color = "lightgrey") +
        geom_line(aes(win_start, age_1_mean), color = "cornflowerblue") +
        geom_point(aes(win_start, age_1_mean, color = outliers), size = 0.1) + 
        scale_color_manual(values = c('cornflowerblue', "red")) +
        theme_clean() +
        ylab("Proportion of Population with ROH") +
        xlab("Bp")

running_roh %>% filter(outliers == 1)

p_running_roh <- running_roh %>% 
        filter(CHR %in% 1:26) %>% 
        ggplot(aes(window.start, ninds)) +
        geom_line() +
        scale_y_continuous(breaks = c(369, 1844, 5531, 7005), labels = c("5%", "25%", "75%", "95%"),
                           limits = c(0,7374)) +
        theme_clean() + 
        ggtitle("Average ROH across the genome of 7374 Soay sheep") +
        xlab("Window start (in Base Pairs) of 200 BP window") +
        ylab("Individuals with ROH") + 
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=1844, ymax=5531, alpha=0.1) +
        geom_hline(yintercept = 369, size = 0.3) +
        geom_hline(yintercept = 7005, size = 0.3) +
        facet_grid(CHR~.) 
p_running_roh 



# add to df
# hom_sum_sub$roh <- num_roh_per_snp







# check
any(!hom_sum_sub$UNAFF == hom_sum_sub$roh)











# create vector with all SNP names
# snps <- hom_sum %>% 
#         filter(CHR == 26) %>% 
#         .[["SNP"]]


test <- all_roh %>% 
        as_tibble() %>% 
        mutate(index1 = match(SNP1, snps),
               index2 = index1 + NSNP - 1) %>% 
        mutate(all_ind = seq2(from = index1, to = index2))

num_roh_per_snp <- table(unlist(test$all_ind))
num_roh_per_snp <- tabulate(unlist(test$all_ind))
        
hom_sum$roh <- 0
hom_sum[as.numeric(names(num_roh_per_snp)), "roh"] <- as.numeric(num_roh_per_snp)


Vectorize(seq.default, vectorize.args = c("from", "to"))

eval("553:951")

seq()






hom_sum <- hom_sum %>%
        mutate(MB = BP / 100000,
               index = 1:nrow(.))


# Here's how you can do this with the dplyr functions group_by and do:
window_width <- 200
jumps <- 200
running_roh <- hom_sum %>% 
        group_by(CHR) %>% 
        do(
                data.frame(
                        window.start = rollapply(.$BP, width= window_width, by=jumps, FUN=min, align="left"),
                        window.end = rollapply(.$BP, width=window_width, by=jumps, FUN=max, align="left"),
                        ninds = rollapply(.$UNAFF, width=window_width, by=jumps, FUN=mean, align="left")
                )
        )

# running_roh$index = 1:nrow(x)
p_running_roh <- running_roh %>% 
        filter(CHR %in% 1:26) %>% 
        ggplot(aes(window.start, ninds)) +
        geom_line() +
        scale_y_continuous(breaks = c(369, 1844, 5531, 7005), labels = c("5%", "25%", "75%", "95%"),
                           limits = c(0,7374)) +
        theme_clean() + 
        ggtitle("Average ROH across the genome of 7374 Soay sheep") +
        xlab("Window start (in Base Pairs) of 200 BP window") +
        ylab("Individuals with ROH") + 
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=1844, ymax=5531, alpha=0.1) +
        geom_hline(yintercept = 369, size = 0.3) +
        geom_hline(yintercept = 7005, size = 0.3) +
        facet_grid(CHR~.) 
p_running_roh 

