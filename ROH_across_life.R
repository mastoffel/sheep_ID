library(dplyr)
library(tidyverse)
library(ggridges)
# ROH across the lifetime
load("model_in/fitness_roh_df.RData")

fitness_data %>% 
        group_by(age) 

fitness_data %>% 
        filter(Age<11) %>% 
        mutate(Age = as.factor(Age)) %>% 
        pivot_longer(cols = starts_with("FROH"),
                     names_to = "roh_class",
                     values_to = "prop_roh") %>% 
        ggplot(aes(prop_roh)) +
        geom_histogram(bins = 200) +
        #facet_wrap(Age ~ roh_class,  ncol = 4) +
        facet_grid(Age ~ roh_class, scales = "free_x") +
        scale_y_log10()


fitness_data %>% 
        filter(Age<11) %>% 
        mutate(Age = as.factor(Age)) %>% 
        group_by(Age) %>% 
        summarise_at(vars(FROH_short:FROH_long), funs(median, mean))

fitness_data %>% 
        filter(Age == 0 | Age == 1) %>% 
        filter(!is.na(Survival)) %>% 
        mutate(Survival = as.factor(Survival)) %>% 
        ggplot(aes(Survival, FROH_long)) +
        geom_boxplot() 



fitness_data %>% 
        filter(Age<11) %>% 
        mutate(Age = as.factor(Age)) %>% 
        ggplot(aes(FROH_long, Age)) +
        stat_binline() 


fitness_data %>% 
        group_by(Age) %>% 
        tally()

fitness_data %>% 
        group_by(Age) %>% 
        summarise_at(c("FROH_long", "FROH_medium", "FROH_short", "hom"), max, na.rm = TRUE)

fitness_data %>% 
        rename(FROH_not_roh = hom) %>%
        dplyr::select(-FROH_all) %>% 
        filter(FROH_not_roh < 0.55) %>% 
        dplyr::select(starts_with("FROH"), Age, SEX, ID) %>% 
        filter(SEX %in% c("M", "F")) %>% 
        pivot_longer(
                cols = starts_with("FROH"),
                names_to = "ROH",
                values_to = "Prop",
                values_drop_na = TRUE
        ) %>% 
        ggplot(aes(as.factor(Age), Prop, color = as.factor(ROH))) + 
        geom_point(size = 2, alpha = 0.3) + 
        geom_smooth()  +
        facet_wrap(SEX~ROH) + 
        theme_minimal()

fitness_data %>% 
        rename(FROH_not_roh = hom) %>% 
        dplyr::select(FROH_all, Age, SEX, ID) %>% 
        #dplyr::select(starts_with("FROH"), Age, SEX, ID) %>% 
        filter(SEX %in% c("M", "F")) %>% 
        pivot_longer(
                cols = starts_with("FROH"),
                names_to = "ROH",
                values_to = "Prop",
                values_drop_na = TRUE
        ) %>% 
        ggplot(aes(Age, Prop, color = SEX)) + geom_smooth(method = "lm") + 
                geom_point() 
        



#~~~ ROH density
hom_sum <- fread("output/ROH/roh_nofilt.hom.summary")

hom_sum <- hom_sum %>%
        mutate(MB = BP / 100000,
               index = 1:nrow(.))


# Here's how you can do this with the dplyr functions group_by and do:
window_width <- 500
jumps <- 500
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
        filter(CHR %in% 1) %>% 
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
