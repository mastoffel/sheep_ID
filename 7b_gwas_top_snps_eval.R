library(tidyverse)
library(windowscanr)
library(data.table)
source("theme_simple.R")
library(snpStats)
library(viridis)
library(gt)

# prepare data -----------------------------------------------------------------
chr_info <- read_delim("../sheep/data/sheep_genome/chromosome_info_oar31.txt", "\t") %>% 
        .[-1, ] %>% 
        rename(chromosome = Part) %>% 
        mutate(chromosome = str_replace(chromosome, "Chromosome ", "")) %>% 
        mutate(chromosome = as.integer(chromosome)) %>% 
        filter(!is.na(chromosome))

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map <- full_sample$map 
table(full_sample$map$chromosome, useNA = "always")

load("data/survival_mods_data.RData")
# survival data preprocessing
annual_survival <- fitness_data %>% 
        # filter na rows
        filter_at(vars(survival, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
        mutate(age_cent = age - mean(age, na.rm = TRUE),
               age_cent2 = age_cent^2,
               age_std = as.numeric(scale(age)),
               age_std2 = age_std^2,
               # times 10 to estimate a 10% percent increase
               froh_all10 = froh_all * 10,
               froh_all10_cent = froh_all10 - mean(froh_all10, na.rm = TRUE),
               lamb = ifelse(age == 0, 1, 0),
               lamb_cent = lamb - mean(lamb, na.rm = TRUE),
               lamb = as.factor(lamb)) %>% 
        as.data.frame() 

# detailed look at significant regions -----------------------------------------
gwas_res <- read_rds("output/gwas_res_oar_long.rds")
# put into df
gwas_full <- gwas_res %>%
        rename(snp.name = term) %>%
        filter(!str_detect(snp.name, "sd")) %>% 
        mutate(groups = ifelse(str_detect(snp.name, "roh"), "roh", "add")) %>% 
        mutate(snp.name = str_replace(snp.name, "roh_", "")) %>%
        left_join(snps_map) 

# extract roh
gwas_roh <- gwas_full %>% filter(groups == "roh") 

# get cumsums
chr_info2 <- chr_info %>% 
        mutate(tot=cumsum(Length)-Length) %>% 
        dplyr::select(chromosome, tot)

gwas_p <- gwas_roh %>% 
        left_join(chr_info2) %>% 
        # Add a cumulative position of each SNP
        arrange(chromosome, position) %>%
        mutate(positive_cum = position + tot) 

axisdf <- gwas_p %>% group_by(chromosome) %>% 
        summarize(center = (max(positive_cum) + min(positive_cum)) / 2 )

# roh info
hom_sum <- fread("output/ROH/roh.hom.summary") %>% 
        rename(snp.name = SNP, roh_count = UNAFF) %>% 
        dplyr::select(snp.name, roh_count) 

quants <- quantile(hom_sum$roh_count, probs = c(0.1, 0.9))
hom_sum <- hom_sum %>% mutate(roh_prevalence = case_when(
        roh_count < quants[1] ~ "<10% of inds",
        roh_count > quants[2] ~ ">90% of inds",
        TRUE ~ "in between")
)
gwas_plot <- gwas_p %>% left_join(hom_sum)
gwas_plot %>% arrange(p.value)


# plot 1: Regional estimates, log10p and ROH counts ---------------------------- 
get_genome_region <- function(chr, pos) {
        gwas_tbl <- gwas_plot %>% filter(
                (chromosome == chr) & (position >= (pos - plusminus  * 1000000)) & (position <= (pos +  plusminus * 1000000))
        )
        out <- gwas_tbl %>% left_join(hom_sum) %>%
                mutate(log_p = -log10(p.value), estimate = estimate,
                       roh_count = roh_count/5952 * 100) %>%
                pivot_longer(names_to = "var", values_to = "vals",
                             cols = c("log_p", "roh_count", "estimate"))
        
}
plusminus <- 3# Mb

top_snps <- gwas_plot %>% filter(-log10(p.value) > -log10(0.05/39184)) %>%  # -log10(0.05/39149)
        group_by(chromosome) %>% 
        top_n(-1, wt = p.value)

# check MAF of top snps
maf <- col.summary(full_sample$genotypes[, top_snps$snp.name])$maf

# check that GWAS hits are legit using linkage map -----------------------------
lmap <- read_delim("data/Oar3.1_Interpolated.txt", "\t") %>% 
        setNames(c("chr", "snp_name", "bp", "cM")) %>% 
        mutate(mb_pos = bp/1000000,
               kb_pos = bp/1000)

win_area <- function(position, chromosome, ...) {
        diffs <- lmap %>% 
                filter(chr == chromosome) %>% 
                mutate(diff = abs(position/1000 - kb_pos)) %>% 
                arrange(diff) %>% 
                top_n(-1000)
}

all_top <- pmap(top_snps, win_area) %>% 
        map(as_tibble) %>% 
        bind_rows(.id = "top_snps") %>% 
        mutate(top_snps2 = case_when(
                top_snps == 1 ~ "DU289160_204.1, Chr 3, 177.23Mb",
                top_snps == 2 ~ "oar3_OAR10_85579446, Chr 10, 85.59Mb",
                top_snps == 3 ~ "s23340.1, Chr 23, 36.16Mb"
        )) %>% 
        mutate(Mb_pos = kb_pos / 1000)
all_top$top_snp <- rep(top_snps$position, each = 1000)
ggplot(all_top, aes(Mb_pos, cM)) +
        geom_point() +
        theme_simple(grid_lines = FALSE) +
        facet_wrap(~top_snps2, scales = "free") +
        geom_vline(aes(xintercept = top_snp/1e06)) 


# produce table for supplementary ----------------------------------------------
top_snps %>% select(snp.name, chromosome, position, estimate, p.value,   allele.1, allele.2) %>% 
        mutate(estimate = round(estimate, 2)) %>% 
        setNames(c("SNP","Chromosome",  "Position (Bp)", "Estimate (log-odds)", "p-value",  "Allele1", "Allele2")) %>% 
        ungroup() %>% 
        gt() %>% 
        fmt_scientific(
                columns = vars(`p-value`),
                decimals = 2
        ) %>% 
        gtsave("figs/tables/top_hits_gwas.png")

top_snps_region <- map2_df(top_snps$chromosome, top_snps$position, get_genome_region, .id = "snp")

df_plot <- top_snps_region %>% 
        #mutate(top_snp = ifelse(-log10(p.value) >  -log10(0.05/39149), 1, 0)) %>% 
        mutate(top_snp = ifelse(snp.name %in% top_snps$snp.name, 1, 0)) %>% 
        mutate(top_snp = as.factor(top_snp)) %>% 
        mutate(pos_Mb = position/1000000) 

# c("estimate (log-odds)", "-log10(p)", "ROH count")
# "log_p", "roh_count", "estimate"
snp_stats <- c(estimate = "Estimate\n(log-odds)", log_p = "-log10(p-value)", # og_p = "-log10(p-value)", 
               roh_count = "% of sheep\n with ROH")

#levels(df_plot$var) <- c("Estimate~(log-odds)", "-log[10]p-value", "%~sheep~with~ROH")


chrs <- c(`3` = "Chr. 3", `10` = "Chr. 10", `23` = "Chr. 23")
# mean ROH
mean(hom_sum$roh_count)
hlines <- data.frame(y_val = c(NA, NA, 1397/5952 * 100), var = c("estimate", "log_p", "roh_count"))


top_snps_regions_p <- ggplot(df_plot) +
        geom_line(data= df_plot, aes(pos_Mb, y = vals), size =0.05, color = "#4C566A") +
        geom_point(data= df_plot, aes(pos_Mb, y = vals), size =0.2, color = "#D8DEE9") +
        geom_point(data=df_plot %>% filter(top_snp == 1), 
                   aes(pos_Mb, y = vals), shape = 21, fill = viridis(1), size = 1, # viridis(7)[7]
                   stroke = 0.2) +
       geom_hline(data = hlines, aes(yintercept = y_val), size = 0.08,
                   linetype = "dashed") +
        #facet_wrap(var~chromosome, scales = "free", nrow = 3,ncol = 5) +
        facet_grid(var~chromosome, scales = "free", switch = "y",
                   labeller = labeller(.rows = snp_stats,     #snp_stats,
                                        .cols = chrs)) +
                   #labeller = L) +
        scale_x_continuous(breaks = scales::pretty_breaks(3)) +
        scale_y_continuous(position = "right", breaks = scales::pretty_breaks(4)) + # breaks = breaks_fun, 
        # scale_color_manual(values = viridis(2)) +
        #facet_wrap(~var, nrow = 3, scales = "free") + 
        ylab("") +
        xlab("Position in Mb") +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) + 
        theme(axis.line.x=element_line(),
              panel.border = element_rect(size = 0.2, fill = NA))
top_snps_regions_p
#ggsave("figs/top_snps_regions.jpg", plot = top_snps_regions_p,
#       width = 6, height = 3.5)

# get heterozygosity
snps <- unique(df_plot$snp.name)
geno_sub <- as(full_sample$genotypes[, snps ], "numeric")
het <- colSums(geno_sub == 1, na.rm = TRUE) / colSums(!is.na(geno_sub)) %>% 
        as.data.frame() %>% 
        setNames("het") 
het$snp.name <- rownames(het)
het2 <- het %>% 
        left_join(df_plot) 

# smoothing across 30 SNPs
span <- 100
het$heterozygosity <- ksmooth(1:nrow(het), het$het, kernel = "normal", bandwidth = span)$y
hlines <- data.frame(y_val = c(NA, -log10(0.05/39184), 1397/5952 * 100, 0.3083), var = c("estimate", "log_p", "roh_count", "heterozygosity"))

df_plot2 <- df_plot %>% 
                pivot_wider(names_from = var, values_from = vals) %>% 
                left_join(het) %>% 
                pivot_longer(names_to = "var", values_to = "vals",
                        cols = c("log_p", "roh_count", "estimate", "heterozygosity")) %>% 
                mutate(var = factor(var, levels = c("estimate", "log_p", "roh_count", "heterozygosity")))

snp_stats <- c(estimate = "Estimate\n(log-odds)", log_p = "-log10(p-value)", # og_p = "-log10(p-value)", 
               roh_count = "% of sheep\n with ROH", heterozygosity = "SNP\nheterozygosity")
df_plot2 <- df_plot2 %>% 
        mutate(type = ifelse(var %in% c("estimate", "log_p"), "gwas", "diversity")) %>% 
        mutate(type = as.factor(type))

p_final <- ggplot(df_plot2) +
        geom_hline(data = hlines, aes(yintercept = y_val), size = 0.08,
                   linetype = "dashed") +
        geom_line(data= df_plot2, aes(pos_Mb, y = vals, color = type), size =0.05) + #  color = "#4C566A"
        geom_point(data= df_plot2, aes(pos_Mb, y = vals, color = type), alpha = 0.5, size =0.3) + # "#D8DEE9"
        geom_point(data=df_plot2 %>% filter(top_snp == 1), 
                   aes(pos_Mb, y = vals), # fill = type 
                   shape = 21, 
                   #fill = "green",
                   #fill = viridis(2)[2], 
                   #fill = viridis(10)[10],
                   fill = "white",
                   #color = "white",
                   size = 1.5, # viridis(7)[7]
                   stroke = 0.6) +
        #facet_wrap(var~chromosome, scales = "free", nrow = 3,ncol = 5) +
        facet_grid(var~chromosome, scales = "free", switch = "y",
                   labeller = labeller(.rows = snp_stats,     #snp_stats,
                                       .cols = chrs)) +
        #labeller = L) +
        scale_x_continuous(breaks = scales::pretty_breaks(3)) +
        scale_y_continuous(position = "right", breaks = scales::pretty_breaks(4)) + # breaks = breaks_fun, 
        scale_color_manual(values = c("#404788FF", "#238A8DFF")) +
        #scale_fill_manual(values = rev(viridis(2))) +
        # scale_color_manual(values = viridis(2)) +
        #facet_wrap(~var, nrow = 3, scales = "free") + 
        ylab("") +
        xlab("Position in Mb") +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) + 
        theme(axis.line.x=element_line(),
              panel.border = element_rect(size = 0.2, fill = NA),
              legend.position = "none")
p_final
ggsave("figs/Sup_gwas_and_diversity_oar_long.jpg", p_final, width = 6, height = 5)







# plot 2: check whether top SNP roh decreases with increasing age---------------
# this whole analysis doesn't really make sense
ids_per_age <- readRDS("output/ids_per_age.rds")
roh <- fread("output/ROH/roh_ram.hom")

get_roh_per_age <- function(ageclass) {
        ageclass <- ageclass + 1
        
        roh_filt <- roh %>% 
                filter(IID %in% ids_per_age[[ageclass]]) 
        
        #chrs <- unique(top_snps_region$chromosome)
        map_int(1:nrow(top_snps_region), function(x) {
                pos <- as.numeric(top_snps_region[x, "position"])
                chr <- as.numeric(top_snps_region[x, "chromosome"])
                sum((roh_filt$CHR == chr) & (roh_filt$POS1 <= pos) & (roh_filt$POS2 >= pos))
        })
}

top_snps_per_age <- map(c(0:10), get_roh_per_age)

#saveRDS(object = top_snps_per_age, file = "output/roh_at_top_snps_per_age.rds")
top_snps_per_age <- readRDS("output/roh_at_top_snps_per_age.rds")

roh_df <- map_dfc(0:10, function(x) {
        out <- tibble(top_snps_per_age[[x+1]])
        names(out) <- paste0("age",x)
        out
})

plot_per_age <- top_age_roh <- top_snps_region %>% 
        bind_cols(roh_df) %>% 
        pivot_longer(cols = starts_with("age"), names_to = "age", values_to = "roh") %>% 
       # filter(age %in% c("age0", "age3", "age6", "age9")) %>% 
        filter(age %in% c("age0", "age1", "age4", "age8")) %>% 
        mutate(top_snp = ifelse(snp.name %in% top_snps$snp.name, 1, 0)) %>% 
        mutate(top_snp = as.factor(top_snp)) %>% 
        mutate(pos_Mb = position/1000000) 

# ages_labs <- chrs <- c(age0 = "Age 0", age3 = "Age 3", age6 = "Age 6", age9 = "Age 9")
ages_labs <- chrs <- c(age0 = "Age 0", age1 = "Age 1", age4 = "Age 4", age8 = "Age 8")
top_age_roh <- ggplot(plot_per_age, aes(pos_Mb, roh)) +
        geom_point(size =0.15, color = "lightgrey") +
        geom_line(size =0.15, color = "lightgrey") +
        theme_simple(grid_lines = FALSE) +
        facet_grid(age ~ chromosome, scales = "free",
                   labeller = labeller(.rows = ages_labs)) +
        scale_x_continuous(breaks = scales::pretty_breaks(3)) +
        scale_y_continuous(breaks = scales::pretty_breaks(3)) +
        # scale_color_manual(values = viridis(2)) +
        #facet_wrap(~var, nrow = 3, scales = "free") + 
        ylab("") +
        xlab("Position in Mb") +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) + 
        theme(axis.line.x=element_line(),
              panel.border = element_rect(size = 0.2, fill = NA),
              strip.text.x = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 12)) +
        ggtitle("ROH count around GWAS hits in different age classes")

top_age_roh


# make composite plot 1 and 2
top_snps_regions_p2 <- top_snps_regions_p +
        ggtitle("GWAS estimates, p-values and ROH count around the top GWAS hits") +
        theme(plot.title = element_text(hjust = 0.5, size = 12))
sup_plot_gwas_regions <- top_snps_regions_p2 / top_age_roh + plot_annotation(tag_levels = 'A')

#ggsave("figs/age_roh_gwas.jpg", sup_plot_gwas_regions , width = 8, height = 9)

# make composit plot 1 and 2, different style
comp_plot_df <- bind_rows(df_plot, plot_per_age %>% filter(var == "roh_count")) %>% 
                        select(snp, var, vals, roh, pos_Mb, age, chromosome) %>% 
                        mutate(age = ifelse(is.na(age), "all", age)) %>% 
                        mutate(vals = ifelse( (age != "all") & (var == "roh_count"), roh, vals)) %>% 
                        ungroup()

test <- comp_plot_df %>% filter(var == "roh_count")

table(test$var)

ggplot(test, aes(pos_Mb, roh, color = age)) +
        geom_point(alpha = 1, size = 0.1) +
        facet_wrap(var~chromosome, scales = "free")

ggplot(comp_plot_df, aes(pos_Mb, vals, colour = as.factor(age))) +
        geom_point(size =1, alpha = 0.1) +
        #geom_line(size =0.15) +
        theme_simple(grid_lines = FALSE) +
        facet_grid(var ~ chromosome, scales = "free") 
        #facet_grid(rows = vars(var), cols = vars(chromosome), scales = "free") +
       # scale_x_continuous(breaks = scales::pretty_breaks(3)) +
        #scale_y_continuous(breaks = scales::pretty_breaks(3)) +
        # scale_color_manual(values = viridis(2)) +
        #facet_wrap(~var, nrow = 3, scales = "free") + 
        ylab("") +
        xlab("Position in Mb") +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) + 
        theme(axis.line.x=element_line(),
              panel.border = element_rect(size = 0.2, fill = NA),
              strip.text.x = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 12)) +
        ggtitle("ROH count around GWAS hits in different age classes")

comp_plot_df %>% 
        filter(var == "roh_count") %>% 
        ggplot() +
                geom_point(aes(pos_Mb, y = vals, color = as.factor(age)), size =0.15) 

# plot 3: ROH for every hit across all individuals -----------------------------
get_genome_region_roh <- function(x) {
        pos <- as.numeric(top_snps[x, "position"])
        chr <- as.numeric(top_snps[x, "chromosome"])
        gwas_tbl <- roh %>% filter(
                (CHR == chr) & (pos >= POS1) & (pos <= POS2)
        ) 
}
#plusminus <- 2# Mb
roh <- fread("output/ROH/roh.hom")
roh <- as_tibble(roh)
all_roh_regional <- map_dfr(1:nrow(top_snps), get_genome_region_roh, .id = "snp") %>% 
        mutate(POS1 = POS1/1000000,
               POS2 = POS2/1000000)

top_snp_pos <- top_snps %>% 
        ungroup() %>% 
        mutate(snp = c(1:3)) %>% 
        mutate(pos_Mb = position/1000000)

all_roh_regional %>% 
        group_by(snp) %>% 
        mutate(number = 1) %>% 
        mutate(ticker = cumsum(number)) %>% 
        ggplot(aes(y = ticker)) +
        geom_linerange(aes(xmin = POS1, xmax = POS2), alpha = 0.8, size = 0.1 ) +
        theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
        theme(axis.line.y = element_blank()) +
        ylab("Individuals") + 
        facet_wrap(~snp, scales = "free", ncol = 5) +
        geom_vline(data = top_snp_pos, mapping = aes(xintercept = pos_Mb),
                   size = 0.5, color = "#eceff4", linetype="longdash") -> all_roh_regional_p
all_roh_regional_p
ggsave("figs/roh_regional.jpg", width = 10, height = 3.5)



# genes 
gene_files <- list.files("data/ncbi", full.names = TRUE)
all_genes <- map(gene_files, read_csv, skip = 1, col_types = "cddccdc") %>% 
                bind_rows(.id = "chr")

# check maf in regions
top_snps_region




# plot proportion of survival vs. detailed hom/het classes ---------------------
top_snps_lab <- top_snps$snp.name
inds <- colnames(full_sample$genotypes) %in% top_snps_lab
survival <- annual_survival %>% 
                select(id, survival, sheep_year)
geno_info <- as(full_sample$genotypes[, inds], "numeric") %>% 
                as_tibble(rownames = "id") %>% 
                left_join(survival)

rohs <- fread("output/ROH/roh_ram.hom") %>% as_tibble()
rohs

get_roh_per_snp <- function(x) {
        pos <- as.numeric(top_snps[x, "position"])
        chr <- as.numeric(top_snps[x, "chromosome"])
        gwas_tbl <- rohs %>% 
                mutate(snp = (CHR == chr) & (pos >= POS1) & (pos <= POS2)) %>% 
                group_by(IID) %>% 
                summarise(snp = sum(snp)) %>% 
                rename(id = "IID")
        #names(gwas_tbl)[names(gwas_tbl) == "snp"] <- paste0("snp", x)
        #as.numeric(gwas_tbl[[paste0("snp", x)]])
}
#plusminus <- 2# Mb
rohs_per_top_snp <- map(1:nrow(top_snps), get_roh_per_snp) %>% 
                        map2(top_snps$snp.name, function(x, y) {
                                names(x)[2] <- y
                                x
                                }) %>% 
                        reduce(left_join) %>% 
                        mutate(id = as.character(id))
names(rohs_per_top_snp)[2:5] <- paste0("roh_", names(rohs_per_top_snp)[2:5])

# add to additive genotypes
geno_full <- geno_info %>% 
        left_join(rohs_per_top_snp) %>% 
        select(id, survival, sheep_year, everything()) 
geno_simple <- setNames(geno_full, 
                        c("id", "survival", "sheep_year", 
                          paste0("snp", 1:4), paste0("roh_snp", 1:4)))


# check additive effects
gwas_full %>% filter(snp.name %in% top_snps$snp.name)

# check MAF
map(geno_simple[4:7], function(x) {
        freqs <- as.numeric(table(x))
        maf <- (freqs[1]*2 + freqs[1]) / (sum(freqs) * 2)
})


geno_simple <- geno_simple %>%
        mutate(class1 = case_when(
                snp1 == 0 & roh_snp1 == 0 ~ "0",
                snp1 == 0 & roh_snp1 == 1 ~ "0roh",
                snp1 == 2 & roh_snp1 == 0 ~ "2",
                snp1 == 2 & roh_snp1 == 1 ~ "2roh",
                snp1 == 1 ~ "1",
                is.na(snp1) | is.na(roh_snp1) ~ "NA"
        )) %>% 
        mutate(class2 = case_when(
                snp2 == 0 & roh_snp2 == 0 ~ "0",
                snp2 == 0 & roh_snp2 == 1 ~ "0roh",
                snp2 == 2 & roh_snp2 == 0 ~ "2",
                snp2 == 2 & roh_snp2 == 1 ~ "2roh",
                snp2 == 1 ~ "1",
                is.na(snp2) | is.na(roh_snp2) ~ "NA"
        )) %>% 
        mutate(class3 = case_when(
                snp3 == 0 & roh_snp3 == 0 ~ "0",
                snp3 == 0 & roh_snp3 == 1 ~ "0roh",
                snp3 == 2 & roh_snp3 == 0 ~ "2",
                snp3 == 2 & roh_snp3 == 1 ~ "2roh",
                snp3 == 1 ~ "1",
                is.na(snp3) | is.na(roh_snp3) ~ "NA"
        )) %>% 
        mutate(class4 = case_when(
                snp4 == 0 & roh_snp4 == 0 ~ "0",
                snp4 == 0 & roh_snp4 == 1 ~ "0roh",
                snp4 == 2 & roh_snp4 == 0 ~ "2",
                snp4 == 2 & roh_snp4 == 1 ~ "2roh",
                snp4 == 1 ~ "1",
                is.na(snp4) | is.na(roh_snp4) ~ "NA"
        )) 

geno_simple %>% group_by(class4) %>% summarise(mean(survival), n())
geno_simple %>% group_by(id) %>% sample_n(1) %>% group_by(class1) %>% summarise(mean(survival), n())

ggplot(test, aes(class1)) +
        geom_bar(aes(fill = survival))

