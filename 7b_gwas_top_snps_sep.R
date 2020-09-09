# This script is to plot genetic diversity and model estimates
# in regions around GWAS peaks. 
# Contains some exploratory analyses which are not in the manuscript.

library(tidyverse)
library(windowscanr)
library(data.table)
source("theme_simple.R")
library(snpStats)
library(viridis)
library(gt)

# prepare data -----------------------------------------------------------------
chr_info <- read_delim("data/chromosome_info_oar31.txt", "\t") %>% 
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

# fitness data
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
gwas_res <- read_rds("output/gwas_res_oar_roh_sep.rds")
#add map positions
gwas_out <- gwas_res %>% 
                rename(snp.name = snp) %>% 
                left_join(snps_map) %>% 
                rename(snp = snp.name)

# plot 1: Regional estimates, log10p and ROH counts ---------------------------- 

# filter top snp from every peak
top_snps <- gwas_out %>% 
        filter(state != "add") %>% 
        filter(-log10(p.value) > -log10(0.05/(2*39184))) %>%  # -log10(0.05/39149)
        mutate(peak = ifelse((chromosome == 3) & (position > 13000000 & position < 14000000), "3a", chromosome)) %>% 
        mutate(peak = ifelse((chromosome == 3) & (position > 170000000 & position < 180000000), "3b", chromosome)) %>% 
        group_by(peak) %>% 
       # mutate(pos_diff = position - lag(position, n = 1, default = position[1]))
        top_n(-1, wt = p.value) %>% 
        # two snps have same pvalue
        filter(!(snp %in% c("oar3_OAR10_85722487", "oar3_OAR3_177455437")))

# check MAF of top snps
maf <- col.summary(full_sample$genotypes[, top_snps$snp])$MAF



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
                top_n(-100)
}

all_top <- pmap(top_snps, win_area) %>% 
        map(as_tibble) %>% 
        bind_rows(.id = "top_snps") %>% 
        # mutate(top_snps2 = case_when(
        #         top_snps == 1 ~ "DU289160_204.1, Chr 3, 177.23Mb",
        #         top_snps == 2 ~ "oar3_OAR10_85579446, Chr 10, 85.59Mb",
        #         top_snps == 3 ~ "s23340.1, Chr 23, 36.16Mb"
        # )) %>% 
        mutate(Mb_pos = kb_pos / 1000)

all_top$top_snp <- rep(top_snps$position, each = 100)

p_top <- ggplot(all_top, aes(Mb_pos, cM)) +
        geom_point(shape = 21, size = 2, fill = "#eceff4", stroke = 0.1) +
        theme_simple(grid_lines = FALSE) +
        facet_wrap(~top_snps, scales = "free") +
        geom_vline(aes(xintercept = top_snp/1e06)) 

# produce table for supplementary ----------------------------------------------
top_snps %>% 
        ungroup() %>% 
        dplyr::select(snp, chromosome, position, state, estimate, p.value,  allele.1, allele.2) %>% 
        mutate(estimate = round(estimate, 2)) %>% 
        mutate(state = case_when(
                state == "roh_2" ~ paste0(allele.1, allele.1),
                state == "roh_0" ~ paste0(allele.2, allele.2))) %>% 
        arrange(chromosome) %>% 
        setNames(c("SNP","Chromosome",  "Position (Bp)", "ROH status", "Estimate (log-odds)", "p-value", "Allele1", "Allele2")) %>% 
        ungroup() %>% 
        gt() %>% 
        cols_align(
                align = "center") %>% 
        fmt_scientific(
                columns = vars(`p-value`),
                decimals = 2) %>% 
        gtsave("figs/tables/top_hits_gwas.png")

hom_sum <- fread("output/ROH/roh.hom.summary") %>% 
        rename(snp = SNP, roh_count = UNAFF) %>% 
        dplyr::select(snp, roh_count) 

get_genome_region <- function(chr, pos) {
        gwas_tbl <- gwas_out %>% filter(
                (chromosome == chr) & (position >= (pos - plusminus  * 1000000)) & (position <= (pos +  plusminus * 1000000))
        )
        out <- gwas_tbl %>% left_join(hom_sum) %>%
                mutate(direction = case_when(
                        estimate < 0 ~ "negative",
                        estimate >= 0 ~ "positive"
                )) %>% 
                mutate(direction = ifelse(p.value > (0.05/(2*39184)), "non_sig", direction)) %>% 
             #   mutate(direction = ifelse((estimate < 0 & p.value < 0.05/(2*39184)), "negative", "positive")) %>% 
                mutate(log_p = -log10(p.value), estimate = abs(estimate),                #abs(estimate),
                       roh_count = roh_count/5952 * 100) %>%
                #filter(p.value < 0.05) %>% 
                pivot_longer(names_to = "var", values_to = "vals",
                             cols = c("log_p", "roh_count", "estimate"))
        
}
plusminus <- 1# Mb
top_snps <- arrange(top_snps, chromosome)
top_snps_region <- map2_df(top_snps$chromosome, top_snps$position, get_genome_region, .id = "peak")

df_plot <- top_snps_region %>% 
        #mutate(top_snp = ifelse(-log10(p.value) >  -log10(0.05/39149), 1, 0)) %>% 
        mutate(top_snp = ifelse(snp %in% top_snps$snp, 1, 0)) %>% 
        mutate(top_snp = as.factor(top_snp)) %>% 
        mutate(pos_Mb = position/1000000) %>% 
        group_by(snp) %>%
        top_n(-1, p.value)

# c("estimate (log-odds)", "-log10(p)", "ROH count")
# "log_p", "roh_count", "estimate"
snp_stats <- c(estimate = "Estimate\n(log-odds)", log_p = "-log10(p-value)", # og_p = "-log10(p-value)", 
               roh_count = "% of sheep\n with ROH")

#levels(df_plot$var) <- c("Estimate~(log-odds)", "-log[10]p-value", "%~sheep~with~ROH")

chrs <- c(`1` = "Chr. 3", `2` = "Chr. 3", `3` = "Chr. 10", `4` = "Chr. 14", `5` = "Chr. 18", `6` = "Chr. 19",
          `7` = "Chr. 23")
# mean ROH
mean(hom_sum$roh_count)
hlines <- data.frame(y_val = c(NA, NA, 1397/5952 * 100), var = c("estimate", "log_p", "roh_count"))


# top_snps_regions_p <- ggplot(df_plot) +
#         geom_line(data= df_plot, aes(pos_Mb, y = vals), size =0.05, color = "#4C566A") +
#         geom_point(data= df_plot, aes(pos_Mb, y = vals), size =0.2, color = "#D8DEE9") +
#         geom_point(data=df_plot %>% filter(top_snp == 1), 
#                    aes(pos_Mb, y = vals), shape = 21, fill = viridis(1), size = 1, # viridis(7)[7]
#                    stroke = 0.2) +
#         geom_hline(data = hlines, aes(yintercept = y_val), size = 0.08,
#                    linetype = "dashed") +
#         #facet_wrap(var~chromosome, scales = "free", nrow = 3,ncol = 5) +
#         facet_grid(var~peak, scales = "free", switch = "y",
#                    labeller = labeller(.rows = snp_stats,     #snp_stats,
#                                        .cols = chrs)) +
#         #labeller = L) +
#         scale_x_continuous(breaks = scales::pretty_breaks(3)) +
#         scale_y_continuous(position = "right", breaks = scales::pretty_breaks(4)) + # breaks = breaks_fun, 
#         # scale_color_manual(values = viridis(2)) +
#         #facet_wrap(~var, nrow = 3, scales = "free") + 
#         ylab("") +
#         xlab("Position in Mb") +
#         theme_simple(grid_lines = FALSE, axis_lines = TRUE) + 
#         theme(axis.line.x=element_line(),
#               panel.border = element_rect(size = 0.2, fill = NA))
# top_snps_regions_p
#ggsave("figs/top_snps_regions.jpg", plot = top_snps_regions_p,
#       width = 6, height = 3.5)


# heterozygosity
snps <- unique(df_plot$snp)
geno_sub <- as(full_sample$genotypes[, snps ], "numeric")
het <- colSums(geno_sub == 1, na.rm = TRUE) / colSums(!is.na(geno_sub)) %>% 
        as.data.frame() %>% 
        setNames("het") 
het$snp <- rownames(het)
het2 <- het %>% 
        left_join(df_plot) 

# smoothing across 30 SNPs
span <- 100
het$heterozygosity <- ksmooth(1:nrow(het), het$het, kernel = "normal", bandwidth = span)$y
hlines <- data.frame(y_val = c(NA, -log10(0.05/(2*39184)), 1397/5952 * 100, 0.3083), var = c("estimate", "log_p", "roh_count", "heterozygosity"))

df_plot2 <- df_plot %>% 
        pivot_wider(names_from = var, values_from = vals) %>% 
        left_join(het) %>% 
        pivot_longer(names_to = "var", values_to = "vals",
                     cols = c("log_p", "roh_count", "estimate", "heterozygosity")) %>% 
        mutate(var = factor(var, levels = c("estimate", "log_p", "roh_count", "heterozygosity")))

snp_stats <- c(estimate = "|Estimate|\n(log-odds)", log_p = "-log10(p-value)", # og_p = "-log10(p-value)", 
               roh_count = "% of sheep\n with ROH", heterozygosity = "SNP\nheterozygosity")
df_plot2 <- df_plot2 %>% 
        mutate(type = ifelse(var %in% c("estimate", "log_p"), "gwas", "diversity")) %>% 
        mutate(type = as.factor(type))

df_plot2 <- df_plot2 %>% 
        mutate(direction = ifelse(type == "diversity", "non_sig", direction)) 

cols <- c(viridis(2)[1], "#d8dee9", viridis(2)[2])
p_final <- ggplot(df_plot2) +
        geom_hline(data = hlines, aes(yintercept = y_val), size = 0.08,
            linetype = "dashed") +
        geom_line(data= df_plot2, aes(pos_Mb, y = vals, color = direction), size =0.05) + #  color = "#4C566A"
        geom_point(data= df_plot2, aes(pos_Mb, y = vals, color = direction, size = direction), alpha = 0.5) + # "#D8DEE9"
        # geom_point(data=df_plot2 %>% filter(top_snp == 1), 
        #            aes(pos_Mb, y = vals), # fill = type 
        #            shape = 21, 
        #            #fill = "green",
        #            #fill = viridis(2)[2], 
        #            #fill = viridis(10)[10],
        #            fill = "white",
        #            #color = "white",
        #            size = 1, # viridis(7)[7]
        #            stroke = 0.6) +
        #facet_wrap(var~chromosome, scales = "free", nrow = 3,ncol = 5) +
        facet_grid(var~peak, scales = "free", switch = "y",
                   labeller = labeller(.rows = snp_stats,     #snp_stats,
                                       .cols = chrs)) +
        #facet_wrap(var~peak, scales = "free") +
        #labeller = L) +
        scale_x_continuous(breaks = scales::pretty_breaks(2)) +
        scale_y_continuous(position = "right", breaks = scales::pretty_breaks(4)) + # breaks = breaks_fun, 
        #scale_fill_manual(values = cols) + #c("#404788FF", "#238A8DFF")
        scale_color_manual(values = cols) +
        scale_size_manual(values = c(1.5, 0.1, 1.5)) +
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
ggsave("figs/Sup_gwas_and_diversity_oar_long.jpg", p_final, width = 8, height = 5)


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
        mutate(snp = c(1:7)) %>% 
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
#ggsave("figs/roh_regional.jpg", width = 10, height = 3.5)

all_roh_regional %>% 
        group_by(snp) %>% 
        summarise(mean(KB))
mean(roh$KB)

all_roh_regional %>% 
        group_by(snp) %>% 
        tally() %>% 
        mutate(prop_roh = n/5952)

