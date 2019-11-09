library(tidyverse)
library("corpcor")
library(snpStats)
library(tidyimpute)

chr <- 20

# general data
# data
load("data/fitness_roh_df.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing.txt", delim = " ")

# genotypes
# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019_pruned"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map_sub <- full_sample$map %>%  filter(chromosome == chr) 
# prepare additive genotypes subset
snps_sub <- full_sample$map %>% filter(chromosome == chr) %>%  .$snp.name
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_sub], Class = "numeric"), rownames = "id")

# fitness
# survival data
annual_survival <- fitness_data %>% 
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
                      survival = Survival) %>% 
        # some individuals arent imputed well and should be discarded 
        filter(!(id %in% IDs_lots_missing$id)) %>% 
        #filter(age == 0) %>% 
        filter(!is.na(survival)) %>% 
        filter(!is.na(froh_all)) %>% 
        filter(!(is.na(birth_year) | is.na(sheep_year))) %>%  # no mum_id here
        mutate_at(c("id", "birth_year", "sex", "sheep_year", "survival"), as.factor) %>% 
        mutate(age2 = age^2) %>% 
        mutate(age_std = as.numeric(scale(age)),
               age2_std = as.numeric(scale(age2))) %>% 
        as.data.frame() 


# testset
annual_survival_sub <- annual_survival %>% 
                         sample_frac(0.05) 
# genos
geno_mat <- annual_survival_sub %>% 
                select(id) %>% 
                left_join(geno_sub) %>% 
                select(-id) %>% 
                impute_mean() %>% 
                # filter snps without variation
                .[(map_lgl(., function(x) length(table(x)) != 1))] %>% 
                # standardise genos
               # mutate_all(function(x) (x-mean(x, na.rm = TRUE))/ sd(x, na.rm = TRUE)) %>% 
                as.matrix()

rownames(geno_mat) <- annual_survival_sub$id

annual_survival_sub <- annual_survival_sub %>% mutate(survival = as.numeric(as.character(survival)))
library(bWGR)
mod <- mixed(y = survival, random = ~birth_year + sheep_year + id, 
             fixed = ~sex + twin, data = annual_survival_sub,
             X=list(id = geno_mat), maxit = 100)
plot(mod$Structure$id)


svd_geno <- svd(geno_mat)
U <- svd_geno$u # dim N * N -> eigenvectors of GRM
S <- diag(svd_geno$d) # dim N * N singular values matrix
V <- svd_geno$v # dim k * N -> describes LD structure

effs <- V %*% solve(S) %in% t(U) %*% y

sigma2a <- 0.01
sigma2e <- 0.6
lambda <- sigma2e/sigma2a

# blup of marker effects
sigma2a <- 0.05
sigma2e <- 0.4
lambda <- sigma2e/sigma2a
y <- as.numeric(as.character(annual_survival_sub$survival))
W <- geno_mat
Wt <- t(W)
WWt <- W%*%Wt
Ivar <- lambda*(diag(nrow(W)))
WWt2 <- WWt + Ivar
WWt2_inv <- solve(WWt2)
blup_a <- (Wt %*% WWt2_inv) %*% y
plot(blup_a ^ 2)

blup_a <- (Wt %*% solve(W%*%Wt + lambda*(diag(nrow(W))))) %*% y


# eq 2
# inverse
shat_0 <- solve(S%*%S + (diag(nrow(S)) * lambda))
shat <- shat_0 * sigma2e
b_markers <- V %*% shat
b_markers %>% rowMeans() %>% plot()

V %*% solve(S)

y <- as.numeric(as.character(annual_survival_sub$survival))
X <- model.matrix(~annual_survival_sub$twin)
        
WtW <- geno_mat %*% t(geno_mat)
varcompSNP <- regress(y ~ X, ~WtW)

sigma2a <- varcompSNP$sigma[1]
sigma2e <- varcompSNP$sigma[2]  # residual variance -> sigma^2_e

lambda <- sigma2e/sigma2a


