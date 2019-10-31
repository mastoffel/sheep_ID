# evaluate whether additive genetic variation changes things
library(tidyverse)
library(data.table)
library(snpStats)
library(rrBLUP)
library(tidyimpute)
library(regress)
library(asreml)
library(snpReady)
library(MASS)

# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_27092019_pruned"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_chr20 <- full_sample$map %>% filter(chromosome == 20) %>% .$snp.name
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_chr20], Class = "numeric"),
                      rownames = "id")

# surv dat
sheep_dat <- read_delim("output/early_survival_top_snps_pca.txt", delim = " ") %>% 
                sample_frac(1) %>% 
                filter(!is.na(sex) & !is.na(twin)) %>% 
                mutate(id = as.factor(id)) %>% 
                mutate_at(vars(sex:mum_id), as.factor) 

sheep_dat2 <- sheep_dat[, 1:9] %>% 
        left_join(geno_sub, by = "id")


genos_pre <- sheep_dat2 %>% 
         dplyr::select(id, oar3_OAR20_121503:oar3_OAR20_51006360) %>% 
         group_by(id) %>% 
         #top_n(n = 1)
         summarise_all(first)

not_zero_var <- map_lgl(genos_pre, function(x) !(length(unique(x)) == 1))
genos_pre <- genos_pre[not_zero_var]

ids <- genos_pre$id
genos <- genos_pre %>% 
         dplyr::select(-id) %>%
         impute_mean() %>% 
         as.matrix()
grm <- G.matrix(genos, method = "VanRaden", format = "wide")$Ga
ginv <- ginv(grm)

# prepare for asreml
rownames(ginv) <- ids
colnames(ginv) <- ids
attr(ginv, "INVERSE") <- TRUE
head(ginv)

modelGBLUP <- asreml(fixed = survival ~ sex + twin + 1,
                     random = ~vm(id, ginv),
                     workspace = 128e06, family = asr_binomial(),
                     data = sheep_dat)

plot(modelGBLUP)
summary(modelGBLUP)$varcomp
out <- summary(modelGBLUP, coef=TRUE)
out$coef.fixed
BLUP <- summary(modelGBLUP, coef=TRUE)$coef.random
cor(sheep_dat$survival, BLUP, method = "pearson", use = "complete.obs")

# calculate SNP effects
# Incidence matrix X, I, and phenotype y
X <- model.matrix(~ sheep_dat2$sex + sheep_dat2$twin)
W <- sheep_dat2 %>% 
        dplyr::select(oar3_OAR20_121503:oar3_OAR20_51006360) %>% 
        impute_mean() %>% 
        as.matrix() 
# n*n diagonal mat I
Inxn <- diag(nrow(W))
y <- sheep_dat2$survival
p <- colMeans(W)/2
# Quality control of genotype matrix / scale geno mat
W2 <- scale(W, center = TRUE, scale = FALSE)
# m*m diagonal matrix
Imxm <- diag(ncol(W2))
#
WtW <- W2 %*% t(W2)
# varcomps
sigma2a <- summary(modelGBLUP)$varcomp$component[1]  # additive SNP variance -> sigma^2_a
sigma2e <- summary(modelGBLUP)$varcomp$component[2]   # residual variance -> sigma^2_e
lambda <- sigma2e/sigma2a
# inverse of V
V <- W2 %*% Imxm %*% t(W2) * sigma2a + Inxn * sigma2e
dim(V)
Vinv <- solve(V)
dim(Vinv)

# snp effects
fix_eff_est <- out$coef.fixed[c(5,4,2), 1]
W2 <- as.matrix(W2)
ahat2 <- sigma2a * Imxm %*% t(W2) %*% Vinv %*% (matrix(y) - matrix(X %*% fix_eff_est)) 

plot(ahat2)







# Incidence matrix X, I, and phenotype y
X <- model.matrix(~ sheep_dat$sex + sheep_dat$twin)
W <- sheep_dat %>% dplyr::select(oar3_OAR15_63743364:oar3_OAR15_64194091) %>% 
        impute_mean()

# n*n diagonal mat I
Inxn <- diag(nrow(W))
y <- sheep_dat$survival

p <- colMeans(W)/2

# Quality control of genotype matrix / scale geno mat
W2 <- scale(W)
# m*m diagonal matrix
Imxm <- diag(ncol(W2))

#Variance components estimation REML
library(regress)
WtW <- W2 %*% t(W2)
varcompSNP <- regress(y ~ -1 + X, ~WtW)

sigma2a <- varcompSNP$sigma[1]  # additive SNP variance -> sigma^2_a
sigma2e <- varcompSNP$sigma[2]  # residual variance -> sigma^2_e

lambda <- sigma2e/sigma2a

# Inverse of V
V <- W2 %*% Imxm %*% t(W2) * sigma2a + Inxn * sigma2e
dim(V)
Vinv <- solve(V)
dim(Vinv)

# Ordinary least squares
fit <- lm(survival ~ twin + sex, data = sheep_dat)
summary(fit)
residuals(fit)

# snp effects
# ahat <- sigma2a * Imxm %*% t(W2) %*% Vinv %*% (matrix(y) - matrix(predict(fit)))
ahat2 <- sigma2a * Imxm %*% t(W2) %*% Vinv %*% (matrix(y) - matrix(X %*% fit$coefficients))  # alternative
table(ahat == ahat2)
head(ahat)
tail(ahat)

library(rrBLUP)
rr <- mixed.solve(y, X = X, Z = W2)
names(rr)

rr$beta  # fixed effect

rr$Vu  # additive SNP variance 
rr$Ve  # residual variance 

head(rr$u)  # BLUP of SNP marker effects 
tail(rr$u)


# ASREML
load("data/sheep_ped.RData")
?ainverse
sheep_ped_asr <- sheep_ped[c(1,3,2)] %>% 
        as.data.frame() %>% 
        rename(id = ID)
sheep_ainv <- asreml::ainverse(sheep_ped_asr)

sheep_dat

m_basic_1<-asreml(fixed=survival ~ 1 ,
                  random= ~ ~vm(id, sheep_ainv),
                  data=sheep_dat,
                  na.action = na.method(x=c("omit")))
summary(m_basic_1)$varcomp

plot(modelGBLUP)
summary(modelGBLUP)$varcomp
(h2<-nadiv:::pin(modelGBLUP,h2_1~V1/(V1+V2)))

## Obtaining Predictions - BLUP
predGBLUP<-predict(m_basic_1,classify="id", workspace = "20gb")
predGBLUP$predictions$pvals
head(predGBLUP)
View(predGBLUP)

preds<-as.matrix(predGBLUP[,2])
head(preds)
(corr_pearson<-cor(datag$WT_MOT,preds,method='pearson',use="complete.obs"))
plot(datag$WT_MOT,preds)



