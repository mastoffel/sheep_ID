# # run on server
library(MCMCglmm)
library(tidyverse)
library(asreml)
library(nadiv)
library(broom.mixed)
library(AnimalINLA)
library(ggregplot)
library(INLAutils)
library(INLAOutputs)
load("model_in/fitness_roh_df.RData")
load("model_in/sheep_ped.RData")

# ped
# inv.phylo <- MCMCglmm::inverseA(sheep_ped)
# sheep_ped_inv <- solve(inv.phylo$Ainv)
# rownames(sheep_ped_inv) <- rownames(inv.phylo$Ainv)

LBS <- fitness_data %>% 
        group_by(ID) %>% 
        summarise(LBS = sum(OffspringBorn, na.rm = TRUE),
                  birth_year = first(BIRTHYEAR),
                  death_year = first(DeathYear),
                  sex = first(SEX),
                  MumID = first(MOTHER),
                  FROH_short = first(FROH_short),
                  FROH_medium = first(FROH_medium),
                  FROH_long = first(FROH_long),
                  FROH_all = first(FROH_all)) %>% 
        filter(!is.na(FROH_all)) %>% 
        # filter(sex == "M") %>% 
        filter(!(is.na(birth_year) | is.na(death_year) | is.na(MumID))) %>% 
        mutate(olre = factor(1:nrow(.))) %>% 
        mutate(ID = as.factor(ID)) %>% 
        mutate_at(c("birth_year", "death_year", "MumID"), as.factor) %>% 
        as.data.frame() 

traits <- fitness_data %>% 
        filter(CapMonth == 8) %>% 
        filter(SheepYear == BIRTHYEAR) %>% 
        mutate(ID = as.factor(ID), SEX = as.factor(SEX)) %>% 
        mutate(MOTHER = as.factor(MOTHER),
        BIRTHYEAR = as.factor(BIRTHYEAR)) %>% 
        dplyr::rename(birth_year = BIRTHYEAR,
               sex = SEX,
               MumID = MOTHER,
               weight = Weight,
               hindleg = Hindleg)


####################### ASREML ######################
# pedigree format: Needs to be sorted, and individual, father, mother
sheep_ped_asr <- read_delim("../sheep/data/SNP_chip/20190208_Full_Pedigree.txt", delim = "\t")[c(1,3,2)] %>% 
        as.data.frame() 
    #    filter(ID %in% hindleg $ID) %>% 
#        orderPed() 
# inverse relationship mat
sheep_ainv <- asreml::ainverse(sheep_ped_asr)

mod_asr <- asreml(fixed = weight ~ 1 + FROH_long, random = ~vm(ID, sheep_ainv) + idv(MumID) + idv(birth_year), 
              data = traits, na.action = na.method(x=c("omit"))) # "omit" "include"

mod_sum <- summary(mod_asr, coef = TRUE)
mod_sum$varcomp[3,1] / sum(mod_sum$varcomp[,1]) # var(fitted(mod_asr), na.rm = TRUE))



######################## INLA ######################

sheep_ped_inla <- sheep_ped %>% 
                        mutate(Father = ifelse(is.na(Father), 0, Father)) %>% 
                        mutate(Mother = ifelse(is.na(Mother), 0, Mother)) %>% 
                        as.data.frame()

xx = compute.Ainverse(sheep_ped_inla)
Ainv = xx$Ainverse
map = xx$map
Cmatrix = sparseMatrix(i=Ainv[,1],j=Ainv[,2],x=Ainv[,3])
Ndata = dim(traits)[1]

## Mapping the same index number for "Individual" as in Ainv
## The IndexA column is the index in the A inverse matrix
traits$IndexA <- rep(0,Ndata)
for(i in 1:Ndata) traits$IndexA[i] = which(map[,1]==traits$ID[i])

#Including an extra column for individual effect
#sparrowGaussian$IndexA.2=sparrowGaussian$IndexA

traits <- traits %>% mutate(
  FROH_short_std = scale(FROH_short),
  FROH_medium_std = scale(FROH_medium),
  FROH_long_std = scale(FROH_long)
)

formula <- hindleg ~ 1 + FROH_long_std + sex + 
                f(birth_year, model = "iid") +
                f(MumID, model="iid") + 
                f(IndexA ,model="generic0", Cmatrix=Cmatrix,
                constr=TRUE,param = c(0.5, 0.5)) #+
      #  f(IndexA.2,model="iid",param = c(1,0.001), constr=TRUE)

# y, IndexA and IndexA.2 is the individuals in the
# data (these have to be given different names) where IndexA is the additive genetic effect and IndexA.2 is
# the individual random effect.
mod_inla3 <- inla(formula=formula, family="gaussian",
             data=traits,
             control.predictor=list(compute=TRUE), 
             control.family=list(hyper = list(theta = list(param = c(0.5, 0.5), fixed = FALSE))),
             only.hyperparam =FALSE,control.compute=list(dic=T))

summary(mod_inla1)
bri.fixed.plot(mod_inla1)
brinla::bri.hyperpar.plot(mod_inla)


bri.fixed.plot(mod_inla1)
bri.fixed.plot(mod_inla3)
bri.fixed.plot(mod_inla2)

invsqrt <- function(x) 1/sqrt(x)
sdt <- invsqrt(mod_inla$summary.hyperpar[,-2]) # no sd
row.names(sdt) <- c("SD of epsilson","SD of birth year", "SD of MumID", "SD of additive genetic")
prec.birth_year <- mod_inla$marginals.hyperpar$"Precision for birth_year"
prec.epsilon <- imod$marginals.hyperpar$"Precision for the Gaussian observations"
c(epsilon=inla.emarginal(invsqrt,prec.epsilon),
  site=inla.emarginal(invsqrt,prec.birth_year))
sigma.birth_year <- inla.tmarginal(invsqrt, prec.birth_year)
sigma.epsilon <- inla.tmarginal(invsqrt, prec.epsilon)
sampvars <- 1/inla.hyperpar.sample(1000,mod_inla)
sampicc <- sampvars[,2]/(rowSums(sampvars))
quantile(sampicc, c(0.025,0.5,0.975))
bri.hyperpar.summary(mod_inla)
alpha <- data.frame(mod_inla$marginals.fixed[[1]])

mod_inla$internal.marginals.hyperpar
mod_sum <-INLARep(mod_inla)
mod_sum[3,1] / sum(mod_sum$Mean) 
mod_inla$summary.random
mod_inla$summary.fixed
mod_inla$marginals.hyperpar

plot(mod_inla, plot.fixed.effects = TRUE, plot.lincomb = FALSE, plot.random.effects = FALSE, plot.hyperparameters = TRUE,
     plot.predictor = FALSE, plot.q = FALSE, plot.cpo = FALSE, single = FALSE)




# inla LBS
LBS_male <- LBS %>% filter(sex == 'M')
Ndata = dim(LBS_male)[1]
LBS_male$IndexA <- rep(0,Ndata)
for(i in 1:Ndata) LBS_male$IndexA[i] = which(map[,1]==LBS_male$ID[i])

formula <- LBS ~ 1 + FROH_short +
  f(olre, model = "iid") +
  f(birth_year, model = "iid") +
  f(MumID, model="iid") + 
  f(IndexA, model = "generic0", constr = TRUE, 
    Cmatrix = Cmatrix)
  # f(IndexA ,model="generic0", Cmatrix=Cmatrix,
  #   constr=TRUE,param = c(0.5, 0.5)) #+

# the individual random effect.
mod_inla_lbs3 <- inla(formula=formula, family="poisson",
                  data=LBS_male,
                  control.predictor=list(compute=TRUE), 
                 # control.family=list(hyper = list(theta = list(param = c(0.5, 0.5), fixed = FALSE))),
                  only.hyperparam =FALSE,control.compute=list(dic=T))

bri.fixed.plot(mod_inla_lbs3)
brinla::bri.hyperpar.plot(mod_inla_lbs3)

## rstanarm
library(rstanarm)

traits_sub <- traits %>% sample_frac( size = 0.05)
mod_rstan <- stan_lmer( weight ~ 1 + FROH_long + (1|MumID) + (1|birth_year),
                        data = traits_sub,
                        seed = 355)
summary(mod_rstan)
tidy(mod_rstan)



#Example finding the posterior marginal distribution and mean (95% CI)
#for additive genetic variance and individual random variance
# sigma.IndexA = inla.marginal.transform(function(x) 1/x,
#                                        model$marginals.hyperpar$`Precision for IndexA`)
# e.IndexA=inla.expectation(function(x) x, sigma.IndexA)
# ci.IndexA=inla.qmarginal(c(0.025, 0.975), sigma.IndexA)
# model$marginals.hyperpar$`Precision for IndexA`

mod_inla <- animal.inla(response = Weight,
                   fixed = c("FROH_long"),
                   #random = c("BIRTHYEAR"),
                   genetic = c("ID"),
                   data = hindleg,
                   Ainverse =sheep_ainv_inla,
                   type.data = "gaussian")

mod_inla$formula
mod_inla$inla.call
summary(mod_inla)

Ainv <- sheep_ainv_inla$Ainverse
Cmatrix = sparseMatrix(i = Ainv[, 1], j = Ainv[, 2],
                      x = Ainv[, 3])

mod_inla$prec.A
mod_inla_org <- inla(formula = Weight ~ 1 + FROH_long + 
                             f(MOTHER, model = "iid") + 
                             f(BIRTHYEAR, model = "iid") + 
                             f(ID, model = "generic2", 
                               constr = TRUE, Cmatrix = Cmatrix),
                             family = "gaussian",
                             data = hindleg)

mod2 <- inla(formula = Weight ~ 1 + FROH_long +
                     f(BIRTHYEAR, model = "iid", constr = TRUE, hyper = list(theta = list(param = c(1, 0.001), fixed = FALSE))) + 
                     f(MOTHER, model = "iid", constr = TRUE, hyper = list(theta = list(param = c(1, 0.001), fixed = FALSE))) +
                     f(ID, model = "generic0", hyper = list(theta = list(param = c(0.5, 0.5), fixed = FALSE)), constr = TRUE, 
                       Cmatrix = Cmatrix), 
             family = "gaussian",
             data = hindleg)



out <- animal.inla(response = "LBS",
                   fixed = "FROH_long",
                   random = "birth_year",
                   genetic = c("ID"),
                   data = LBS_male,
                   Ainverse = sheep_ainv_inla,
                   type.data = "poisson")

out$formula
out$inla.call
out$model.covariates
out$summary.covariates
out$name.covariates
out$prec.fixed.generic2
summary(out)
summary(out)$Model.summary
plot.Animalinla(out)

library(INLA)
inla_form <- LBS ~ f(FROH_long, model = "iid", constr = TRUE, hyper = list(theta = prec.Fixed)) + 
        f(birth_year, model = "iid", constr = TRUE, hyper = list(theta = prec.Random)) + 
        f(ID, model = "generic0", hyper = list(theta = prec.A), constr = TRUE, 
          Cmatrix = Cmatrix) + 1

prec.Fixed <- list(initial = -10, fixed = TRUE)
out$Cmatrix
LBS_male$ID

dat <- out$dataset
out2 <- inla(formula = LBS ~ 1 + FROH_long + FROH_medium + FROH_short +
                      f(olre, model = "iid", constr = TRUE, hyper = list(theta = list(param = c(1, 0.001), fixed = FALSE))) + 
                     # f(death_year = "iid", constr = TRUE, hyper = list(theta = list(param = c(1, 0.001), fixed = FALSE))) + 
                      f(birth_year, model = "iid", constr = TRUE, hyper = list(theta = list(param = c(1, 0.001), fixed = FALSE))) + 
                      f(MumID, model = "iid", constr = TRUE, hyper = list(theta = list(param = c(1, 0.001), fixed = FALSE))) +
                      f(ID, model = "generic0", hyper = list(theta = list(param = c(0.5, 0.5), fixed = FALSE)), constr = TRUE, 
                        Cmatrix = Cmatrix), 
              family = "poisson",
              data = dat)

summary(out2)
################################ asreml ########################################

# inverse for asreml
sheep_ainv <- asreml::ainverse(sheep_ped)
sheep_dinv <- nadiv::makeD(sheep_ped)
sheep_dinv_asreml <- sheep_dinv$Dinv
table(sheep_dinv$listDinv$Dinverse)
LBS_male$IDD <- LBS_male$ID

#attr(sheep_ainv, "rowNames")
# 1 + FROH, random = ~vm(sheep_ped) + idv(MumID) + idv(BirthYear) + idv(DeathYear)
asreml.options(maxit = 30)
mod_all <- asreml(fixed = LBS ~ 1 + FROH_long + FROH_medium , random = ~vm(ID, sheep_ainv) +  idv(MumID) + idv(birth_year) + idv(death_year) , # + idv(olre)
                  data = LBS_male,
                  family = asr_gaussian(),
                  na.action = na.method(x=c("omit")))

summary(mod_all, coef = TRUE)$coef.fixed

mod_sums <- function(mod) {
        mod_sum <- summary(mod, coef = TRUE)
        out <- list(random = mod_sum$varcomp, fixed = mod_sum$coef.fixed)
}

all_mods <- map(list(mod_all, mod_short, mod_medium, mod_long), mod_sums)
map(all_mods, function(x) x$fixed)





# MCMCglmm
LBS_male <- LBS_male %>% 
        rename(animal = ID)




# model0 MCMCglmm
sheep_ped_MCMCglmm <- sheep_ped %>% 
        rename(animal = ID)

Ainv <- inverseA(sheep_ped)$Ainv 
Dinv <- makeD(sheep_ped)$Dinv 
LBS_male$IDD <- LBS_male$ID

#
LBS_dat <- LBS_male %>% 
        sample_frac(0.1)

#phenvar <- var(lrt_roh_df_mod1$n_offspring)
prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002),
                                                G2 = list(V=1, nu=0.002),
                                                G3 = list(V=1, nu=0.002),
                                                G4 = list(V=1, nu=0.002)))

prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002),
                                                G2 = list(V=1, nu=0.002)))

nitt <- 11000
burnin <- 1000
thin <- 10

mod1 <- MCMCglmm(LBS ~ 1 + FROH_long + FROH_medium + FROH_short, random = ~ ID + IDD, #+ MumID + birth_year + death_year
                 family = "poisson",
                 ginverse =  list(ID = Ainv, IDD = Dinv),
                 prior = prior, 
                 data = LBS_dat, 
                 nitt = nitt, burnin = burnin, thin = thin)

save(mod1, file = "model_out/all_roh_poisson.RData")

plot(mod1$Sol)
plot(mod1$VCV)
summary(mod1)
tidy(mod1, conf.int = TRUE)



LBS_male_scaled <- LBS_male %>% 
                mutate_at(c("FROH_short", "FROH_long", "FROH_medium"), scale)

mod2 <- MCMCglmm(LBS ~ 1 + FROH_long + FROH_medium + FROH_short, random = ~animal + MumID + birth_year + death_year, #+ MumID + BirthYear + DeathYear, 
                 family = "poisson",
                 prior = prior, pedigree = sheep_ped_MCMCglmm, 
                 data = LBS_male_scaled, 
                 nitt = nitt, burnin = burnin, thin = thin)

tidy(mod2, conf.int = TRUE)

save(mod2, file = "model_out/all_roh_poisson2.RData")


mod3 <- MCMCglmm(LBS ~ 1 + FROH_long + FROH_medium + FROH_short, random = ~animal + MumID + birth_year, #+ MumID + BirthYear + DeathYear, 
                 family = "poisson",
                 prior = prior, pedigree = sheep_ped_MCMCglmm, 
                 data = LBS_male, 
                 nitt = nitt, burnin = burnin, thin = thin)

save(mod3, file = "model_out/all_roh_poisson3.RData")

mod4 <- MCMCglmm(LBS ~ 1 + FROH_long + FROH_medium + FROH_short, random = ~animal, #+ MumID + BirthYear + DeathYear, 
                 family = "poisson",
                 prior = prior, pedigree = sheep_ped_MCMCglmm, 
                 data = LBS_male, 
                 nitt = nitt, burnin = burnin, thin = thin)

save(mod4, file = "model_out/all_roh_poisson4.RData")

x



