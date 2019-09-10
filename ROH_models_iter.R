# # run on server
library(MCMCglmm)
library(tidyverse)
library(asreml)
library(tidyr)
library(nadiv)
library(broom.mixed)
library(AnimalINLA)
library(ggregplot)
library(INLA)
library(brinla)
source("theme_clean.R")

# data
load("model_in/fitness_roh_df.RData")
load("model_in/sheep_ped.RData")
IDs_lots_missing <- read_delim("data/ids_more_than_5perc_missing_imputation.txt", delim = " ")

LBS <- fitness_data %>% 
        # some individuals arent imputed well and should be discarded 
        rename(id = ID) %>% 
        filter(!(id %in% IDs_lots_missing)) %>% 
        group_by(id) %>% 
        mutate(counts = n()) %>% 
        summarise(lbs = sum(OffspringBorn, na.rm = TRUE),
                  lifespan = max(Age),
                  birth_year = first(BIRTHYEAR),
                  twin = first(TWIN),
                  sex = first(SEX),
                  mum_id = first(MOTHER),
                  froh_short =  first(FROH_short),
                  froh_medium = first(FROH_medium),
                  froh_long =   first(FROH_long),
                  froh_all =    first(FROH_all),
                  froh_not_roh = first(hom),
                  counts = first(counts)) %>% 
        filter(!is.na(froh_all)) %>% 
        # exclude under 1 year olds
        filter(lifespan > 0) %>% 
        # filter(sex == "M") %>% 
        filter(!(is.na(birth_year) | is.na(mum_id))) %>% 
        mutate(olre = factor(1:nrow(.))) %>% 
        mutate_at(c("id", "birth_year", "mum_id", "sex"), as.factor) %>% 
        mutate(froh_long_std = scale(froh_long),
               froh_medium_std = scale(froh_medium),
               froh_short_std = scale(froh_short),
               froh_not_roh_std = scale(froh_not_roh)) %>% 
        as.data.frame() 

# general plots
ggplot(LBS, aes(lifespan)) + geom_histogram() + facet_wrap(sex ~ ., scales = "free")
ggplot(LBS, aes(lbs)) + geom_histogram() + facet_wrap(sex ~ ., scales = "free")
ggplot(LBS, aes(lifespan, lbs, color = sex)) + geom_jitter(alpha = 0.3) + geom_smooth()

traits <- fitness_data %>% 
        # some individuals arent imputed well and should be discarded 
        filter(!(ID %in% IDs_lots_missing)) %>% 
        dplyr::rename(birth_year = BIRTHYEAR,
                      sheep_year = SheepYear,
                      age = Age,
                      id = ID,
                      twin = TWIN,
                      sex = SEX,
                      mum_id = MOTHER,
                      weight = Weight,
                      hindleg = Hindleg,
                      hornlen = HornLen, 
                      froh_short = FROH_short,
                      froh_medium = FROH_medium,
                      froh_long = FROH_long,
                      froh_all = FROH_all,
                      froh_not_roh = hom) %>% 
        mutate_at(c("id", "birth_year", "mum_id", "sheep_year"), as.factor) %>% 
        dplyr::select(id, sheep_year, age, birth_year, sex, mum_id, twin, weight,
                      hindleg, froh_short, froh_medium, froh_long, froh_all,
                      froh_not_roh, CapMonth) %>% 
        # estimate from yearlings onwards
        filter(age > 0) %>% 
        filter(CapMonth == 8) %>% 
        filter(!is.na(sex) & !is.na(mum_id) & !is.na(birth_year) & !is.na(froh_long),
               !is.na(sheep_year), !is.na(age)) %>% 
        mutate_at(c("id", "birth_year", "mum_id", "sex"), as.factor) %>% 
        mutate(froh_long_std =   scale(froh_long),
               froh_medium_std = scale(froh_medium),
               froh_short_std =  scale(froh_short),
               froh_not_roh_std =scale(froh_not_roh),
               age2 = age^2) 


#~~~~~~~~~~~ Response: 4-month weight, Gaussian. Comparison ASREML vs INLA ~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ROH vs weight
traits %>% 
  dplyr::select(froh_short, froh_medium, froh_long, froh_all, weight, hindleg, sex) %>% 
  tidyr::pivot_longer(cols = starts_with("froh"), names_to = "froh", values_to = "prop") %>% 
  ggplot(aes(prop, weight)) + geom_point() + geom_smooth(method = "lm") + 
  facet_wrap(sex~froh, scales = "free", nrow = 2) +
  theme_clean()

# ROH vs hindleg
traits %>% 
  dplyr::select(froh_short, froh_medium, froh_long, froh_all, weight, hindleg, sex) %>% 
  tidyr::pivot_longer(cols = starts_with("froh"), names_to = "froh", values_to = "prop") %>% 
  ggplot(aes(prop, hindleg)) + geom_point() + geom_smooth(method = "lm") + 
  facet_wrap(sex~froh, scales = "free", nrow = 2) +
  theme_clean()

# predictors: 
# fixed: sex, froh (maybe Litter size, Maternal age quadratic), random: birth_year, mum_id, additive genetic

#~~~~~~~~~~~ Gaussian weight models as animal models ~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Gaussian models for hindleg and weight
# inla
sheep_ped_inla <- sheep_ped %>% 
                        as_tibble() %>% 
                        rename(id = ID,
                               mother = MOTHER,
                               father = FATHER) %>% 
                        mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "F", "888")) %>% 
                        mutate_at(c("id", "mother", "father"), function(x) str_replace(x, "M", "999")) %>% 
                        #filter(!is.na(id)) %>% 
                        mutate(father = ifelse(is.na(father), 0, father)) %>% 
                        mutate(mother = ifelse(is.na(mother), 0, mother)) %>% 
                        mutate_if(is.character, list(as.numeric)) %>% 
                        as.data.frame() 

comp_inv <- AnimalINLA::compute.Ainverse(sheep_ped_inla)
ainv <- comp_inv$Ainverse
ainv_map <- comp_inv$map
Cmatrix <- sparseMatrix(i=ainv[,1],j=ainv[,2],x=ainv[,3])



## Mapping the same index number for "Individual" as in Ainv
## The IndexA column is the index in the A inverse matrix
traits_hindleg <- traits %>% 
  dplyr::filter(!is.na(hindleg))

traits_weight <- traits %>% 
  dplyr::filter(!is.na(weight))

add_index_inla <- function(traits) {
  Ndata <- dim(traits)[1]
  traits$IndexA <- rep(0, times = Ndata)
  for(i in 1:Ndata) {
    traits$IndexA[i] <- which(ainv_map[,1]==traits$id[i])
  }
  traits
}

traits_hindleg <- add_index_inla(traits_hindleg)
traits_weight<- add_index_inla(traits_weight)

# which prior? Same as for additive genetic variance seems to work well
# for all other random effects too
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))

# fit a few models 
froh_combinations <- c("froh_long + froh_medium + froh_short",
                       "froh_long_std + froh_medium_std + froh_short_std",
                       "froh_long",
                       "froh_medium",
                       "froh_short",
                       "froh_all")

# this part of the formula stays the same
formula_base <- paste('sex + age + age2 + twin + 1', 
                'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                'f(sheep_year, model = "iid", hyper = list(prec = prec_prior))',
                'f(id, model = "iid", hyper = list(prec = prec_prior))',
                'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + ")

# make all formulas for hindleg models
all_formulas_hindleg <- map(paste('hindleg ~', froh_combinations, ' + ', formula_base), as.formula) # env=globalenv()

# model inla for gaussian traits
inla_gaussian_traits <- function(formula, data) {
    mod_inla <- inla(formula=formula, family="gaussian",
                     data=data, 
                     control.compute = list(dic = TRUE))
}

# run models 
all_inla_mods_hindleg <- map(all_formulas_hindleg, inla_gaussian_traits, traits_hindleg)

# make all formulas for weight models
all_formulas_weight  <- map(paste('weight ~', froh_combinations, ' + ', formula_base), as.formula)

# run models
all_inla_mods_weight <- map(all_formulas_weight, inla_gaussian_traits, traits_weight)

# combine into df
gaussian_mods <- tibble(fitted = c(all_inla_mods_hindleg, all_inla_mods_weight), 
                        formulas = c(all_formulas_hindleg, all_formulas_weight),
                        response = c(rep("hindleg", 6), rep("weight", 6)))
saveRDS(gaussian_mods, file = "output/models/gaussian_mods.rds")


map(all_inla_mods_hindleg, function(x) x$summary.fixed)
map(all_inla_mods_weight, function(x) x$summary.fixed)

summary(mod_inla0)
bri.hyperpar.summary(mod_inla2)^2
hrtblt <- (bri.hyperpar.summary(mod_inla0)[6, ] / sum(bri.hyperpar.summary(mod_inla0)[, 1])) 




#~~~~~~~~~~~ LBS models ~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ROH vs LBS
LBS %>% 
  dplyr::select(froh_short, froh_medium, froh_long, froh_all, lbs, lifespan, sex) %>% 
  tidyr::pivot_longer(cols = starts_with("froh"), names_to = "froh", values_to = "prop") %>% 
  ggplot(aes(lbs, prop)) + geom_point() + geom_smooth(method = "lm") + 
  facet_wrap(sex~froh, scales = "free", nrow= 2) 

LBS_m <- LBS %>% 
  dplyr::filter(!is.na(lbs)) %>% 
  dplyr::filter(sex == "M") %>% 
  mutate(olre = 1:nrow(.))

LBS_f <- LBS %>% 
  dplyr::filter(!is.na(lbs)) %>% 
  dplyr::filter(sex == "F") %>% 
  mutate(olre = 1:nrow(.))

add_index_inla <- function(traits) {
  Ndata <- dim(traits)[1]
  traits$IndexA <- rep(0, times = Ndata)
  for(i in 1:Ndata) {
    traits$IndexA[i] <- which(ainv_map[,1]==traits$id[i])
  }
  traits
}

LBS_m <- add_index_inla(LBS_m)
LBS_f <- add_index_inla(LBS_f)

# prior same for all random effects
prec_prior <- list(prior = "loggamma", param = c(0.5, 0.5))

# fit a few models 
froh_combinations <- c("froh_long + froh_medium + froh_short",
                       "froh_long_std + froh_medium_std + froh_short_std",
                       "froh_long",
                       "froh_medium",
                       "froh_short",
                       "froh_all")

# this part of the formula stays the same
formula_base <- paste('lifespan + twin + 1', 
                      'f(birth_year, model = "iid", hyper = list(prec = prec_prior))',
                      'f(mum_id, model="iid",  hyper = list(prec = prec_prior))', 
                      'f(olre, model="iid",  hyper = list(prec = prec_prior))', 
                      'f(IndexA, model="generic0", hyper = list(theta = list(param = c(0.5, 0.5))),Cmatrix=Cmatrix)', sep = " + ")

# make all formulas for hindleg models
all_formulas_lbs <- map(paste('lbs ~', froh_combinations, ' + ', formula_base), as.formula) # env=globalenv()




formula <- lbs ~ 1  + froh_long + froh_medium + froh_short + lifespan + #  + froh_long_std + froh_medium_std + froh_short_std + froh_not_roh_std 
  f(birth_year, model = "iid", hyper = list(prec = prec_prior)) +
  f(mum_id, model="iid", hyper = list(prec = prec_prior)) + 
  f(olre, model="iid", hyper = list(prec = prec_prior)) +
  f(IndexA ,model="generic0", Cmatrix=Cmatrix,
    constr=TRUE, param = c(0.5, 0.5)) #+

mod_inla1 <- inla(formula=formula, family="poisson",
                 data=LBS_m, 
                 control.compute = list(dic = TRUE))


summary(mod_inla1)

exp(0.21)
EY <- mean(LBS_female$LBS)
estdv_link = log(1/EY+1)
library(QGglmm)
QGparams(mu=-1.5232, var.a=0.274242132^2, var.p = sum(bri.hyperpar.summary(mod_inla1)[, 1])^2,
                              model = "Poisson.log")

(bri.hyperpar.summary(mod_inla1)[4, 1] / (sum(bri.hyperpar.summary(mod_inla1)[, 1][-1]) + 1/ +  estdv_link))

mreff <- mod_inla1$summary.random$IndexA$mean
qqnorm(mreff)
qqline(mreff)

mod_inla$marginals.fixed[[3]] %>% 
  as.data.frame() %>% 
  ggplot(aes(x,y)) + geom_line() 
bri.hyperpar.plot(mod_inla)

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

formula <- LBS ~ 1 + froh_all + 
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
                  control.family=list(hyper = list(theta = list(param = c(0.5, 0.5), fixed = FALSE))),
                  only.hyperparam =FALSE,control.compute=list(dic=T))

bri.fixed.plot(mod_inla_lbs3)
brinla::bri.hyperpar.plot(mod_inla_lbs3)

## rstanarm
library(rstanarm)

traits_sub <- traits %>% sample_frac( size = 0.05)
mod_rstan <- stan_lmer( weight ~ 1 + froh_long + (1|MumID) + (1|birth_year),
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
                   fixed = c("froh_long"),
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
mod_inla_org <- inla(formula = Weight ~ 1 + froh_long + 
                             f(MOTHER, model = "iid") + 
                             f(BIRTHYEAR, model = "iid") + 
                             f(ID, model = "generic2", 
                               constr = TRUE, Cmatrix = Cmatrix),
                             family = "gaussian",
                             data = hindleg)

mod2 <- inla(formula = Weight ~ 1 + froh_long +
                     f(BIRTHYEAR, model = "iid", constr = TRUE, hyper = list(theta = list(param = c(1, 0.001), fixed = FALSE))) + 
                     f(MOTHER, model = "iid", constr = TRUE, hyper = list(theta = list(param = c(1, 0.001), fixed = FALSE))) +
                     f(ID, model = "generic0", hyper = list(theta = list(param = c(0.5, 0.5), fixed = FALSE)), constr = TRUE, 
                       Cmatrix = Cmatrix), 
             family = "gaussian",
             data = hindleg)



out <- animal.inla(response = "LBS",
                   fixed = "froh_long",
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
inla_form <- LBS ~ f(froh_long, model = "iid", constr = TRUE, hyper = list(theta = prec.Fixed)) + 
        f(birth_year, model = "iid", constr = TRUE, hyper = list(theta = prec.Random)) + 
        f(ID, model = "generic0", hyper = list(theta = prec.A), constr = TRUE, 
          Cmatrix = Cmatrix) + 1

prec.Fixed <- list(initial = -10, fixed = TRUE)
out$Cmatrix
LBS_male$ID

dat <- out$dataset
out2 <- inla(formula = LBS ~ 1 + froh_long + froh_medium + froh_short +
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
# 1 + froh, random = ~vm(sheep_ped) + idv(MumID) + idv(BirthYear) + idv(DeathYear)
asreml.options(maxit = 30)
mod_all <- asreml(fixed = LBS ~ 1 + froh_long + froh_medium , random = ~vm(ID, sheep_ainv) +  idv(MumID) + idv(birth_year) + idv(death_year) , # + idv(olre)
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

mod1 <- MCMCglmm(LBS ~ 1 + froh_long + froh_medium + froh_short, random = ~ ID + IDD, #+ MumID + birth_year + death_year
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
                mutate_at(c("froh_short", "froh_long", "froh_medium"), scale)

mod2 <- MCMCglmm(LBS ~ 1 + froh_long + froh_medium + froh_short, random = ~animal + MumID + birth_year + death_year, #+ MumID + BirthYear + DeathYear, 
                 family = "poisson",
                 prior = prior, pedigree = sheep_ped_MCMCglmm, 
                 data = LBS_male_scaled, 
                 nitt = nitt, burnin = burnin, thin = thin)

tidy(mod2, conf.int = TRUE)

save(mod2, file = "model_out/all_roh_poisson2.RData")


mod3 <- MCMCglmm(LBS ~ 1 + froh_long + froh_medium + froh_short, random = ~animal + MumID + birth_year, #+ MumID + BirthYear + DeathYear, 
                 family = "poisson",
                 prior = prior, pedigree = sheep_ped_MCMCglmm, 
                 data = LBS_male, 
                 nitt = nitt, burnin = burnin, thin = thin)

save(mod3, file = "model_out/all_roh_poisson3.RData")

mod4 <- MCMCglmm(LBS ~ 1 + froh_long + froh_medium + froh_short, random = ~animal, #+ MumID + BirthYear + DeathYear, 
                 family = "poisson",
                 prior = prior, pedigree = sheep_ped_MCMCglmm, 
                 data = LBS_male, 
                 nitt = nitt, burnin = burnin, thin = thin)

save(mod4, file = "model_out/all_roh_poisson4.RData")

x



