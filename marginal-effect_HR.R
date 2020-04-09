library(INLA)
library(tidyverse)
# data
load("data/survival_mods_data.RData") 
load("data/sheep_ped.RData")
ped <- sheep_ped
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Annual survival~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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


froh_all10_cent
invlogit <- function(x) exp(x)/(1+exp(x))

# survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin
table(annual_survival$age_cent)
mean(annual_survival$froh_all10)
hist(annual_survival$froh_all10_cent)

annual_survival %>% 
    summarise(mean(froh_all), min(froh_all), max(froh_all))
annual_survival %>% 
    summarise(mean(froh_all10_cent), min(froh_all10_cent), max(froh_all10_cent))


mod_inla <- readRDS("output/AS_mod_INLA.rds")
fun <- function(...) {
    ## 'x' is here the regression coefficient: beta.x, not the covariate itself. see
    ## ?inla.posterior.sample.eval
    test <- -0.5
    one <-     invlogit(Intercept + 1 * lamb1 + (test * froh_all10_cent) + (-2.4 * age_cent) + 0.15 * twin1 + 0.4 * sexM + (-0.5 * 1) * `froh_all10_cent:lamb1` + (-0.5 * -2.4) * `froh_all10_cent:age_cent`)
    two <-     invlogit(Intercept + 1 * lamb1 +    (0 * froh_all10_cent) + (-2.4 * age_cent) + 0.15 * twin1 + 0.4 * sexM +    (0 * 1) * `froh_all10_cent:lamb1` + (0 * -2.4) * `froh_all10_cent:age_cent`)
    three <-   invlogit(Intercept + 1 * lamb1 +  (1.5 * froh_all10_cent) + (-2.4 * age_cent) + 0.15 * twin1 + 0.4 * sexM +  (1.5 * 1) * `froh_all10_cent:lamb1` + (1.5 * -2.4) * `froh_all10_cent:age_cent`)
  
    four <-    invlogit(Intercept + 0 * lamb1 + (-0.5 * froh_all10_cent) + (-1.4 * age_cent) + 0.15 * twin1 + 0.4 * sexM + (-0.5 * 0) * `froh_all10_cent:lamb1` + (-0.5 * -1.4) * `froh_all10_cent:age_cent`)
    five <-    invlogit(Intercept + 0 * lamb1 +    (0 * froh_all10_cent) + (-1.4 * age_cent) + 0.15 * twin1 + 0.4 * sexM +    (0 * 0) * `froh_all10_cent:lamb1` + (0 * -1.4) * `froh_all10_cent:age_cent`)
    six <-     invlogit(Intercept + 0 * lamb1 +  (1.5 * froh_all10_cent) + (-1.4 * age_cent) + 0.15 * twin1 + 0.4 * sexM +  (1.5 * 0) * `froh_all10_cent:lamb1` + (1.5 * -1.4) * `froh_all10_cent:age_cent`)
  
    seven <-   invlogit(Intercept + 0 * lamb1 + (-0.5 * froh_all10_cent) + (2.6 * age_cent)  + 0.15 * twin1 + 0.4 * sexM + (-0.5 * 0) * `froh_all10_cent:lamb1` + (-0.5 * 2.6) * `froh_all10_cent:age_cent`)
    eight <-   invlogit(Intercept + 0 * lamb1 +    (0 * froh_all10_cent) + (2.6 * age_cent)  + 0.15 * twin1 + 0.4 * sexM +    (0 * 0) * `froh_all10_cent:lamb1` + (0 * 2.6) * `froh_all10_cent:age_cent`)
    nine <-    invlogit(Intercept + 0 * lamb1 +  (1.5 * froh_all10_cent) + (2.6 * age_cent)  + 0.15 * twin1 + 0.4 * sexM +  (1.5 * 0) * `froh_all10_cent:lamb1` + (1.5 * 2.6) * `froh_all10_cent:age_cent`)
    
    return (list(one, two, three, four, five, six, seven, eight, nine))
}

fun <- function(...) {
    one <-  invlogit(Intercept + 
                            df1$x1 * froh_all10_cent + 
                            df1$x2 * age_cent + 
                            df1$x3 * lamb1 + 
                            df1$x4 * twin1 + 
                            df1$x5 * sexM + 
                            df1$x6 * `froh_all10_cent:lamb1` + 
                            df1$x7 * `froh_all10_cent:age_cent`)

    return (list(one))
}

froh <- seq(from = min(annual_survival$froh_all10_cent), to = (max(annual_survival$froh_all10_cent)), by = 0.1)
age <- c(-2.4, -1.4, 2.6, 5.6)
combined_df <- expand_grid(froh, age) %>% 
                mutate(lamb = ifelse(age == -2.4, 1, 0),
                       twin = 0.15,
                       sex = 0.4,
                       frohxlamb = froh*lamb,
                       frohxage = froh*age) 
names(combined_df) <- paste0("x", 1:7)
                       
xx <- inla.posterior.sample(10, mod_inla)
marg_means <- purrr::map(1:nrow(combined_df), function(x) {
    df1 <<- combined_df[x, ]
    out <- inla.posterior.sample.eval(fun, xx)
}) 

d <- marg_means %>% 
        map(as_tibble) %>% 
        bind_rows() %>% 
        pmap_df(function(...) {
          samp <- as.numeric(unlist(list(...)))
          c(mean = mean(samp), quantile(samp, probs = c(0.025, 0.975)))
        }) %>% 
        bind_cols(combined_df)

ggplot(d, aes(x1, mean)) +
  geom_line(aes(color = factor(x2))) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = factor(x2)))

d %>% 
    t() %>% 
    as_tibble() %>% 
    map_df(function(x) c(mean = mean(x), quantile(x, probs = c(0.025, 0.975)))) %>% 
    mutate(froh = rep(c(-0.5, 0, 1.5), 3), age = rep(c(-2.5, -1.4, 2.6), each = 3)) %>% 
    mutate(age = as.factor(age)) %>% 
    ggplot(aes(froh, mean)) +
    geom_line(aes(color = age)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = age)) 

geom_line(data = df_full, aes(x = x, predicted, color = group), size = 1) +
    geom_ribbon(data= df_full, aes(x=x, ymin=conf.low, ymax=conf.high, fill = group, color = group), alpha= 0.2,
                linetype = 2, size = 0.1) +
    scale_color_viridis_d("Age", labels = c(0, 1, 4, 7)) +
    scale_fill_viridis_d("Age", labels = c(0, 1, 4, 7)) +
    theme_simple(grid_lines = FALSE, axis_lines = TRUE) 


n = 100
ng = 10
g = rep(1:ng, each = n %/% ng)
x = rnorm(n)
eta = 1 + x + rnorm(g)
y = eta + rnorm(n, sd=0.2)

r = inla(y ~ 1 + x + f(g, model="iid"),
         family = "gaussian",
         data = data.frame(y, x, g), 
         control.compute = list(config=TRUE))


fun = function(...) {
    ## 'x' is here the regression coefficient: beta.x, not the covariate itself. see
    ## ?inla.posterior.sample.eval
    eta.1 = Intercept + x * 1
    eta.0 = Intercept + x * 0
    diff = eta.1 - eta.0
    return (diff)
}

xx = inla.posterior.sample(1000,  r)
d = inla.posterior.sample.eval(fun, xx)
hist(d, n = 100, prob=TRUE)



## x1 and x2 are continuous covariates, x3 is a factor with 3

n  = 100
ns = 1E5
x1 = rnorm(n)
x2 = rnorm(n)
x3 = sample(c("A", "B", "C"), n,  replace=TRUE)

eta = 1 + x1 + x2 + 1*(x3 == "A") + 2*(x3 == "B") + 3*(x3 == "C")
y = eta + rnorm(n,  sd = 0.001)

r = inla(y ~ 1 + x1 + x2 + x3,
         data = data.frame(x1, x2, x3),
         control.compute = list(config=TRUE))
print(rownames(r$summary.fixed))

xx = inla.posterior.sample(ns, r)
fun = function(cov.x1, cov.x2, cov.x3B, cov.x3C)
{
    ## see the documentation...
    return (Intercept + x1 * cov.x1 + x2 * cov.x2 +
                x3B * cov.x3B + x3C * cov.x3C)
}
zz = inla.posterior.sample.eval(fun,
                                xx,
                                cov.x1 = x1,
                                cov.x2 = x2,
                                cov.x3B = as.numeric(x3 == "B"),
                                cov.x3C = as.numeric(x3 == "C"))

plot(rowMeans(zz),  eta)
abline(a=0, b=1)



