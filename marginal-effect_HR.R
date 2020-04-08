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

froh_all10_cent
invlogit <- function(x) exp(x)/(1+exp(x))

# survival ~ froh_all10_cent * age_cent + froh_all10_cent * lamb + sex + twin

fun <- function(...) {
    ## 'x' is here the regression coefficient: beta.x, not the covariate itself. see
    ## ?inla.posterior.sample.eval
    one <-   invlogit(Intercept + 1 * lamb1 + (-0.5 * froh_all10_cent) + (-2.4 * age_cent) + (-0.5 * 1) * `froh_all10_cent:lamb1` + (-0.5 * -2.4) * `froh_all10_cent:age_cent`)
    two <-   invlogit(Intercept + 1 * lamb1 +    (0 * froh_all10_cent) + (-2.4 * age_cent) +    (0 * 1) * `froh_all10_cent:lamb1` + (0 * -2.4) * `froh_all10_cent:age_cent`)
    three <- invlogit(Intercept + 1 * lamb1 +  (1 * froh_all10_cent) + (-2.4 * age_cent) +  (1 * 1) * `froh_all10_cent:lamb1` + (0.5 * -2.4) * `froh_all10_cent:age_cent`)
    
    four <-   invlogit(Intercept + 0 * lamb1 + (-0.5 * froh_all10_cent) + (-1.4 * age_cent) + (-0.5 * 0) * `froh_all10_cent:lamb1` + (-0.5 * -1.4) * `froh_all10_cent:age_cent`)
    five <-   invlogit(Intercept + 0 * lamb1 +    (0 * froh_all10_cent) + (-1.4 * age_cent) +    (0 * 0) * `froh_all10_cent:lamb1` + (0 * -1.4) * `froh_all10_cent:age_cent`)
    six <-    invlogit(Intercept + 0 * lamb1 +  (1 * froh_all10_cent) + (-1.4 * age_cent) +  (1 * 0) * `froh_all10_cent:lamb1` + (0.5 * -1.4) * `froh_all10_cent:age_cent`)
    
    seven <-   invlogit(Intercept + 0 * lamb1 + (-0.5 * froh_all10_cent) + (2.6 * age_cent) + (-0.5 * 0) * `froh_all10_cent:lamb1` + (-0.5 * 2.6) * `froh_all10_cent:age_cent`)
    eight <-   invlogit(Intercept + 0 * lamb1 +    (0 * froh_all10_cent) + (2.6 * age_cent) +    (0 * 0) * `froh_all10_cent:lamb1` + (0 * 2.6) * `froh_all10_cent:age_cent`)
    nine <-    invlogit(Intercept + 0 * lamb1 +  (1 * froh_all10_cent) + (2.6 * age_cent) +  (1 * 0) * `froh_all10_cent:lamb1` + (0.5 * 2.6) * `froh_all10_cent:age_cent`)
    
    return (list(one, two, three, four, five, six, seven, eight, nine))
}

mod_inla <- readRDS("output/AS_mod_INLA.rds")
xx <- inla.posterior.sample(10, mod_inla)
d = inla.posterior.sample.eval(fun, xx)
d
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



