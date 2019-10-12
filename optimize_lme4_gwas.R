run_gwas <- function(snp, data) {
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + ", snp, "+ ", paste0("roh_", snp), "+ (1|birth_year) + (1|sheep_year) + (1|id)"))
        mod <- glmer(formula = formula_snp,
                     data = data, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod)
        out
}

snps <- c("oar3_OAR20_137796", "oar3_OAR20_149702", "oar3_OAR20_200365")

#data_mod <- data %>% select(id:sheep_year, oar3_OAR20_137796, roh_oar3_OAR20_137796)
surv <- early_survival_gwas %>% select(id:mum_id, snps) %>% 
        mutate(survival = as.factor(survival)) %>% 
        mutate(id = as.factor(id), twin = as.factor(twin), mum_id = as.factor(mum_id)) %>% 
        drop_na() %>% 
        as.data.frame()


fitGene <- function(snp) {
        f <- reformulate(c(snp, "sex", "twin",  "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
        glmer(f,data=surv, family = "binomial",
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
}

fitGene0 <- function(snp) {
        f <- reformulate(c(snp, "sex", "twin",  "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
        glmer(f,data=surv, family = "binomial",
              nAGQ=0, control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
}

mod1 <- fitGene("oar3_OAR20_137796")
mod2 <- fitGene0("oar3_OAR20_137796")

summary(mod1)
summary(mod2)
glmod0 <- glFormula(survival~ 1 + sex + twin + (1|birth_year) + (1|sheep_year) + (1|id), 
                   data = surv, family = binomial) ## basic structure

refitGene <- function(snp,retmod=TRUE) {
        ## set up new fixed and full formulas
        f0 <- reformulate(c(snp, "sex", "twin"),response="survival")
        f <- reformulate(c(snp, "sex", "twin", "(1|birth_year)", "(1|sheep_year)", "(1|id)"),response="survival")
        ## copy baseline structure and replace relevant pieces
        glmod <- glmod0
        glmod$formula <- f
        glmod$X <- model.matrix(f0,data=surv)
        ## now finish the fit (construct dev fun, optimize,
        ##  optionally return the full model)
        # 1. Create the deviance function for optimizing over theta
        devfun <- do.call(mkGlmerDevfun, glmod)
        # 2. optimize over theta using a rough approximation (i.e. nAGQ = 0)
        opt <- optimizeGlmer(devfun)
        # 3. Update the deviance function for optimizing over theta and beta:
        devfun <- updateGlmerDevfun(devfun, glmod$reTrms)
        # 4. Optimize over theta and beta:
        opt <- optimizeGlmer(devfun, stage=2, optimizer = "nloptwrap", calc.derivs = FALSE)
        
        if (!retmod) opt else {
                mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
        }
        
}

fitAll <- function(...) {
        lapply(grep("^oar",names(surv),value=TRUE),...)
}



out <- benchmark(fitGene(snp = "oar3_OAR20_137796"), replications = 5)
out2 <- benchmark(refitGene(snp = "oar3_OAR20_137796"), replications = 5)


sumfun <- function(x) {
        c(logLik(x),fixef(x),getME(x,"theta"))
}
all.equal(sumfun(fitGene( "oar3_OAR20_137796")),sumfun(refitGene( "oar3_OAR20_137796")))


benchmark(fitAll(fitGene),
          fitAll(refitGene),
          #fitAll(refitGene,retmod=FALSE),
          replications=3,
          columns = c("test", "replications", "elapsed", 
                      "relative"))








