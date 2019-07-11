# ############################################################################
# Variance decomposition workshop
# Bielefeld, 18-22 March 2019
# Holger Schielzeth & Tim Schmoll
#
# Day 5:     Continuous random effects and animal models
# Script 5A: Animal models
# ############################################################################

#_____________________________________________________________________________
# [1] Animal model for modelling neophilia in zebra finches: Preparation
#_____________________________________________________________________________

# Routine clean-up upon start
rm(list=ls())

# Loading of required packages
require(pedantics)
require(pedigreemm)
require(MCMCglmm)

# [1a] The pedigree
#___________________

pedfile = file.choose()
ped = read.table(pedfile, header=TRUE)
names(ped) = c("id", "dam", "sire")

drawPedigree(ped)
round(pedStatSummary(pedigreeStats(ped, graphicalReport="n")),3)

# pedigree format required for pedigreemm
pedS4 = pedigree(sire=as.character(ped$sire), dam=as.character(ped$dam), label=as.character(ped$id))  

# [1b] The data
#___________________

md = read.table(file.choose(), header=TRUE)
str(md)
# Some reformatting required to match individual names with the pedigree file
md$Animal = paste("AA", md$IndID, sep="")
str(md)
md$ObjectType = factor(md$ObjectType)
levels(md$ObjectType) = c("Flower", "Apple")
str(md)

# [1c] The repeatability model
#______________________________

mod = lmer(NeophiliaScore ~ 1 + (1|IndID), data=md)
summary(mod)
# Fitting the covariate
mod = lmer(NeophiliaScore ~ ObjectType + (1|IndID), data=md)
summary(mod)
# Removal of intercept
mod = lmer(NeophiliaScore ~ ObjectType -1 + (1|IndID), data=md)
summary(mod)

data.frame(summary(mod)$varcor)
varcomp = data.frame(summary(mod)$varcor)[,4]
varcomp
names(varcomp) = data.frame(summary(mod)$varcor)[,1]
varcomp
varcomp/sum(varcomp)

#_____________________________________________________________________________
# [2] Animal model with pedigreemm
#_____________________________________________________________________________

# [2a] The animal model
#______________________________

#md$Animal = paste("AA", md$Animal, sep="")
str(md)
# Removal of intercept
mod = pedigreemm(NeophiliaScore ~ ObjectType -1 + (1|Animal) + (1|IndID), pedigree=list(Animal=pedS4), data=md)
summary(mod)

varcomp = data.frame(summary(mod)$varcor)[,4]
names(varcomp) = data.frame(summary(mod)$varcor)[,1]
varcomp/sum(varcomp)

# [2b] The animal model with additional random effects
#______________________________________________________

# Adding maternal effects
mod = pedigreemm(NeophiliaScore ~ ObjectType -1 + (1|Animal) + (1|MotherID) + (1|IndID),
                 pedigree=list(Animal=pedS4), data=md)
summary(mod)
varcomp = data.frame(summary(mod)$varcor)[,4]
names(varcomp) = data.frame(summary(mod)$varcor)[,1]
varcomp/sum(varcomp)

# Adding further early life effects
mod = pedigreemm(NeophiliaScore ~ ObjectType -1 + (1|Animal) + (1|MotherID) + (1|FosterID) + (1|PeerGroupID) + (1|IndID),
                 pedigree=list(Animal=pedS4), data=md)
summary(mod)
varcomp = data.frame(summary(mod)$varcor)[,4]
names(varcomp) = data.frame(summary(mod)$varcor)[,1]
varcomp/sum(varcomp)


#_____________________________________________________________________________
# [3] Animal model with MCMCglmm
#_____________________________________________________________________________

phenvar = var(md$NeophiliaScore)

# The repeatability model
prior = list(G = list(G1 = list(V = phenvar/2, nu = 0.002)),
             R = list(V = phenvar/2, nu = 0.002))
mod = MCMCglmm(NeophiliaScore ~ ObjectType - 1, random = ~ IndID, data=md, 
               nitt=15000, thin=5, burnin=5000, 
               start = list(QUASI=FALSE),
               prior=prior)
plot(mod$Sol)
plot(mod$VCV)

# The animal model
md$animal = md$Animal  # MCMC expects a column "animal" (with small caps) for linking to the pedigree
prior = list(G = list(G1 = list(V = phenvar/3, nu = 0.002),
                      G2 = list(V = phenvar/3, nu = 0.002)),
             R = list(V = phenvar/3, nu = 0.002))
mod = MCMCglmm(NeophiliaScore ~ ObjectType - 1, random = ~ animal + IndID, 
               data=md, pedigree = ped,
               nitt=15000, thin=5, burnin=5000, 
               start = list(QUASI=FALSE),
               prior=prior)
plot(mod$Sol)
plot(mod$VCV)

# The animal model
md$animal = md$Animal  # MCMC expects a column "animal" (with small caps) for linking to the pedigree
prior = list(G = list(G1 = list(V = phenvar/6, nu = 0.002),
                      G2 = list(V = phenvar/6, nu = 0.002),
                      G3 = list(V = phenvar/6, nu = 0.002),
                      G4 = list(V = phenvar/6, nu = 0.002),
                      G5 = list(V = phenvar/6, nu = 0.002)),
             R = list(V = phenvar/6, nu = 0.002))
mod.chain1 = MCMCglmm(NeophiliaScore ~ ObjectType - 1, random = ~ animal + MotherID + FosterID + PeerGroupID + IndID, 
                      data=md, pedigree = ped,
                      nitt=15000, thin=5, burnin=5000, 
                      start = list(QUASI=FALSE),
                      prior=prior)
mod.chain2 = MCMCglmm(NeophiliaScore ~ ObjectType - 1, random = ~ animal + MotherID + FosterID + PeerGroupID + IndID, 
                      data=md, pedigree = ped,
                      nitt=15000, thin=5, burnin=5000, 
                      start = list(QUASI=FALSE),
                      prior=prior)

# Some post processing
mcmcmod = function(Sol, VCV) {
        mychain = mcmc(data=cbind(Sol, 
                                  Vp = VCV[,"animal"] + VCV[,"MotherID"] + VCV[,"FosterID"] + VCV[,"PeerGroupID"] + VCV[,"IndID"] + VCV[,"units"]),
                       start=attr(VCV, "mcpar")[1],
                       end=attr(VCV, "mcpar")[2],
                       thin=attr(VCV, "mcpar")[3])
        mychain = mcmc(data=cbind(mychain, 
                                  AnimalComp = VCV[,"animal"]   / mychain[,"Vp"],
                                  MotherComp = VCV[,"MotherID"] / mychain[,"Vp"],
                                  FosterComp = VCV[,"FosterID"] / mychain[,"Vp"],
                                  PeerComp   = VCV[,"PeerGroupID"] / mychain[,"Vp"],
                                  IndComp    = VCV[,"IndID"]    / mychain[,"Vp"],
                                  ResidComp  = VCV[,"units"]    / mychain[,"Vp"]),
                       start=attr(VCV, "mcpar")[1],
                       end=attr(VCV, "mcpar")[2],
                       thin=attr(VCV, "mcpar")[3])
}
chain1 = mcmcmod(mod.chain1$Sol, mod.chain1$VCV)
chain2 = mcmcmod(mod.chain2$Sol, mod.chain2$VCV)
chains = mcmc.list(chain1, chain2)
plot(chains, ask=TRUE)
posterior.mode(chains)
HPDinterval(chains)
summary(chains)

#_____________________________________________________________________________
# [3] Phylogenetic models with MCMCglmm (see example ?MCMCglmm)
#_____________________________________________________________________________

?MCMCglmm
data(bird.families) 
bird.families
str(bird.families)

phylo.effect = rbv(bird.families, 1, nodes="TIPS") 
phenotype = phylo.effect+rnorm(dim(phylo.effect)[1], 0, 1)  
# simulate phylogenetic and residual effects with unit variance
test.data = data.frame(phenotype=phenotype, taxon=row.names(phenotype))
Ainv = inverseA(bird.families)$Ainv

# inverse matrix of shared phyloegnetic history
prior = list(G=list(G1=list(V=1, nu=0.002)),
             R=list(V=1, nu=0.002))
mod = MCMCglmm(phenotype ~ 1, random = ~ taxon, ginverse = list(taxon=Ainv),
               data=test.data, 
               prior=prior)

plot(mod$Sol)
plot(mod$VCV)
