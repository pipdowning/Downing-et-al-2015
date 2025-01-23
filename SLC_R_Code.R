# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #
#                                                                                                             #
# Supplementary R code for:                                                                                   #
#                                                                                                             #
#                 SEX, LONG LIFE AND THE EVOLUTIONARY TRANSITION TO COOPERATIVE BREEDING IN BIRDS             #
#                                                                                                             #
# code written by P.A.D. and C.K.C. (philip.downing@zoo.ox.ac.uk)                                             #
#                                                                                                             #
#                                                                                                             #
# NOTE - each of the analyses below is repeated 10 times, each time using a different MCC phylogenetic tree   #
#        posterior distributions from the 10 repeats of each model were combined for parameter estimation     #
#        we included B priors to improve mixing of fixed effects in models with categorical responses         #
#        the raw data in supplementary table 2 was manipulated prior to analysis to get one value of each     #
#        parameter (survival and promiscuity) per species                                                     #
#                                                                                                             #
# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #

# packages

library(ape)
library(coda)
library(MCMCglmm)
library(QuantPsyc)

# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

### DATA

## Phenotypic Data
sxData <- read.csv("")      # dataset with 238 species.
# 6 columns: species, cooperation (2 level factor), latitude (continuous), mass (continuous), promiscuity (%), survival (%)
# these variables are transformed in various ways - see 'Methods: Analyses' for details
breedData <- read.csv("")   # dataset with 40 species
# 6 columns: species, cooperation (2 level factor), age (discrete), mass (continuous)
# most species appear in the dataset multiple times, depending the number of years for which reproduction is delayed (see supplementary tables)
reproData <- read.csv("")   # dataset with 40 species
# 3 columns: species, the number of young sired by dominant individuals within groups, the number of young sired by subordinate individuals within groups

## Trees
sxTree <- read.nexus("")       # 1 of 10 phylogenetic trees with 238 species matching those in SxData
breedTree <- read.nexus("")    # 1 of 10 phylogenetic trees with 40 species matching those in breedData
reproTree <- read.nexus("")    # 1 of 10 phylogenetic trees with 40 species matching those in reproData

# note that the species in 'breedData' and 'reproData' are not necessarily included in sxData
# only species for which both a measure of survival and promiscuity are included in 'sxData'
# we include priors for fixed effects (B) to improve mixing of fixed effects in models with categorical response variables

# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

## PHYLOGENETIC HERITABILITY

# survival - phylogenetic heritability
pr.PH.1 <- list(R = list(V=1, nu=0.002), G = list(G1=list(V=1, nu=0.002)))
PHSx <- MCMCglmm(cbind(alive, dead) ~ 1, random= ~animal, pedigree=sxTree, family="multinomial2", nodes="ALL", data=sxData, prior=pr.PH.1, verbose=FALSE, nitt=4100000, burnin=100000, thin=1000, pr=TRUE, slice=TRUE)
posterior.mode(PHSx$VCV[,1] / (PHSx$VCV[,1] + PHSx$VCV[,2]))
HPDinterval(PHSx$VCV[,1] / (PHSx$VCV[,1] + PHSx$VCV[,2]))

# promiscuity - phylogenetic heritability
pr.PH.2 <- list(R = list(V=1, nu=0.002), G = list(G1=list(V=1, nu=0.002)))
PHprom <- MCMCglmm(cbind(extraGroup, withinGroup) ~ 1, random= ~animal, pedigree=sxTree, family="multinomial2", nodes="ALL", data=sxData, prior=pr.PH.2, verbose=FALSE, nitt=4100000, burnin=100000, thin=1000, pr=TRUE, slice=TRUE)
posterior.mode(PHprom$VCV[,1] / (PHprom$VCV[,1] + PHprom$VCV[,2]))
HPDinterval(PHprom$VCV[,1] / (PHprom$VCV[,1] + PHprom$VCV[,2]))

# cooperation - intra-class correlation coefficient
pr.PH.3 <- list(B = list(mu=0, V=1*(1+pi^2/3)), R = list(V=1, fix=1), G = list(G1=list(V=1, nu=0.002)))
PHcoop <- MCMCglmm(coop ~ 1, random= ~animal, pedigree=sxTree, family="categorical", nodes="ALL", data=sxData, prior=pr.PH.3, verbose=FALSE, nitt=4100000, burnin=100000, thin=1000, pr=TRUE)
posterior.mode(PHcoop$VCV[,1] / ((PHcoop$VCV[,1] + PHcoop$VCV[,2])+pi^2/3))
HPDinterval(PHcoop$VCV[,1] / ((PHcoop$VCV[,1] + PHcoop$VCV[,2])+pi^2/3))

# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------

## STATISTICS

# QUESTION 1: IS REPRODUCTION DELAYED IN COOPERATIVE BREEDERS?

# A. Mean Age At First Breeding
pr.Q1.A <- list(R = list(V=1, nu=0.002), G = list(G1=list(V=1, nu=0.002)))
mod1A <- MCMCglmm(meanFirstLogZ ~ mass + cooperation, random= ~animal, family = "gaussian", pedigree=breedTree, data=meanFirst, prior=pr.Q1.A, verbose=FALSE, nitt=4100000, burnin=100000, thin=1000, pr=TRUE, slice=TRUE)
# the data used to run this model are calculated from breedData so that each species has just one entry in the dataframe

# B. First Breeding Through Time
pr.Q1.B <- list(R = list(V=1, nu=0.002), G = list(G1=list(V=diag(2), nu=0.002), G2=list(V=diag(2), nu=0.002)))
mod1B <- MCMCglmm(cbind(breeds, notbreeds) ~ mass + age*cooperation, random= ~us(1 + age):animal + us(1 + age):ID, family = "multinomial2", pedigree=breedTree, data=breedData, prior=pr.Q1.B, verbose=FALSE, nitt=4100000, burnin=100000, thin=1000, pr=TRUE, slice=TRUE)

# C. Direct Fitness Benefits for Helpers
pr.Q1.C <- list(R = list(V=1, nu=0.002), G = list(G1=list(V=1, nu=0.002)))
mod1C <- MCMCglmm(cbind(dominant, subordinate) ~ 1, random= ~animal, pedigree=reproTree, family="multinomial2", nodes="ALL", data=reproData, prior=pr.Q1.C, verbose=FALSE, nitt=4100000, burnin=100000, thin=1000, pr=TRUE, slice=TRUE)
mcmc(inv.logit(diffMod$Sol[,1]))
# to transform coefficient
posterior.mode(mcmc(inv.logit(mod1C$Sol[,1])))
HPDinterval(mcmc(inv.logit(mod1C$Sol[,1])))

# -----------------------------------------------------------------------------------------------------------

# QUESTION 2: DO COOPERATIVE BREEDERS LIVE LONGER THAN NON-COOPERATIVE BREEDERS?

pr.Q2 <- list(R = list(V=1, nu=0.002), G = list(G1=list(V=1, nu=0.002)))
mod2 <- MCMCglmm(cbind(alive, dead) ~ mass.gramsT + mean.latitudeT  + promT + coop, random=~animal, pedigree=sxTree, family ="multinomial2", nodes="ALL", data = sxData, prior=pr.Q2, verbose=FALSE, nitt=4100000, burnin=100000, thin=1000, pr=TRUE, slice=TRUE)

# -----------------------------------------------------------------------------------------------------------

# QUESTION 3: DOES HIGH SURVIVAL MAKE THE EVOLUTION OF COOPERATIVE BREEDING MORE LIKELY?

# A. Phylogenetic Correlation
pr.Q3.A <- list(R = list(V=diag(2), nu=1.002, fix=2), G = list(G1=list(V=diag(2), nu=1.002)))
mod3A <- MCMCglmm(cbind(cbind(alive, dead), coop) ~ trait - 1,
                                    random = ~ us(trait):animal, rcov = ~ us(trait):units,
                                    family = c("multinomial2", "categorical"), data = sxData, pedigree = sxTree,
                                    prior = pr.Q3.A, verbose = FALSE, nitt = 4100000, burnin = 100000, thin = 1000, pr = TRUE, slice = TRUE)

phylogenetic.corr3A <- mod3A$VCV[,2] / sqrt(mod3A$VCV[,1] * mod3A$VCV[,4])
posterior.mode(phylogenetic.corr3A)
HPDinterval(phylogenetic.corr3A)

# B. Causality - reconstruct ancestral states
pr.Q3.B <- list(B = list(mu=0, V=1*(1+pi^2/3)), R = list(V=1, fix=1), G = list(G1=list(V=1, nu=0.002)))
recon <- MCMCglmm(coop ~ 1, random = ~ animal, pedigree = sxTree, family = "categorical", nodes="ALL", data= sxData, prior=pr.Q3.B, verbose=FALSE, nitt=4100000, burnin=400000, thin=1000, pr=TRUE)

sxData$nodecode <- paste("animal", sxData$animal,sep=".")

# create transition dataset
Coop <- data.frame(animal = colnames(recon$Sol), coop = posterior.mode(mcmc(inv.logit(recon$Sol[sample(500),]))))
Coop$coop2[Coop$coop > 0.95] <- "Cooperative"
Coop$coop2[Coop$coop < 0.05] <- "Noncooperative"
Coop$coop2[Coop$coop < 0.95 & Coop$coop > 0.05] <- NA

TransDat <- as.data.frame(sxTree$edge)
TransDat$V1name <- paste("animal.Node", (TransDat$V1-length(sxTree$tip.label)), sep="")
TransDat$V2name <- paste("animal.Node", (TransDat$V2-length(sxTree$tip.label)), sep="")
treesp <- data.frame(tip.label=paste("animal.", sxTree$tip.label,sep=""), no=1:length(sxTree$tip.label))
TransDat <- data.frame(TransDat, ancestors=treesp$tip.label[match(TransDat$V1, treesp$no)], descendents=treesp$tip.label[match(TransDat$V2, treesp$no)])
TransDat$ancestors <- as.character(TransDat$ancestors)
TransDat$descendents<-as.character(TransDat$descendents)
TransDat$ancestors <- as.character(ifelse(!is.na(TransDat$ancestors),TransDat$ancestors,TransDat$V1name))
TransDat$descendents <- as.character(ifelse(!is.na(TransDat$descendents),TransDat$descendents,TransDat$V2name))
TransDat <- data.frame(TransDat, ancCoop=Coop$coop2[match(TransDat$ancestors,Coop$animal)], desCoop=Coop$coop2[match(TransDat$descendents,Coop$animal)])

obs <- data.frame(table(TransDat$desCoop,TransDat$ancestors))
obs$CAT <- obs$Freq
obs$CAT[obs$Freq == 2 & obs$Var1 == "Cooperative"] <- "Only Cooperative"
obs$CAT[obs$Freq == 0 & obs$Var1 == "Noncooperative"] <- "Only Cooperative"
obs$CAT[obs$Freq == 2 & obs$Var1 == "Noncooperative"] <- "Only Noncooperative"
obs$CAT[obs$Freq == 0 & obs$Var1 == "Cooperative"] <- "Only Noncooperative"
obs$CAT[obs$Freq == 1] <- "Both"

TransDat <- TransDat[!is.na(TransDat$ancCoop),]
TransDat <- TransDat[!is.na(TransDat$desCoop),]
TransDat <- data.frame(TransDat, CAT = obs$CAT[match(TransDat$ancestors,obs$Var2)])
TransDat$CAT2 <- paste(TransDat$ancCoop, TransDat$CAT, sep = ".")
TransDat$CAT2[which(TransDat$CAT2 == "Noncooperative.Only Cooperative")] <- "Noncooperative.Both"
TransDat$CAT2[which(TransDat$CAT2 == "Cooperative.Only Noncooperative")] <- "Cooperative.Both"

TransDat <- data.frame(TransDat, trait1 = sxData$alive[match(TransDat$descendents, sxData$nodecode)], trait2 = sxData$dead[match(TransDat$descendents, sxData$nodecode)])
tree <- makeNodeLabel(sxTree, method = "number")
TransDat$animal <- gsub("animal.", "", TransDat$ancestors)
TransDat <- TransDat[TransDat$animal !='Node1',]
TransDat$CAT2 <- as.factor(TransDat$CAT2)
TransDat$CAT2 <- relevel(TransDat$CAT2, ref = "Noncooperative.Both")
rownames(TransDat) <- NULL

pr.Q3.B.1 <- list(B = list(mu=rep(0, length(table(TransDat$CAT2))), V=diag(length(table(TransDat$CAT2)))*(1+pi^2/3)), R = list(V=1,nu=0.002), G = list(G1=list(V=1,nu=0.002)))
M1 <- MCMCglmm(cbind(trait1, trait2) ~ CAT2, random = ~ animal, pedigree=tree, family = "multinomial2", nodes="ALL", data = TransDat, prior=pr.Q3.B.1, nitt=510000, burnin=10000, thin=500, pr=TRUE,verbose = FALSE)

TransSolRes <- data.frame(matrix(data=0,nrow=0,ncol=length(M1$Sol[1,])))
colnames(TransSolRes) <- colnames(M1$Sol[,order(colnames(M1$Sol))])

TransVCVRes <- data.frame(matrix(data=0,nrow=0,ncol=length(M1$VCV[1,])))
colnames(TransVCVRes) <- colnames(M1$VCV)

## sample 500 times from the posterior distribution of nodal estimates. repeat 100 times 

for(i in 1:100) {
  
  Coop <- data.frame(animal = colnames(recon$Sol), coop = posterior.mode(mcmc(inv.logit(recon$Sol[sample(500),]))))
  Coop$coop2[Coop$coop > 0.95] <- "Cooperative"
  Coop$coop2[Coop$coop < 0.05] <- "Noncooperative"
  Coop$coop2[Coop$coop < 0.95 & Coop$coop > 0.05] <- NA
  
  TransDat <- as.data.frame(sxTree$edge)
  TransDat$V1name <- paste("animal.Node", (TransDat$V1-length(sxTree$tip.label)), sep="")
  TransDat$V2name <- paste("animal.Node", (TransDat$V2-length(sxTree$tip.label)), sep="")
  treesp <- data.frame(tip.label=paste("animal.", sxTree$tip.label,sep=""), no=1:length(sxTree$tip.label))
  TransDat <- data.frame(TransDat, ancestors=treesp$tip.label[match(TransDat$V1, treesp$no)], descendents=treesp$tip.label[match(TransDat$V2, treesp$no)])
  TransDat$ancestors <- as.character(TransDat$ancestors)
  TransDat$descendents<-as.character(TransDat$descendents)
  TransDat$ancestors <- as.character(ifelse(!is.na(TransDat$ancestors),TransDat$ancestors,TransDat$V1name))
  TransDat$descendents <- as.character(ifelse(!is.na(TransDat$descendents),TransDat$descendents,TransDat$V2name))
  TransDat <- data.frame(TransDat, ancCoop=Coop$coop2[match(TransDat$ancestors,Coop$animal)], desCoop=Coop$coop2[match(TransDat$descendents,Coop$animal)])
  
  obs <- data.frame(table(TransDat$desCoop,TransDat$ancestors))
  obs$CAT <- obs$Freq
  obs$CAT[obs$Freq == 2 & obs$Var1 == "Cooperative"] <- "Only Cooperative"
  obs$CAT[obs$Freq == 0 & obs$Var1 == "Noncooperative"] <- "Only Cooperative"
  obs$CAT[obs$Freq == 2 & obs$Var1 == "Noncooperative"] <- "Only Noncooperative"
  obs$CAT[obs$Freq == 0 & obs$Var1 == "Cooperative"] <- "Only Noncooperative"
  obs$CAT[obs$Freq == 1] <- "Both"
  
  TransDat <- TransDat[!is.na(TransDat$ancCoop),]
  TransDat <- TransDat[!is.na(TransDat$desCoop),]
  TransDat <- data.frame(TransDat, CAT = obs$CAT[match(TransDat$ancestors,obs$Var2)])
  TransDat$CAT2 <- paste(TransDat$ancCoop, TransDat$CAT, sep = ".")
  TransDat$CAT2[which(TransDat$CAT2 == "Noncooperative.Only Cooperative")] <- "Noncooperative.Both"
  TransDat$CAT2[which(TransDat$CAT2 == "Cooperative.Only Noncooperative")] <- "Cooperative.Both"
  
  TransDat <- data.frame(TransDat, trait1 = sxData$alive[match(TransDat$descendents, sxData$nodecode)], trait2 = sxData$dead[match(TransDat$descendents, sxData$nodecode)])
  tree <- makeNodeLabel(sxTree, method = "number")
  TransDat$animal <- gsub("animal.", "", TransDat$ancestors)
  TransDat <- TransDat[TransDat$animal !='Node1',]
  TransDat$CAT2 <- as.factor(TransDat$CAT2)
  TransDat$CAT2 <- relevel(TransDat$CAT2, ref = "Noncooperative.Both")
  rownames(TransDat) <- NULL
  
  pr.Q3.B.2 <- list(B = list(mu=rep(0, length(table(TransDat$CAT2))), V=diag(length(table(TransDat$CAT2)))*(1+pi^2/3)), R = list(V=1,nu=0.002), G = list(G1=list(V=1,nu=0.002)))
  M2 <- MCMCglmm(cbind(trait1, trait2) ~ CAT2, random= ~animal, pedigree=tree, family = "multinomial2", nodes="ALL", data = TransDat, prior=pr.Q3.B.2, nitt = 5100000, burnin = 100000, thin = 5000, pr=TRUE,verbose = FALSE)
  
  TransSolRes <- data.frame(rbind(as.mcmc(TransSolRes), M2$Sol[,order(colnames(M2$Sol))]))
  TransVCVRes <- data.frame(rbind(as.mcmc(TransVCVRes), M2$VCV))
} 

results <- list(TransSolRes, TransVCVRes)

# -----------------------------------------------------------------------------------------------------------

# QUESTION 4: IS LONG LIFE MORE PRONOUNCED IN PROMISCUOUS COOPERATIVE BREEDERS?

# A. Does the Relationship between Promiscuity and Survival Differ between Cooperative and Non-Cooperative Species
pr.Q4.A <- list(R = list(V=1, nu=0.002), G = list(G1=list(V=1, nu=0.002)))
mod4A <- MCMCglmm(cbind(extraGroup, withinGroup) ~ (mass.gramsT + mean.survivalT) * coop, random= ~animal, pedigree=sxTree, family ="multinomial2", nodes="ALL", data=sxData, prior=pr.Q4.A, verbose=FALSE, nitt=4100000, burnin=100000, thin=1000, pr=TRUE, slice=TRUE)

# B. Phylogenetic Correlation between Promiscuity and Survival for Different Levels of Cooperation
pr.Q4.B <- list(R = list(V=diag(2), nu=1.002), G = list(G1=list(V=diag(2), nu=1.002), G2=list(V=diag(2), nu=1.002)))
mod4B <- MCMCglmm(cbind(cbind(alive, dead), cbind(extraGroup, withinGroup)) ~ at.level(coop,"1"):trait+at.level(coop,"0"):trait-1,
                          random = ~us(at.level(coop,"1"):trait):animal+us(at.level(coop,"0"):trait):animal, rcov = ~us(trait):units,
                          family = c("multinomial2", "multinomial2"), data = sxData, pedigree = sxTree,
                          prior = pr.Q4.B, verbose = FALSE, nitt = 4100000, burnin = 100000, thin = 1000, pr = TRUE, slice = TRUE)

phylogenetic.corr4Bcoop <- mod4B$VCV[,2] / sqrt(mod4B$VCV[,1] * mod4B$VCV[,4])
posterior.mode(phylogenetic.corr4Bcoop)
HPDinterval(phylogenetic.corr4Bcoop)
phylogenetic.corr4Bnoncoop <- mod4B$VCV[,6] / sqrt(mod4B$VCV[,5] * mod4B$VCV[,8])
posterior.mode(phylogenetic.corr4Bnoncoop)
HPDinterval(phylogenetic.corr4Bnoncoop)

# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #
#                                                                                                             #
#                                             END - thanks for reading! p.                                    #
#                                                                                                             #
# ----------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #