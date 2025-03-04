
############################################################
##    Simulation study to evaluate species co-occurrence  ##
##  model performance under different sample size regimes ##
############################################################

### Set up libary
library(dplyr, quietly = T)
library(unmarked, quietly = T)
library(ggplot2, quietly = T)

# Source simulation functions
source("MSOM_SimFun_Revised.R")


# General characteristics (constant across scenarios) ------------
# Sample size scenarios (nsites or N):  30 to 3,000 in steps of 10
ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)

# Number of simulations (nsim) per sample size in N: 100
nsim <- 1000

# Number of visits (J): 3
J <- 3
### All parameters ----
## First order
beta1 <- -0.4  #f1
beta2 <- 0.3   #f2
beta3 <- 0.8   #f3
beta4 <- -0.2  #f4
beta5 <- 0.5   #f5

## Second order
beta12 <- 1    #f12
beta13 <- -0.6 #f13
beta14 <- 0.8  #f14
beta15 <- 0.6  #f15
beta23 <- -1   #f23
beta24 <- 0.6  #f24
beta25 <- -1   #f25  
beta34 <- -0.8 #f34
beta35 <- 0.8  #f35
beta45 <- -0.6 #f45




# Network Models ---------
## 3 species simulation ------

ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)


# Number of simulations (nsim) per sample size in N: 100
nsim <- 1000

# Number of visits (J): 3
J <- 3

nspecies <- 3
p_true <- rep(0.5, nspecies) # 50% detection probability

occ_formulas <-  rep("~1", nspecies + nspecies*(nspecies-1)/2)
det_formulas <- rep("~1", nspecies)

beta <-  list(beta1, beta2, beta3, 
              beta12, beta13, beta23)

# Check marginal occupancy
og.psi <- gen.param.fun(beta, nspecies = nspecies)
Gen.Par <- colnames(og.psi)
c(sum(og.psi[substr(Gen.Par,1,1) ==1]),
  sum(og.psi[substr(Gen.Par,2,2) ==1]),
  sum(og.psi[substr(Gen.Par,3,3) ==1]))

## Run
seed <- 1337
scen.3sp <- MSOM_simfit.fun.v2(beta = beta, p_true = p_true,
                               nsites = nsites, 
                               nspecies = nspecies,
                               store.data = F,
                               seed = seed, 
                               nsim = nsim, J = J)


scenario <- "Scenario1_3sp"
model.type = "Null"

saveRDS(object = scen.3sp,file = paste0("RevisionSimulations/Results/",scenario, model.type, "seed", seed,
                                         "_simResults.rds" ))

## 4 species simulation ----


ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)


# Number of simulations (nsim) per sample size in N: 1000
nsim <- 1000
# Number of visits (J): 3
J <- 3

nspecies <- 4
p_true <- rep(0.5, nspecies) # 50% detection probability

occ_formulas <-  rep("~1", nspecies + nspecies*(nspecies-1)/2)
det_formulas <- rep("~1", nspecies)

beta <-  list(beta1, beta2, beta3, beta4,
              beta12, beta13, beta14,
              beta23, beta24,
              beta34)
# 
# og.psi <- gen.param.fun(beta, nspecies = nspecies)
# Gen.Par <- colnames(og.psi)
# c(sum(og.psi[substr(Gen.Par,1,1) ==1]),
#   sum(og.psi[substr(Gen.Par,2,2) ==1]),
#   sum(og.psi[substr(Gen.Par,3,3) ==1]),
#   sum(og.psi[substr(Gen.Par,4,4) ==1]))



seed <- 1337
scen.4sp <- MSOM_simfit.fun.v2(beta = beta, 
                               p_true = p_true,
                               nsites = nsites, 
                               nspecies = nspecies,
                               store.data = F,
                               seed = seed, 
                               nsim = nsim, J = J)


scenario <- "Scenario1_4sp"
model.type = "Null"

saveRDS(object = scen.4sp,file = paste0("RevisionSimulations/Results/",scenario, model.type, "seed", seed,
                                        "_simResults.rds" ))
## 5 species ----


ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)


# Number of simulations (nsim) per sample size in N: 100
nsim <- 1000

# Number of visits (J): 3
J <- 3

nspecies <- 5
p_true <- rep(0.5, nspecies) # 50% detection probability

occ_formulas <-  rep("~1", nspecies + nspecies*(nspecies-1)/2)
det_formulas <- rep("~1", nspecies)

beta <-  list(beta1, beta2, beta3, beta4, beta5,
              beta12, beta13, beta14, beta15,
              beta23, beta24, beta25, 
              beta34, beta35,
              beta45)

## Check marginal occ.probs
og.psi <- gen.param.fun(beta, nspecies = nspecies)
Gen.Par <- colnames(og.psi)
c(sum(og.psi[substr(Gen.Par,1,1) ==1]),
  sum(og.psi[substr(Gen.Par,2,2) ==1]),
  sum(og.psi[substr(Gen.Par,3,3) ==1]),
  sum(og.psi[substr(Gen.Par,4,4) ==1]),
  sum(og.psi[substr(Gen.Par,5,5) ==1]))

## Run

seed <- 1337
scen.5sp <- MSOM_simfit.fun.v2(beta = beta, 
                               p_true = p_true,
                               nsites = nsites, 
                               nspecies = nspecies,
                               seed = seed, 
                               store.data = F,
                               nsim = nsim, J = J)


# scen.5sp$time.ellapsed

scenario <- "Scenario1_5sp"
model.type = "Null"

saveRDS(object = scen.5sp,file = paste0("RevisionSimulations/Results/",scenario, model.type, "seed", seed,
                                        "_simResults.rds" ))

