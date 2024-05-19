
############################################################
##    Simulation study to evaluate species co-occurrence  ##
##  model performance under different sample size regimes ##
############################################################

### Set up libary
library(dplyr, quietly = T)
library(unmarked, quietly = T)
library(ggplot2, quietly = T)

# Source simulation function
source("MSOM_SimFun.R")


# General characteristics (constant across scenarios) ------------
# Sample size scenarios (nsites or N):  30 to 3,000 in steps of 10

#nsites <- seq(20, 3000, by =10) # this is the goal

#nsites <- c(30,50,100,300, 500,1000,3000) # to test (13/11/2023)

ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)


# Number of simulations (nsim) per sample size in N: 100
nsim <- 100

# Number of visits (J): 3
J <- 3

## Covariate models --------------


# Consistent through all scenarios:
nspecies <- 2
p_true <- rep(0.5, nspecies)
### det formulas: formulas for the nspecieslinear predictors of detection.
###               Null models specified with "~1" in each linear predictor
det_formulas <- rep("~1", nspecies)

# Scenario cov4: Covariate model shared effect in interaction and f1  ------

# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1.0 <-  -0.2 # baseline prob(~40%)
beta1.1 <- 0.4 # occ_cov1 effect 
beta1.2 <- -0.3 # 

beta1 <- c(beta1.0, beta1.1, beta1.2)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2.0 <-  -0.4 # (~35%)
beta2.1 <- 0.6

beta2 <- c(beta2.0, beta2.1)
# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2

beta12.0 <- -0.6  
beta12.1 <- 0.8
beta12.2 <- -0.5

beta12 <- c(beta12.0, beta12.1, beta12.2)

beta = list(beta1, beta2, beta12)


occ_formulas <- c("~occ_cov1 + occ_cov3",          #f1
                  "~occ_cov2",          #f2
                  "~occ_cov1 + occ_cov4")   #f12

nocc_covs = 4

seed = 1337

scencov.4 <-  MSOM_simfit.fun.v2(beta = beta, 
                                nsites = nsites, 
                                occ_formulas = occ_formulas,
                                det_formulas = det_formulas, 
                                nocc_covs = nocc_covs,
                                nspecies = nspecies, seed = seed, 
                                nsim = nsim, J = J)

scencov.4$State.params.pl

scenario <- "Scenario4_2sp"
model.type = "NCovariate4_SharedCovs"


saveRDS(object = scencov.4,file = paste0(scenario, model.type, "seed", seed,
                                         "_simResults_v2.rds" ))


