
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

nsites <- seq(20, 3000, by =10) # this is the goal

nsites <- c(30,50,100,300, 500,1000,3000) # to test (13/11/2023)


# Number of visits (J): 3
J <- 3
# Number of species (nspecies): 2
nspecies <- 2
# Number of simulations (nsim) per sample size in N: 100
nsim <- 100
# Detection probabilities are equal for both spp & kept constant across all scenarios
p_true <- c(0.5,0.5) # 50% detection probability
# Number of scenarios: 8
## Null model scenarios: 4
## Covariate model scenarios: 4
## Detection probabilities are kept constant

## Null models --------------

### occ_formulas: formulas for the 2^nspecies-1 linear predictors of occupancy & 
###               co-occurrence. Null models specified with "~1" in each linear predictor
occ_formulas = rep("~1", 2^nspecies-1) 
               
### det formulas: formulas for the nspecieslinear predictors of detection.
###               Null models specified with "~1" in each linear predictor
det_formulas = rep("~1", nspecies)

# Scenario 1: Null model with strong positive interaction  ------
# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1 <-  -0.4 # (~40%)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2 <-  -0.4 # (~40%)

# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12 <- 1  # 2.71 log-odds



scen1 <- MSOM_SimandFit.fun(beta1 = beta1 , beta2, beta12, 
                            nsites = nsites, occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                           nspecies = nspecies, seed = 1337, nsim = nsim, J = 3)


saveRDS(object = scen1,file = "scen1_seed1337_simResults.rds")


# Scenario 2: Null model with weak positive interaction  ------
# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1 <-  -0.4 # (~40%)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2 <-  -0.4 # (~40%)

# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12 <- 0.2 # increase log-odds by ~22%



scen2 <- MSOM_SimandFit.fun(beta1 = beta1 , beta2, beta12, 
                            nsites = nsites, occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, seed = 1337, nsim = nsim, J = 3)




psi.fun(f= cbind(beta1,  beta2, beta12), nspecies = 2) #

saveRDS(object = scen2,file = "scen2_seed1337_simResults.rds")

# Scenario 3: Null model with weak positive interaction  ------
# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1 <-  0 # (50%)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2 <-  0 # (50%)

# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12 <- -1 # decrease log-odds by 63%



scen3 <- MSOM_SimandFit.fun(beta1 = beta1 , beta2, beta12, 
                            nsites = nsites, occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, seed = 1337, nsim = nsim, J = 3)


saveRDS(object = scen3,file = "scen3_seed1337_simResults.rds")

# Scenario 4: Null model with weak positive interaction  ------
# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1 <-  0 # (~50%)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2 <-  0 # (~50%)

# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12 <- -0.2 # decrease log-odds by 18%



scen4 <- MSOM_SimandFit.fun(beta1 = beta1 , beta2, beta12, 
                            nsites = nsites, occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, seed = 1337, nsim = nsim, J = 3)



saveRDS(object = scen4,file = "scen4_seed1337_simResults.rds")


## Covariate models --------------

# Consistent through all scenarios:

### det formulas: formulas for the nspecieslinear predictors of detection.
###               Null models specified with "~1" in each linear predictor
det_formulas = rep("~1", nspecies)

## nocc_covs: number of occupancy (and co-occurrence covariates) to be generated
##            from standard normal distributions occ_cov~N(0,1) and to include in
##            the linear predictors for occupancy

nocc_covs <- 4
# General formulas to be:
# f1 = beta1.0 + beta1.1*occ_cov1 + beta1.2*occ_cov2
# f2 = beta2.0 + beta2.1*occ_cov2
# f12 = beta12.0 + beta12.1*occ_cov2 + beta12.2*occ_cov4


occ_formulas <- c("~occ_cov1+occ_cov2",     #f1
                  "~occ_cov3",              #f2
                  "~occ_cov2 + occ_cov4")   #f12

# Scenario 5: Covariate model with strong positive baseline interaction  ------

# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1.0 <-  -0.4 # baseline prob(~40%)
beta1.1 <- -0.2 # occ_cov1 effect 
beta1.2 <- 0.4 # occ_cov2 effect

beta1 <- c(beta1.0, beta1.1, beta1.2)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2.0 <-  -0.6 # (~35%)
beta2.1 <- 0.5

beta2 <- c(beta2.0, beta2.1)
# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12.0 <- 0.7  # 2.71 log-odds
beta12.1 <- -0.3
beta12.2 <- 0.2

beta12 <- c(beta12.0, beta12.1, beta12.2)


scen5 <- MSOM_SimandFit.fun(beta1 = beta1 , beta2, beta12, 
                            nsites = 3000, occ_formulas = occ_formulas,
                            det_formulas = det_formulas, nocc_covs = nocc_covs,
                            nspecies = nspecies, seed = 1337, nsim = 2, J = 3)


# Scenario 6: Covariate model with weak positive baseline interaction  ------

# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1.0 <-  -0.4 # baseline prob(~40%)
beta1.1 <- -0.2 # occ_cov1 effect 
beta1.2 <- 0.4 # occ_cov2 effect

beta1 <- c(beta1.0, beta1.1, beta1.2)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2.0 <-  -0.6 # (~35%)
beta2.1 <- 0.5

beta2 <- c(beta2.0, beta2.1)
# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12.0 <- 0.3  # ~35% increased log-odds
beta12.1 <- -0.3
beta12.2 <- 0.2

beta12 <- c(beta12.0, beta12.1, beta12.2)


scen6 <- MSOM_SimandFit.fun(beta1 = beta1 , beta2, beta12, 
                            nsites = 3000, occ_formulas = occ_formulas,
                            det_formulas = det_formulas, nocc_covs = nocc_covs,
                            nspecies = nspecies, seed = 1337, nsim = 2, J = 3)


# Scenario 7: Covariate model with strong negative baseline interaction  ------

# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1.0 <-  0.2 # baseline prob(~40%)
beta1.1 <- -0.2 # occ_cov1 effect 
beta1.2 <- 0.4 # occ_cov2 effect

beta1 <- c(beta1.0, beta1.1, beta1.2)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2.0 <- -0.2 # (~35%)
beta2.1 <- 0.5

beta2 <- c(beta2.0, beta2.1)
# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12.0 <- -0.7  # -50% baseline log-odds
beta12.1 <- -0.3
beta12.2 <- 0.4

beta12 <- c(beta12.0, beta12.1, beta12.2)

scen7 <- MSOM_SimandFit.fun(beta1 = beta1 , beta2, beta12, 
                            nsites = 3000, occ_formulas = occ_formulas,
                            det_formulas = det_formulas, nocc_covs = nocc_covs,
                            nspecies = nspecies, seed = 1337, nsim = 2, J = 3)






# Scenario 8: Covariate model with weak negative baseline interaction  ------

# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1.0 <-  0.2 # baseline prob(~40%)
beta1.1 <- -0.2 # occ_cov1 effect 
beta1.2 <- 0.4 # occ_cov2 effect

beta1 <- c(beta1.0, beta1.1, beta1.2)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2.0 <- -0.2 # (~35%)
beta2.1 <- 0.5

beta2 <- c(beta2.0, beta2.1)
# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12.0 <- -0.3  # baseline log-odds decrease by ~26%
beta12.1 <- -0.3
beta12.2 <- 0.4

beta12 <- c(beta12.0, beta12.1, beta12.2)

scen8 <- MSOM_SimandFit.fun(beta1 = beta1 , beta2, beta12, 
                            nsites = 3000, occ_formulas = occ_formulas,
                            det_formulas = det_formulas, nocc_covs = nocc_covs,
                            nspecies = nspecies, seed = 1337, nsim = 2, J = 3)







