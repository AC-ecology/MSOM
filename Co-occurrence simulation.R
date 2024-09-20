
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

ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)

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


saveRDS(object = scen1,file = "scen1_seed1337_simResults2.rds")


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





saveRDS(object = scen2,file = "scen2_seed1337_simResults2.rds")

# Scenario 3: Null model with weak positive interaction  ------
# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1 <-  0.1 # (~52%)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2 <-  0.1 # (~52%)

# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12 <- -1 # decrease log-odds by 63%



scen3 <- MSOM_SimandFit.fun(beta1 = beta1 , beta2, beta12, 
                            nsites = nsites, occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, seed = 1337, nsim = nsim, J = 3)


saveRDS(object = scen3,file = "scen3_seed1337_simResults2.rds")

# Scenario 4: Null model with weak positive interaction  ------
# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1 <-  0.1 # (~50%)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2 <-  0.1 # (~50%)

# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12 <- -0.2 # decrease log-odds by 18%



scen4 <- MSOM_SimandFit.fun(beta1 = beta1 , beta2, beta12, 
                            nsites = nsites, occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, seed = 1337, nsim = nsim, J = 3)



saveRDS(object = scen4,file = "scen4_seed1337_simResults2.rds")








