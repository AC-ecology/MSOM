
############################################################
##    Simulation study to evaluate species co-occurrence  ##
##  model performance under different sample size regimes ##
############################################################

### Set up libary
library(dplyr, quietly = T)
library(unmarked, quietly = T)
library(ggplot2, quietly = T)

# Source simulation functions
 # source("MSOM_SimFun.R")
source("MSOM_SimFun_Revised.R")

# General characteristics (constant across scenarios) ------------
# Sample size scenarios (nsites or N):  30 to 3,000 in steps of 10
ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)

# Number of simulations (nsim) per sample size in N: 1000
nsim <- 1000

# Number of visits (J): 3
J <- 3

# Scenarios with 2 species --------

# Number of species (nspecies): 2
nspecies <- 2

# Detection probabilities are equal for both spp & kept constant across all scenarios
p_true <- rep(0.5, nspecies) # 50% detection probability

# Number of scenarios: 8
## Null model scenarios: 4
## Covariate model scenarios: 4
## Detection probabilities are kept constant

# Null models --------------

### occ_formulas: formulas for the 2^nspecies-1 linear predictors of occupancy & 
###               co-occurrence. Null models specified with "~1" in each linear predictor
occ_formulas = rep("~1", 2^nspecies-1) 
nocc_covs = 3

### det formulas: formulas for the nspecieslinear predictors of detection.
###               Null models specified with "~1" in each linear predictor
det_formulas = rep("~1", nspecies)
ndet_covs = 2
## Test -----
## beta1: baseline occupancy of sp1 
beta1 <-  -0.4 # (~40%)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2 <-  -0.4 # (~40%)

# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12 <- 1  # 2.71 log-odds
beta = list(beta1, beta2, beta12)
# 
# seed = 39234
# test <- MSOM_simfit.fun.v2(beta = beta,
#                            nsites = 50,
#                            occ_formulas = occ_formulas,
#                            det_formulas = det_formulas,
#                            nspecies = nspecies,
#                            seed = seed,
#                            nsim = 30,
#                            store.data = F,
#                            J = J)
# 
# scenario <- "Test"
# model.type = "Null_HighDet"
# saveRDS(object = test,
#         file = paste0("RevisionSimulations/Results/",scenario, model.type, "seed", seed,
#                                      "NSim", nsim,
#                                      "_simResults.rds" ))

# test$State.params # For occupancy estimates with normal likelihood
## Scenario 1: Null model with strong positive interaction  ------
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

beta = list(beta1, beta2, beta12)

seed = 1337
scen1 <- MSOM_simfit.fun.v2(beta = beta, 
                            nsites = nsites, occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, 
                            seed = seed, 
                            nsim = nsim, 
                            store.data = F,
                            J = J)

scenario <- "Scenario1_StrPos"
model.type = "Null_HighDet"

saveRDS(object = scen1,file = paste0("RevisionSimulations/Results/",scenario, model.type, "seed", seed,
                                     "NSim", nsim,
                                     "_simResults.rds" ))


## Scenario 2: Null model with weak positive interaction  ------
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



beta = list(beta1, beta2, beta12)

seed = 1337
scen2 <- MSOM_simfit.fun.v2(beta = beta, 
                            nsites = nsites, occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, 
                            seed = seed, 
                            nsim = nsim, 
                            store.data = F,
                            J = J)
scenario <- "Scenario2_WkPos"
model.type = "Null_HighDet"

saveRDS(object = scen2,file = paste0("RevisionSimulations/Results/",scenario, model.type, "seed", seed,
                                     "_simResults.rds" ))

## Scenario 3: Null model with strong negative interaction  ------
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

beta = list(beta1, beta2, beta12)


scen3 <- MSOM_simfit.fun.v2(beta = beta, 
                            nsites = nsites,
                            occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, 
                            seed = seed, 
                            nsim = nsim, 
                            store.data = F,
                            J = J)


scenario <- "Scenario3_StrNeg"
model.type = "Null_HighDet"

saveRDS(object = scen3,file = paste0("RevisionSimulations/Results/",scenario, model.type, "seed", seed,
                                     "_simResults.rds" ))

## Scenario 4: Null model with weak negative interaction  ------
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

beta = list(beta1, beta2, beta12)

scen4 <- MSOM_simfit.fun.v2(beta = beta, 
                            nsites = nsites,
                            occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, 
                            seed = seed, 
                            nsim = nsim, 
                            store.data = F,
                            J = J)



scenario <- "Scenario4_WkNeg"
model.type = "Null_HighDet"

saveRDS(object = scen4,file = paste0("RevisionSimulations/Results/",scenario, model.type, "seed", seed,
                                     "_simResults.rds" ))






## Scenario 5: Null model with mid positive interaction  ------

# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1 <-  -0.2 # (~50%)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2 <-  -0.2 # (~50%)

# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12 <- 0.6 # decrease log-odds by 18%

beta = list(beta1, beta2, beta12)

seed = 1337
scen5 <-MSOM_simfit.fun.v2(beta = beta, 
                           nsites = nsites,
                           occ_formulas = occ_formulas,
                           det_formulas = det_formulas,
                           nspecies = nspecies, 
                           seed = seed, 
                           nsim = nsim, 
                           store.data = F,
                           J = J)



scenario <- "Scenario5_MidPos"
model.type = "Null_HighDet"

saveRDS(object = scen5,file = paste0("RevisionSimulations/Results/",scenario, model.type, "seed", seed,
                                     "_simResults.rds" ))


## Scenario 6: Null model with mid negative interaction  ------
# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1 <-  0.2 # (~50%)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2 <-  0.2 # (~50%)

# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12 <- -0.6 # decrease log-odds by 18%

beta = list(beta1, beta2, beta12)
seed = 1337
scen6 <- MSOM_simfit.fun.v2(beta = beta, 
                            nsites = nsites,
                            occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, 
                            seed = seed, 
                            nsim = nsim, 
                            store.data = F,
                            J = J)


scenario <- "Scenario6_MidNeg"
model.type = "Null_HighDet"

saveRDS(object = scen6,file = paste0("RevisionSimulations/Results/",scenario, model.type, "seed", seed,
                                     "_simResults.rds" ))

