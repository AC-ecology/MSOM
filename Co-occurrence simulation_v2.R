
############################################################
##    Simulation study to evaluate species co-occurrence  ##
##  model performance under different sample size regimes ##
############################################################

### Set up libary
library(dplyr, quietly = T)
library(unmarked, quietly = T)
library(ggplot2, quietly = T)

# Source simulation functions
source("MSOM_SimFun.R")


# General characteristics (constant across scenarios) ------------
# Sample size scenarios (nsites or N):  30 to 3,000 in steps of 10
ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)

# Number of simulations (nsim) per sample size in N: 100
nsim <- 100

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
               
### det formulas: formulas for the nspecieslinear predictors of detection.
###               Null models specified with "~1" in each linear predictor
det_formulas = rep("~1", nspecies)

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

seed = 1645
test <- MSOM_simfit.fun.v2(beta = beta, 
                            nsites = 100, 
                            occ_formulas = occ_formulas,
                            det_formulas = det_formulas,
                            nspecies = nspecies, 
                            seed = seed, 
                            nsim = 30, 
                            store.data = F,
                            J = J)

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

scenario <- "Scenario1"
model.type = "Null"

saveRDS(object = scen1,file = paste0("Results/",scenario, model.type, "seed", seed,
                                    "_simResults_v2.rds" ))


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
scenario <- "Scenario2"
model.type = "Null"

saveRDS(object = scen2,file = paste0("Results/",scenario, model.type, "seed", seed,
                                    "_simResults_v2.rds" ))

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


scenario <- "Scenario3"
model.type = "Null"

saveRDS(object = scen3,file = paste0("Results/",scenario, model.type, "seed", seed,
                                     "_simResults_v2.rds" ))

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



scenario <- "Scenario4"
model.type = "Null"

saveRDS(object = scen4,file = paste0("Results/",scenario, model.type, "seed", seed,
                                     "_simResults_v2.rds" ))






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
model.type = "Null"

saveRDS(object = scen5,file = paste0("Results/",scenario, model.type, "seed", seed,
                                     "_simResults_v2.rds" ))


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
model.type = "Null"

saveRDS(object = scen6,file = paste0("Results/",scenario, model.type, "seed", seed,
                                     "_simResults_v2.rds" ))


# Covariate models --------------

# Consistent through all scenarios:
nspecies <- 2
p_true <- rep(0.5, nspecies)
### det formulas: formulas for the nspecieslinear predictors of detection.
###               Null models specified with "~1" in each linear predictor
det_formulas <- rep("~1", nspecies)

## Scenario cov1: Covariate model with strong positive baseline interaction  ------

# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 


beta1 <- -0.4
# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1

beta2 <- -0.4
# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2
## beta12: baseline log-odds of co-occurrence
beta12.0 <- -0.6  
beta12.1 <- 0.8

beta12 <- c(beta12.0, beta12.1)#, beta12.2)

beta = list(beta1, beta2, beta12)

# General formulas to be:
# f1 = beta1.0 + beta1.1*occ_cov1 + beta1.2*occ_cov2
# f2 = beta2.0 + beta2.1*occ_cov2
# f12 = beta12.0 + beta12.1*occ_cov2 + beta12.2*occ_cov4


occ_formulas <- c("~1",          #f1
                  "~1",          #f2
                  "~occ_cov1")   #f12


## nocc_covs: number of occupancy (and co-occurrence covariates) to be generated
##            from standard normal distributions occ_cov~N(0,1) and to include in
##            the linear predictors for occupancy

nocc_covs <- 1

seed = 1337
scencov.1 <- MSOM_simfit.fun.v2(beta = beta, 
                            nsites = nsites,
                            occ_formulas = occ_formulas,
                            nspecies = nspecies, seed = seed, 
                            store.data = F,
                            nsim = nsim, J = J)

scenario <- "Scenario1_2sp"
model.type = "Covariate"

saveRDS(object = scencov.1,file = paste0("Results/",scenario, model.type, "seed", seed,
                                     "_simResults_v2.rds" ))

## Scenario cov2: Covariate model with weak positive baseline interaction  ------

# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1.0 <-  -0.2 # baseline prob(~40%)
beta1.1 <- 0.3 # occ_cov1 effect 

beta1 <- c(beta1.0, beta1.1)

# First order natural parameter f2: log-odds of presence of sp2 
## beta2: baseline occupancy of sp2, assumed equal to sp1
beta2.0 <-  -0.6 # (~35%)
beta2.1 <- 0.5

beta2 <- c(beta2.0, beta2.1)
# Second order natural parameter f12: log-odds of co-occurrene of sp1 & sp2

beta12.0 <- -0.6  
beta12.1 <- 0.8

beta12 <- c(beta12.0, beta12.1)


beta = list(beta1, beta2, beta12)

occ_formulas <- c("~occ_cov1",          #f1
                  "~occ_cov2",          #f2
                  "~occ_cov3")   #f12

nocc_covs = 3

seed = 1337
scencov.2 <- MSOM_simfit.fun.v2(beta = beta, 
                                nsites = nsites, 
                                occ_formulas = occ_formulas,
                                det_formulas = det_formulas,
                                nspecies = nspecies, 
                                store.data = F,
                                seed = seed, 
                                nsim = nsim, J = J)

scenario <- "Scenario2_2sp"
model.type = "NCovariate3"

saveRDS(object = scencov.2,file = paste0("Results/",scenario, model.type, "seed", seed,
                                         "_simResults_v2.rds" ))
## Scenario cov3: Covariate model with S = 2n and 2 covs in interaction  ------

# Parameter definition:
# First order natural parameter f1: log-odds of presence of sp1 
## beta1: baseline occupancy of sp1 
beta1.0 <-  -0.2 # baseline prob(~40%)
beta1.1 <- 0.4 # occ_cov1 effect 

beta1 <- c(beta1.0, beta1.1)

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

## Formulas and occupancy formulas


occ_formulas <- c("~occ_cov1",          #f1
                  "~occ_cov2",          #f2
                  "~occ_cov3 + occ_cov4")   #f12

nocc_covs = 4

seed = 1337
scencov.3 <- MSOM_simfit.fun.v2(beta = beta, 
                                nsites = nsites, 
                                occ_formulas = occ_formulas,
                                det_formulas = det_formulas, 
                                nocc_covs = nocc_covs,
                                nspecies = nspecies, 
                                seed = seed, 
                                store.data = F,
                                nsim = nsim, J = J)

scenario <- "Scenario3_2sp"
model.type = "NCovariate4_2int"

saveRDS(object = scencov.3,file = paste0("Results/",scenario, model.type, "seed", seed,
                                         "_simResults_v2.rds" ))
## Scenario cov4: Covariate model shared effect in interaction and f1  ------

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


occ_formulas <- c("~occ_cov1 + occ_cov4",          #f1
                  "~occ_cov2",          #f2
                  "~occ_cov3 + occ_cov4")   #f12

nocc_covs = 4

seed = 1337

scencov.4 <- MSOM_simfit.fun.v2(beta = beta, 
                                nsites = nsites, 
                                occ_formulas = occ_formulas,
                                det_formulas = det_formulas, 
                                nocc_covs = nocc_covs,
                                nspecies = nspecies, 
                                store.data = F,
                                seed = seed, 
                                nsim = nsim, J = J)

scenario <- "Scenario4_2sp"
model.type = "NCovariate4_SharedCovs"

saveRDS(object = scencov.4,file = paste0("Results/",scenario, model.type, "seed", seed,
                                         "_simResults_v2.rds" ))




# Network Models ---------
## 3 species simulation ------

ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)


# Number of simulations (nsim) per sample size in N: 100
nsim <- 100

# Number of visits (J): 3
J <- 3

nspecies <- 3
p_true <- rep(0.5, nspecies) # 50% detection probability

occ_formulas <-  rep("~1", nspecies + nspecies*(nspecies-1)/2)
det_formulas <- rep("~1", nspecies)

beta1 <- -0.4  #f1
beta2 <- 0.3   #f2
beta3 <- 0.8   #f3
beta12 <- 1    #f12
beta13 <- -0.7 #f13
beta23 <- -1   #f23

beta <-  list(beta1, beta2, beta3, beta12, beta13, beta23)

seed <- 1337
scen.3sp <- MSOM_simfit.fun.v2(beta = beta, p_true = p_true,
                               nsites = nsites, 
                               nspecies = nspecies,
                               store.data = F,
                               seed = seed, 
                               nsim = nsim, J = J)


scenario <- "Scenario1_3sp"
model.type = "Null"

# saveRDS(object = scen.3sp,file = paste0("Results/",scenario, model.type, "seed", seed,
#                                          "_simResults_v2.rds" ))

## 4 species simulation ----


ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)


# Number of simulations (nsim) per sample size in N: 100
nsim <- 100
# Number of visits (J): 3
J <- 3

nspecies <- 4
p_true <- rep(0.5, nspecies) # 50% detection probability

occ_formulas <-  rep("~1", nspecies + nspecies*(nspecies-1)/2)
det_formulas <- rep("~1", nspecies)

beta1 <- -0.4  #f1
beta2 <- 0.3   #f2
beta3 <- 0.8   #f3
beta4 <- -0.6  #f4
beta12 <- 1    #f12
beta13 <- -0.7 #f13
beta14 <- 0.4  #f14
beta23 <- -1   #f23
beta24 <- 0.6  #f24
beta34 <- -0.8 #f34
beta <-  list(beta1, beta2, beta3, beta4,
              beta12, beta13, beta14,
              beta23, beta24,
              beta34)

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

# saveRDS(object = scen.4sp,file = paste0("Results/",scenario, model.type, "seed", seed,
#                                         "_simResults_v2.rds" ))
## 5 species ----


ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
exp(ln.sites)


# Number of simulations (nsim) per sample size in N: 100
nsim <- 100

# Number of visits (J): 3
J <- 3

nspecies <- 5
p_true <- rep(0.5, nspecies) # 50% detection probability

occ_formulas <-  rep("~1", nspecies + nspecies*(nspecies-1)/2)
det_formulas <- rep("~1", nspecies)

beta1 <- -0.4  #f1
beta2 <- 0.3   #f2
beta3 <- 0.8   #f3
beta4 <- -0.6  #f4
beta5 <- 0.5   #
beta12 <- 1    #f12
beta13 <- -0.7 #f13
beta14 <- 0.4  #f14
beta15 <- 0.5  #
beta23 <- -1   #f23
beta24 <- 0.6  #f24
beta25 <- -0.4 #
beta34 <- -0.8 #f34
beta35 <- 0.65 
beta45 <- 0.45
beta <-  list(beta1, beta2, beta3, beta4, beta5,
              beta12, beta13, beta14, beta15,
              beta23, beta24, beta25, 
              beta34, beta35,
              beta45)


seed <- 1337
scen.5sp <- MSOM_simfit.fun.v2(beta = beta, 
                               p_true = p_true,
                               nsites = nsites, 
                               nspecies = nspecies,
                               seed = seed, 
                               store.data = F,
                               nsim = nsim, J = J)


scen.5sp$time.ellapsed

scenario <- "Scenario1_5sp"
model.type = "Null"

# saveRDS(object = scen.5sp,file = paste0("Results/",scenario, model.type, "seed", seed,
#                                         "_simResults_v2.rds" ))



