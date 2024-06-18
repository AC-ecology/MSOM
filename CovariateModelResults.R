##############################################
##### Covariate model results processing #####
##############################################

# Workflow:
## 1. Load result object
## 2. Calculate bias, coverage rate, CV and power (for 2nd order and covariate terms)
## 3. Calculate bias, coverage rate for general parameters
## 4. Summarise results per sample size scenario
## 5. Combine Normal and Penalized likelihood results for each scenario
## 6. Plot with ggplot, faceting for parameter



library(readr)
library(tidyverse)

### Shared features:
ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))

nsim <- 100

nspecies <- 2

# function for retreiving general params
gen.param.fun <- function(beta, nspecies = 2, 
                          occ_formulas = NULL, occ_covs = NULL, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  if(is.null(occ_formulas)){occ_formulas <- rep("~1", nspecies)} else{
    nocc_covs <- max(as.numeric(gsub(pattern = "[a-z]", "",
                                     str_extract(occ_formulas, "cov[0-9]"))), na.rm = T)
    
  }
 
  if(is.null(occ_covs)) {
    occ_covs <- data.frame(matrix(rep(seq(-1, 1, 
                                          length.out  = 100), nocc_covs),
                                  ncol = nocc_covs))
    names(occ_covs) <- paste('occ_cov',1:nocc_covs,sep='')
  } 
  names(occ_covs) <- paste('occ_cov',1:ncol(occ_covs),sep='')
  n.nat.params <- nspecies + nspecies*(nspecies-1)/2
  nat.params <- paste0("f", c(1:nspecies, 
                              apply(combn(1:nspecies, 2),
                                    2,paste0, collapse = "" )))
  
  
  occ_dms <- lapply(X = as.list(occ_formulas), 
                    function(x)model.matrix(as.formula(x), data = occ_covs))
  
  for(p in 1:n.nat.params){
    if(ncol(occ_dms[[p]]) != length(beta[[p]])){
      stop(paste0("Error: Number of regression coefficients for linear predictor ",
                  nat.params[p]," must equal number of covariates in formula + 1"))
    }
    
    if(p == 1){
      f <- occ_dms[[p]]%*%beta[[p]]
      params <- paste0(nat.params[p], ".", seq(0,length(beta[[p]])-1))
    } else{
      
      f <- cbind(f,occ_dms[[p]]%*%beta[[p]] )
      params <- c(params, paste0(nat.params[p], ".", seq(0,length(beta[[p]])-1)))
      if(p > (nspecies) && all(beta[[p]] == 0)){
        occ_formulas[p] <- '0'
      }
    }
    if(p == n.nat.params & nspecies > 2){
      f <-  cbind(f, matrix(rep(0, N*(2^nspecies-1 - n.nat.params)), ncol = (2^nspecies-1 - n.nat.params)))
    }
  }
  if(any(occ_formulas == "0")){
    params <- params[-which(occ_formulas == "0")]
    og.param.val <- unlist(beta[-which(occ_formulas == "0")])
  } else{
    og.param.val <- unlist(beta)
  }
  
  
  #f <- cbind(f1,f2,f12)
  z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
  colnames(z) <- paste('sp',1:nspecies,sep='')
  dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)
  
  psi <- exp(f %*% t(dm))
  psi <- psi/rowSums(psi)
  
  colnames(psi) <- apply(z, 1, function(x) paste0(x, collapse = ""))
  return(psi)
}

# Scenario 1 ----

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

## Read object
scen1 <- scen.cov1 <- read_rds("Results/Scenario1_2spCovariateseed1337_simResults_v2.rds")


## Normal likelihood 
State_simcov1 <- scen.cov1$State.params %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)), 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam))
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

scen.1.conv <-State_simcov1 %>%
  group_by(n.sites, Parameter) %>%
  summarise(Conv.rate = mean(conv)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_simcov1.sum <- State_simcov1 %>%
  group_by(n.sites, Parameter) %>%
  filter(conv == 1) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias), mu.CV = mean(CV, na.rm = T),
            og.val = unique(og_param.val),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered, na.rm = T),                    # Coverage rate 
            PWR = mean(below.alpha, na.rm = T),               # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate (should be 100%)
            ndatasets = n(),                                  # Number datasets in summary
            Lik = "LL") %>% 
  ungroup()





#### General parameters (psi11, psi10, psi01, psi00)
occ_covs = data.frame(occ_cov1 = 1)



GenParam_sim1 <- State_simcov1 %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = list(Estimate[1],
                                                            Estimate[2],
                                                            Estimate[3:4]), 
                                                occ_formulas = occ_formulas,
                                                occ_covs = occ_covs))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               occ_formulas = occ_formulas,
                                               occ_covs = occ_covs))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             occ_formulas = occ_formulas,
                                             occ_covs = occ_covs))) %>%
  ungroup


GenParam_sim1 <-  GenParam_sim1 %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim1.sum <- GenParam_sim1 %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()


#### Derived Parameters 

Marg.prob1 <-  GenParam_sim1 %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                         sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()



  



Cond.prob1 <-  GenParam_sim1 %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1", "Sp2=1|Sp1=0", "Sp2=1|Sp1=1"),
            # Estimated
            Cond.est = c(psi.est[Gen.Par %in% c("10")]/sum(psi.est[Gen.Par %in% c("10", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "01")]),
                         psi.est[Gen.Par %in% c("01")]/sum(psi.est[Gen.Par %in% c("01", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "10")])),
            # Original
            Cond.og =  c(og.psi[Gen.Par %in% c("10")]/sum(og.psi[Gen.Par %in% c("10", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "01")]),
                         og.psi[Gen.Par %in% c("01")]/sum(og.psi[Gen.Par %in% c("01", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "10")])),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()


## Summaries of derived

Marg.prob1.sum <- Marg.prob1 %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,  ndatasets = n(),               # Upper bias CI
            Lik = "LL") %>% 
  ungroup()



Cond.prob1.sum <- Cond.prob1 %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,    ndatasets = n(),             # Upper bias CI
            Lik = "LL") %>% 
  ungroup()



## Penalised likelihood 

State_simcov1.pl <- scen1$State.params.pl %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)), 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_simcov1.pl.sum <- State_simcov1.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),mu.CV = mean(CV, na.rm = T),
            og.val = unique(og_param.val),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),  ndatasets = n(),                    # Mean convergence rate
            Lik = "PL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim1.pl <- State_simcov1.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = list(Estimate[1],
                                                            Estimate[2],
                                                            Estimate[3:4]), 
                                                occ_formulas = occ_formulas,
                                                occ_covs = occ_covs))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               occ_formulas = occ_formulas,
                                               occ_covs = occ_covs))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             occ_formulas = occ_formulas,
                                             occ_covs = occ_covs))) %>%
  ungroup()


GenParam_sim1.pl <-  GenParam_sim1.pl %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim1.pl.sum <- GenParam_sim1.pl %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "PL",  ndatasets = n()) %>% 
  ungroup()
  

#### Derived Parameters 

Marg.prob1.pl <-  GenParam_sim1.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()




Cond.prob1.pl <-  GenParam_sim1.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1", "Sp2=1|Sp1=0", "Sp2=1|Sp1=1"),
            # Estimated
            Cond.est = c(psi.est[Gen.Par %in% c("10")]/sum(psi.est[Gen.Par %in% c("10", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "01")]),
                         psi.est[Gen.Par %in% c("01")]/sum(psi.est[Gen.Par %in% c("01", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "10")])),
            # Original
            Cond.og =  c(og.psi[Gen.Par %in% c("10")]/sum(og.psi[Gen.Par %in% c("10", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "01")]),
                         og.psi[Gen.Par %in% c("01")]/sum(og.psi[Gen.Par %in% c("01", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "10")])),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()


## Summaries of derived

Marg.prob1.pl.sum <- Marg.prob1.pl %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "PL",  ndatasets = n()) %>% 
  ungroup()



Cond.prob1.pl.sum <- Cond.prob1.pl %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "PL",  ndatasets = n()) %>% 
  ungroup()


###### Combining LL and PL into one for plotting

## Natural parameters
full.scen1.nat <- rbind(State_simcov1.sum, State_simcov1.pl.sum)

## General parameters

full.scen1.gen <- rbind(GenParam_sim1.sum, GenParam_sim1.pl.sum)

## Marginal probabilities
full.scen1.mar <- rbind(Marg.prob1.sum, Marg.prob1.pl.sum)

## Conditional probabilities
full.scen1.con <- rbind(Cond.prob1.sum, Cond.prob1.pl.sum)

# Scenario 2 ----

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


occ_covs = data.frame(occ_cov1 = 1,
                      occ_cov2 = 1,
                      occ_cov3 = 1)
## Read object
scen2 <- read_rds("Results/Scenario2_2spNCovariate3seed1337_simResults_v2.rds")




## Normal likelihood 
State_sim2 <- scen2$State.params %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)), 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim2.sum <- State_sim2 %>%
  filter(conv == 1) %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            mu.CV = mean(CV, na.rm = T),
            og.val = unique(og_param.val),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim2 <- State_sim2 %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = list(Estimate[1:2],
                                                            Estimate[3:4],
                                                            Estimate[5:6]), 
                                                occ_formulas = occ_formulas,
                                                occ_covs = occ_covs))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               occ_formulas = occ_formulas,
                                               occ_covs = occ_covs))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             occ_formulas = occ_formulas,
                                             occ_covs = occ_covs))) %>%
  ungroup()


GenParam_sim2 <-  GenParam_sim2 %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim2.sum <- GenParam_sim2 %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL", ndatasets = n()) %>% 
  ungroup()


#### Derived Parameters 

Marg.prob2 <-  GenParam_sim2 %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()




Cond.prob2 <-  GenParam_sim2 %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1", "Sp2=1|Sp1=0", "Sp2=1|Sp1=1"),
            # Estimated
            Cond.est = c(psi.est[Gen.Par %in% c("10")]/sum(psi.est[Gen.Par %in% c("10", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "01")]),
                         psi.est[Gen.Par %in% c("01")]/sum(psi.est[Gen.Par %in% c("01", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "10")])),
            # Original
            Cond.og =  c(og.psi[Gen.Par %in% c("10")]/sum(og.psi[Gen.Par %in% c("10", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "01")]),
                         og.psi[Gen.Par %in% c("01")]/sum(og.psi[Gen.Par %in% c("01", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "10")])),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()



## Summaries of derived

Marg.prob2.sum <- Marg.prob2 %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()



Cond.prob2.sum <- Cond.prob2 %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()


## Penalised likelihood 

State_sim2.pl <- scen2$State.params.pl %>%
  mutate(conv = ifelse(any(is.nan(SE)), 0, conv)) %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim2.pl.sum <- State_sim2.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),mu.CV = mean(CV, na.rm = T),
            og.val = unique(og_param.val),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim2.pl <- State_sim2.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = list(Estimate[1:2],
                                                            Estimate[3:4],
                                                            Estimate[5:6]), 
                                                occ_formulas = occ_formulas,
                                                occ_covs = occ_covs))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               occ_formulas = occ_formulas,
                                               occ_covs = occ_covs))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             occ_formulas = occ_formulas,
                                             occ_covs = occ_covs))) %>%
  ungroup()


GenParam_sim2.pl <-  GenParam_sim2.pl %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim2.pl.sum <- GenParam_sim2.pl %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "PL", ndatasets = n()) %>% 
  ungroup()


#### Derived Parameters 

Marg.prob2.pl <-  GenParam_sim2.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()




Cond.prob2.pl <-  GenParam_sim2.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1", "Sp2=1|Sp1=0", "Sp2=1|Sp1=1"),
            # Estimated
            Cond.est = c(psi.est[Gen.Par %in% c("10")]/sum(psi.est[Gen.Par %in% c("10", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "01")]),
                         psi.est[Gen.Par %in% c("01")]/sum(psi.est[Gen.Par %in% c("01", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "10")])),
            # Original
            Cond.og =  c(og.psi[Gen.Par %in% c("10")]/sum(og.psi[Gen.Par %in% c("10", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "01")]),
                         og.psi[Gen.Par %in% c("01")]/sum(og.psi[Gen.Par %in% c("01", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "10")])),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()

## Summaries of derived

Marg.prob2.pl.sum <- Marg.prob2.pl %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "PL") %>% 
  ungroup()



Cond.prob2.pl.sum <- Cond.prob2.pl %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "PL") %>% 
  ungroup()




###### Combining LL and PL into one for plotting

## Natural parameters
full.scen2.nat <- rbind(State_sim2.sum, State_sim2.pl.sum)
## General parameters

full.scen2.gen <- rbind(GenParam_sim2.sum, GenParam_sim2.pl.sum)

## Marginal probabilities
full.scen2.mar <- rbind(Marg.prob2.sum, Marg.prob2.pl.sum)

## Conditional probabilities
full.scen2.con <- rbind(Cond.prob2.sum, Cond.prob2.pl.sum)


# Scenario 3 ----

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




occ_covs = data.frame(occ_cov1 = 1,
                      occ_cov2 = 1,
                      occ_cov3 = 1,
                      occ_cov4 = 1)
## Read object
scen3 <- read_rds("Results/Scenario3_2spNCovariate4_2intseed1337_simResults_v2.rds")


ok <- scen3$State.params$SE

## Normal likelihood 
State_sim3 <- scen3$State.params %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)) |any(is.na(SE)) , 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim3.sum <- State_sim3 %>%
  filter(conv == 1) %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),mu.CV = mean(CV, na.rm = T),
            og.val = unique(og_param.val),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            Lik = "LL",
            ndatasets = n()) %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim3 <- State_sim3 %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = list(Estimate[1:2],
                                                            Estimate[3:4],
                                                            Estimate[5:7]), 
                                                occ_formulas = occ_formulas,
                                                occ_covs = occ_covs))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               occ_formulas = occ_formulas,
                                               occ_covs = occ_covs))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             occ_formulas = occ_formulas,
                                             occ_covs = occ_covs))) %>%
  ungroup()


GenParam_sim3 <-  GenParam_sim3 %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim3.sum <- GenParam_sim3 %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL", ndatasets = n()) %>% 
  ungroup()


#### Derived Parameters 

Marg.prob3 <-  GenParam_sim3 %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()




Cond.prob3 <-  GenParam_sim3 %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1", "Sp2=1|Sp1=0", "Sp2=1|Sp1=1"),
            # Estimated
            Cond.est = c(psi.est[Gen.Par %in% c("10")]/sum(psi.est[Gen.Par %in% c("10", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "01")]),
                         psi.est[Gen.Par %in% c("01")]/sum(psi.est[Gen.Par %in% c("01", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "10")])),
            # Original
            Cond.og =  c(og.psi[Gen.Par %in% c("10")]/sum(og.psi[Gen.Par %in% c("10", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "01")]),
                         og.psi[Gen.Par %in% c("01")]/sum(og.psi[Gen.Par %in% c("01", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "10")])),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()


## Summaries of derived

Marg.prob3.sum <- Marg.prob3 %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()



Cond.prob3.sum <- Cond.prob3 %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()


## Penalised likelihood 

State_sim3.pl <- scen3$State.params.pl %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)) |any(is.na(SE)) , 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim3.pl.sum <- State_sim3.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),mu.CV = mean(CV, na.rm = T),
            og.val = unique(og_param.val),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim3.pl <- State_sim3.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = list(Estimate[1:2],
                                                            Estimate[3:4],
                                                            Estimate[5:7]), 
                                                occ_formulas = occ_formulas,
                                                occ_covs = occ_covs))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               occ_formulas = occ_formulas,
                                               occ_covs = occ_covs))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             occ_formulas = occ_formulas,
                                             occ_covs = occ_covs))) %>%
  ungroup()


GenParam_sim3.pl <-  GenParam_sim3.pl %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim3.pl.sum <- GenParam_sim3.pl %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()


#### Derived Parameters 

Marg.prob3.pl <-  GenParam_sim3.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()




Cond.prob3.pl <-  GenParam_sim3.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1", "Sp2=1|Sp1=0", "Sp2=1|Sp1=1"),
            # Estimated
            Cond.est = c(psi.est[Gen.Par %in% c("10")]/sum(psi.est[Gen.Par %in% c("10", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "01")]),
                         psi.est[Gen.Par %in% c("01")]/sum(psi.est[Gen.Par %in% c("01", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "10")])),
            # Original
            Cond.og =  c(og.psi[Gen.Par %in% c("10")]/sum(og.psi[Gen.Par %in% c("10", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "01")]),
                         og.psi[Gen.Par %in% c("01")]/sum(og.psi[Gen.Par %in% c("01", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "10")])),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()

## Summaries of derived

Marg.prob3.pl.sum <- Marg.prob3.pl %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()



Cond.prob3.pl.sum <- Cond.prob3.pl %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()


###### Combining LL and PL into one for plotting

## Natural parameters
full.scen3.nat <- rbind(State_sim3.sum, State_sim3.pl.sum)

## General parameters

full.scen3.gen <- rbind(GenParam_sim3.sum, GenParam_sim3.pl.sum)

## Marginal probabilities
full.scen3.mar <- rbind(Marg.prob3.sum, Marg.prob3.pl.sum)

## Conditional probabilities
full.scen3.con <- rbind(Cond.prob3.sum, Cond.prob3.pl.sum)



# Scenario 4 ----


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


occ_covs = data.frame(occ_cov1 = 1,
                      occ_cov2 = 1,
                      occ_cov3 = 1,
                      occ_cov4 = 1)
## Read object
scen4 <- read_rds("Results/Scenario4_2spNCovariate4_SharedCovsseed1337_simResults_v2.rds")




## Normal likelihood 
State_sim4 <- scen4$State.params %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)) |any(is.na(SE)) , 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim4.sum <- State_sim4 %>%
  filter(conv == 1) %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),mu.CV = mean(CV, na.rm = T),
            og.val = unique(og_param.val),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim4 <- State_sim4 %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = list(Estimate[1:3],
                                                            Estimate[4:5],
                                                            Estimate[6:8]), 
                                                occ_formulas = occ_formulas,
                                                occ_covs = occ_covs))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               occ_formulas = occ_formulas,
                                               occ_covs = occ_covs))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             occ_formulas = occ_formulas,
                                             occ_covs = occ_covs))) %>%
  ungroup()


GenParam_sim4 <-  GenParam_sim4 %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim4.sum <- GenParam_sim4 %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()


#### Derived Parameters 

Marg.prob4 <-  GenParam_sim4 %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()




Cond.prob4 <-  GenParam_sim4 %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1", "Sp2=1|Sp1=0", "Sp2=1|Sp1=1"),
            # Estimated
            Cond.est = c(psi.est[Gen.Par %in% c("10")]/sum(psi.est[Gen.Par %in% c("10", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "01")]),
                         psi.est[Gen.Par %in% c("01")]/sum(psi.est[Gen.Par %in% c("01", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "10")])),
            # Original
            Cond.og =  c(og.psi[Gen.Par %in% c("10")]/sum(og.psi[Gen.Par %in% c("10", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "01")]),
                         og.psi[Gen.Par %in% c("01")]/sum(og.psi[Gen.Par %in% c("01", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "10")])),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()

## Summaries of derived

Marg.prob4.sum <- Marg.prob4 %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()



Cond.prob4.sum <- Cond.prob4 %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()

## Penalised likelihood 

State_sim4.pl <- scen4$State.params.pl %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)) |any(is.na(SE)) , 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim4.pl.sum <- State_sim4.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias), mu.CV = mean(CV, na.rm = T),
            og.val = unique(og_param.val),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim4.pl <- State_sim4.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = list(Estimate[1:3],
                                                            Estimate[4:5],
                                                            Estimate[6:8]), 
                                                occ_formulas = occ_formulas,
                                                occ_covs = occ_covs))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               occ_formulas = occ_formulas,
                                               occ_covs = occ_covs))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             occ_formulas = occ_formulas,
                                             occ_covs = occ_covs))) %>%
  ungroup()


GenParam_sim4.pl <-  GenParam_sim4.pl %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim4.pl.sum <- GenParam_sim4.pl %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()


#### Derived Parameters 

Marg.prob4.pl <-  GenParam_sim4.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()




Cond.prob4.pl <-  GenParam_sim4.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1", "Sp2=1|Sp1=0", "Sp2=1|Sp1=1"),
            # Estimated
            Cond.est = c(psi.est[Gen.Par %in% c("10")]/sum(psi.est[Gen.Par %in% c("10", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "01")]),
                         psi.est[Gen.Par %in% c("01")]/sum(psi.est[Gen.Par %in% c("01", "00")]),
                         psi.est[Gen.Par %in% c("11")]/sum(psi.est[Gen.Par %in% c("11", "10")])),
            # Original
            Cond.og =  c(og.psi[Gen.Par %in% c("10")]/sum(og.psi[Gen.Par %in% c("10", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "01")]),
                         og.psi[Gen.Par %in% c("01")]/sum(og.psi[Gen.Par %in% c("01", "00")]),
                         og.psi[Gen.Par %in% c("11")]/sum(og.psi[Gen.Par %in% c("11", "10")])),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()

## Summaries of derived

Marg.prob4.pl.sum <- Marg.prob4.pl %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()



Cond.prob4.pl.sum <- Cond.prob4.pl %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()


###### Combining LL and PL into one for plotting

## Natural parameters
full.scen4.nat <- rbind(State_sim4.sum, State_sim4.pl.sum)

## General parameters

full.scen4.gen <- rbind(GenParam_sim4.sum, GenParam_sim4.pl.sum)

## Marginal probabilities
full.scen4.mar <- rbind(Marg.prob4.sum, Marg.prob4.pl.sum)

## Conditional probabilities
full.scen4.con <- rbind(Cond.prob4.sum, Cond.prob4.pl.sum)



## The plot parade ----
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
### Scenario 1 -----
#### Natural parameters

##### Prop.bias
g1.nat <- ggplot(full.scen1.nat, aes(x = log(n.sites), y = mu.p.bias, 
                                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 1 (Str +ve)", col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")+
  coord_cartesian(ylim =  c(-2.5, 2.5)) +guides(alpha = "none")



##### Power

g1.pwr <- ggplot(full.scen1.nat%>%
                   filter(grepl("f12", Parameter) | grepl(".1", Parameter)),
                 aes(x = log(n.sites), y = PWR,
                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5, pch = 18)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Power", title = "Scenario Cov1",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(0, 1)+
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "red")

#### General parameters
g1.gen <- ggplot(full.scen1.gen, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5, pch = 15,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Gen.Par)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 1 (Str +ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


#### Derived parameters

##### Marginal probs
g1.mar <-  ggplot(full.scen1.mar, aes(x = log(n.sites), y = mu.p.bias,
                                      group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Species)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias",
       title = "Scenario 1 (Str +ve)", col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")

##### conditional probs
g1.con <- ggplot(full.scen1.con, aes(x = log(n.sites), y = mu.p.bias,
                                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5, pch = 15,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Cond.prob)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 1 (Str +ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")



### Scenario 2 ------
#### Natural parameters

##### Prop.bias
g2.nat <- ggplot(full.scen2.nat, aes(x = log(n.sites), y = mu.p.bias,
                                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 2 (1 Cov per predictor)", col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")+
  coord_cartesian(ylim =  c(-2.5, 2.5))


##### Power

g2.pwr <- ggplot(full.scen2.nat%>%
                   filter(grepl("f12", Parameter) |grepl("[0-9].1", Parameter) ),
                 aes(x = log(n.sites), y = PWR,
                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5, pch = 18)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Power", title = "Scenario Cov1",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(0, 1)+
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "red")

#### General parameters
g2.gen <- ggplot(full.scen2.gen, aes(x = log(n.sites), y = mu.p.bias,
                                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5, pch = 15,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Gen.Par)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 2 (Weak +ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


#### Derived parameters

##### Marginal probs
g2.mar <-  ggplot(full.scen2.mar, aes(x = log(n.sites), y = mu.p.bias,
                                      group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Species)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 2 (Weak +ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")

##### conditional probs
g2.con <- ggplot(full.scen2.con, aes(x = log(n.sites), y = mu.p.bias,
                                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5, pch = 15,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Cond.prob)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 2 (Weak +ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


### Scenario 3 -----
#### Natural parameters

##### Prop.bias
g3.nat <- ggplot(full.scen3.nat, aes(x = log(n.sites), y = mu.p.bias,
                                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 3 (Str -ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-5.5, 5.5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


##### Power

g3.pwr <- ggplot(full.scen3.nat%>%
                   filter(Parameter == "beta12.0"), aes(x = log(n.sites), y = PWR,
                                                        group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5, pch = 18)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  #facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Power", title = "Scenario 3 (Str -ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(0, 1)+
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "red")

#### General parameters
g3.gen <- ggplot(full.scen3.gen, aes(x = log(n.sites), y = mu.p.bias,
                                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5, pch = 15,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Gen.Par)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 3 (Str -ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-.8, .8)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


#### Derived parameters

##### Marginal probs
g3.mar <-  ggplot(full.scen3.mar, aes(x = log(n.sites), y = mu.p.bias,
                                      group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Species)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 3 (Str -ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")

##### conditional probs
g3.con <- ggplot(full.scen3.con, aes(x = log(n.sites), y = mu.p.bias,
                                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5, pch = 15,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Cond.prob)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 3 (Str -ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")



### Scenario 4 -----
#### Natural parameters

##### Prop.bias
g4.nat <- ggplot(full.scen4.nat, aes(x = log(n.sites), y = mu.p.bias,
                                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 4 (Weak -ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-3.5, 3.5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")




##### Power

g4.pwr <- ggplot(full.scen4.nat%>%
                   filter(Parameter == "beta12.0"), aes(x = log(n.sites), y = PWR,
                                                        group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  #facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Power", title = "Scenario 4 (Weak -ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(0, 1)+
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "red")

#### General parameters
g4.gen <- ggplot(full.scen4.gen, aes(x = log(n.sites), y = mu.p.bias,group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5, pch = 15,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Gen.Par)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 4 (Weak -ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-.8, .8)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")



#### Derived parameters

##### Marginal probs
g4.mar <-  ggplot(full.scen4.mar, aes(x = log(n.sites), y = mu.p.bias,
                                      group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Species)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 4 (Weak -ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")

##### conditional probs
g4.con <- ggplot(full.scen4.con, aes(x = log(n.sites), y = mu.p.bias,
                                     group = Lik, col = Lik))+
  geom_line()+
  geom_point(size = 2.5, pch = 15,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Cond.prob)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "Scenario 4 (Weak -ve)",
       col = "Likelihood")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")



## All scenarios together -----


full.covscen.nat <- rbind(full.scen1.nat %>% mutate(Scenario = "ScenCov1"),
                          full.scen2.nat %>% mutate(Scenario = "ScenCov2"),
                          full.scen3.nat %>% mutate(Scenario = "ScenCov3"),
                          full.scen4.nat %>% mutate(Scenario = "ScenCov4"))
full.covscen.nat$order <-  ifelse(nchar(full.covscen.nat$Parameter)>4, "2nd", "1st")

full.covscen.nat$CoefReg <- ifelse(grepl("\\.0",full.covscen.nat$Parameter, perl = F), "Intercept", "Slope")

### Power plot(s)


# Colour by Coeficent Regression type (Int v Slope) facet by Scenario
ggplot(full.covscen.nat %>% filter(order == "2nd" | CoefReg == "Slope"), 
       aes(x = log(n.sites), y = PWR, group = interaction(Parameter, Lik)))+
  geom_hline(yintercept = 0.95, col = "red", linetype = "dotted")+
  geom_line(aes( linetype = Lik, col = CoefReg))+
  facet_grid(Scenario~Lik)+
  labs(x = "Log(Number of Sites)", y = "Power (Type I error)")+
  scale_y_continuous(labels = scales::percent_format())+
  scale_colour_manual(values = cbbPalette[c(4, 8)])


# Colour by effect size (Int v Slope) facet by Likelihood
ggplot(full.covscen.nat %>% filter(order == "2nd" | CoefReg == "Slope"), 
       aes(x = log(n.sites), y = PWR, group = interaction(Parameter, Lik)))+
  geom_hline(yintercept = 0.95, col = "red", linetype = "dotted")+
  geom_line(aes( linetype = Lik, col = abs(og.val)))+
  facet_grid(Scenario~Lik)+
  labs(x = "Log(Number of Sites)", y = "Power (Type I error)", col = "Effect Size")+
  scale_y_continuous(labels = scales::percent_format())+
  scale_colour_viridis_c()

# Colour by natural parameter order type (Int v Slope) facet by Scenario
ggplot(full.covscen.nat %>% filter(order == "2nd" | CoefReg == "Slope"), 
       aes(x = log(n.sites), y = PWR, group = interaction(Parameter, Lik)))+
  geom_hline(yintercept = 0.95, col = "red", linetype = "dotted")+
  geom_line(aes( linetype = Lik, col = order))+
  facet_grid(Scenario~Lik)+
  labs(x = "Log(Number of Sites)", y = "Power (Type I error)")+
  scale_y_continuous(labels = scales::percent_format())+
  scale_colour_manual(values = cbbPalette[c(4, 8)])



## Bias plots
### Natural parameters 
ggplot(full.covscen.nat, 
       aes(x = log(n.sites), y = mu.p.bias, group = interaction(Parameter, Lik)))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "black")+
  geom_line(aes(col = interaction(order, Lik)))+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased"),
                 col = interaction(order, Lik)))+
  scale_alpha_manual(values = c(1, 0.3))+
  facet_grid(Scenario~Lik)+
  labs(x = "Log(Number of Sites)", y = "Relative Bias", col = "Order-Likelihood")+ guides(alpha = "none")+
  coord_cartesian(ylim =  c(-1.5, 1.5))+
  scale_y_continuous(labels = scales::percent_format() )+
  scale_colour_manual(values = cbbPalette[c(4:9)])








### Cowplot parade! ----

## Bias plots ----

### Natural parameters 
#### Positive scenarios
cowplot::plot_grid(g1.nat+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank())+
                     guides(alpha = "none"),
                   g2.nat+guides(col = "none"), nrow = 2)

#### Negative scenarios
cowplot::plot_grid(g3.nat+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank()),
                   g4.nat, nrow = 2)


### General parameters 
#### Positive scenarios
cowplot::plot_grid(g1.gen+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank()),
                   g2.gen, nrow = 2)

#### Negative scenarios
cowplot::plot_grid(g3.gen+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank()),
                   g4.gen, nrow = 2)

## Power plots ----

cowplot::plot_grid(g1.pwr+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                legend.position = "none"),
                   g2.pwr+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank(),
                                axis.text.y = element_blank(),
                                axis.title.y = element_blank()),
                   g3.pwr+theme(legend.position = "none"), g4.pwr+theme(axis.text.y = element_blank(),
                                                                        axis.title.y = element_blank()), nrow = 2)


###  Marginal probabilities
#### Positive scenarios
cowplot::plot_grid(g1.mar+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank()),
                   g2.mar, nrow = 2)

#### Negative scenarios
cowplot::plot_grid(g3.mar+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank()),
                   g4.mar, nrow = 2)



###  Conditional probabilities
#### Positive scenarios
cowplot::plot_grid(g1.con+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank()),
                   g2.con, nrow = 2)

#### Negative scenarios
cowplot::plot_grid(g3.con+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank()),
                   g4.con, nrow = 2)





