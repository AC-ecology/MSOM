##############################################
##### Covariate model results processing #####
##############################################

# Workflow:
## 1. Load result object
## 2. Calculate bias, coverage rate, CV and power (for 2nd order and covariate terms)
## 3. Calculate bias, coverage rate for general parameters
## 4. Summarise results per sample size scenario
## 5. Combine Normal and Penalized likelihood results for each scenario
## 6. Plot with ggplot and make relative bias tables

### WARNING: Tables are produced to insert in a LaTeX editor, if you want
###          to visualise in RStudio's viewer change `format = "latex"` to
###          `format = "html` but bear in mind it might not render properly

library(readr)
library(tidyverse)

### Shared features:
ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))

nsim <- 100

nspecies <- 2

# function for retreiving general params

## INPUT:

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
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()



Cond.prob2.sum <- Cond.prob2 %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
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
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()



Cond.prob2.pl.sum <- Cond.prob2.pl %>%
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
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette)

## All scenarios together -----


full.covscen.nat <- rbind(full.scen1.nat %>% mutate(Scenario = "1 Slope"),
                          full.scen2.nat %>% mutate(Scenario = "3 Slopes"),
                          full.scen3.nat %>% mutate(Scenario = "4 Slopes"),
                          full.scen4.nat %>% mutate(Scenario = "5 Slopes"))
full.covscen.nat$order <-  ifelse(nchar(full.covscen.nat$Parameter)>4, "2nd", "1st")

full.covscen.nat$CoefReg <- ifelse(grepl("\\.0",full.covscen.nat$Parameter, perl = F), "Intercept", "Slope")


### General probs

full.covscen.gen <- rbind(full.scen1.gen %>% mutate(Scenario = "1 Slope"),
                          full.scen2.gen %>% mutate(Scenario = "3 Slopes"),
                          full.scen3.gen %>% mutate(Scenario = "4 Slopes"),
                          full.scen4.gen %>% mutate(Scenario = "5 Slopes"))

full.covscen.mar <- rbind(full.scen1.mar %>% mutate(Scenario = "1 Slope"),
                          full.scen2.mar %>% mutate(Scenario = "3 Slopes"),
                          full.scen3.mar %>% mutate(Scenario = "4 Slopes"),
                          full.scen4.mar %>% mutate(Scenario = "5 Slopes"))


### Power plot(s) ------
lik_labels <- c("LL"= "Log-Likelihood", "PL"="Penalised Likelihood")

# Colour by Coeficent Regression type (Int v Slope) facet by Scenario
ggplot(full.covscen.nat %>% filter(order == "2nd" | CoefReg == "Slope")%>% filter(n.sites >33), 
       aes(x = log(n.sites), y = PWR, group = interaction(Parameter, Lik)))+
  geom_hline(yintercept = 0.95, col = "red", linetype = "dotted")+
  geom_line(aes( linetype = Lik, col = CoefReg), lwd = 0.8)+
  facet_grid(Scenario~Lik)+
  labs(x = "Log(Number of Sites)", y = "Power (Type I error)")+
  scale_y_continuous(labels = scales::percent_format())+
  scale_colour_manual(values = cbbPalette[c(5, 6)])+
  theme_bw()


# Colour by effect size (Int v Slope) facet by Likelihood
ggplot(full.covscen.nat %>% filter(order == "2nd" | CoefReg == "Slope")%>%
         filter(n.sites >33), 
       aes(x = log(n.sites), y = PWR, group = interaction(Parameter, Lik)))+
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "black", lwd = 0.15)+
  geom_line(aes(col = abs(og.val)), lwd = 0.3)+
  # geom_point(size = 0.4,
  #            aes(col = abs(og.val)))+
  facet_grid(Scenario~Lik, labeller = labeller(Lik = lik_labels))+
  labs(x = "Number of Sites", y = "Power", col = "Effect Size")+ guides(alpha = "none")+
  coord_cartesian(ylim =  c(0, 1))+
  scale_x_continuous(breaks = log(nsites), labels = nsites)+
  scale_y_continuous(labels = scales::percent_format() )+
  scale_colour_gradientn(colours =safe_colorblind_palette[2:5] )+
  # scale_colour_viridis_c()+
  theme_bw()+
  theme(strip.text.y = element_text(size = 3.5, color = "black",  margin = margin(r = 2, l = 2)),
        strip.text.x = element_text(size = 3.5, color = "black", margin = margin(b = 2, t = 2)),
        axis.text.y = element_text(size = 3.5, color = "black"),
        axis.text.x =element_text(size = 3, color = "black") ,
        axis.title.x = element_text(size = 4, color = "black") ,
        axis.title.y = element_text(size = 4, color = "black"),
        legend.text = element_text(size = 3, margin = margin(l = 2)),
        legend.title = element_text(size = 4, margin = margin(b = 2) ),
        legend.box.spacing = unit(1, "mm"),
        legend.key.spacing.y = unit(0, "mm"),
        legend.key.size = unit(2,"mm"),
        legend.spacing.x = unit(0.1, "mm"),
        panel.spacing=unit(0.2, "mm"))


# ggsave("Figures/CovariatesAllScenarios_Power_EffectSize_54min_v2.jpeg",
#        width = unit(4, "inches"),height = unit(2.5, "inches"), dpi = 600)

# Colour by natural parameter order type (Int v Slope) facet by Scenario
ggplot(full.covscen.nat %>% filter(order == "2nd" | CoefReg == "Slope")%>% filter(n.sites >33), 
       aes(x = log(n.sites), y = PWR, group = interaction(Parameter, Lik)))+
  geom_hline(yintercept = 0.95, col = "red", linetype = "dotted")+
  geom_line(aes( col = order), lwd = 0.8)+
  facet_grid(Scenario~Lik, labeller = labeller(Lik = lik_labels))+
  labs(x = "Number of Sites", y = "Power", col = "Effect Size")+ guides(alpha = "none")+
  coord_cartesian(ylim =  c(0, 1))+
  scale_x_continuous(breaks = log(nsites), labels = nsites)+
  scale_y_continuous(labels = scales::percent_format() )+
  scale_colour_manual(values =  safe_colorblind_palette[8:9] )+
  # scale_colour_viridis_c()+
  theme_bw()+
  theme(strip.text.y = element_text(size = 3.5, color = "black",  margin = margin(r = 2, l = 2)),
        strip.text.x = element_text(size = 3.5, color = "black", margin = margin(b = 2, t = 2)),
        axis.text.y = element_text(size = 3.5, color = "black"),
        axis.text.x =element_text(size = 3, color = "black") ,
        axis.title.x = element_text(size = 4, color = "black") ,
        axis.title.y = element_text(size = 4, color = "black"),
        legend.text = element_text(size = 3, margin = margin(l = 2)),
        legend.title = element_text(size = 4, margin = margin(b = 2) ),
        legend.box.spacing = unit(1, "mm"),
        legend.key.spacing.y = unit(0, "mm"),
        legend.key.size = unit(2,"mm"),
        legend.spacing.x = unit(0.1, "mm"),
        panel.spacing=unit(0.2, "mm"))


## Bias plots -----
### Natural parameters ----
ggplot(full.covscen.nat%>% filter(n.sites >33), 
       aes(x = log(n.sites), y = mu.p.bias, group = interaction(Parameter, Lik)))+
  geom_rect(aes(xmin = 2.5, xmax = 8.5, ymin = -0.05, ymax = 0.05), fill = "grey90",
            inherit.aes = F, alpha = 0.9)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "black", lwd = 0.15)+
  geom_line(aes(col = interaction(order)), lwd = 0.3)+
  # geom_point(size = 0.4,
  #            aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased"),
  #                col = interaction(order)))+
  scale_alpha_manual(values = c(1, 0.3))+
  facet_grid(Scenario~Lik, labeller = labeller(Lik = lik_labels))+
  labs(x = "Number of Sites", y = "Relative Bias", col = "Order")+ guides(alpha = "none")+
  coord_cartesian(ylim =  c(-.5, .5), xlim = c(4, 8))+
  scale_x_continuous(breaks = log(nsites), labels = nsites)+
  scale_y_continuous(labels = scales::percent_format() )+
  scale_colour_manual(values = cbbPalette[c(6:7)])+
  theme_bw()  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 3.5, color = "black",  margin = margin(r = 2, l = 2)),
        strip.text.x = element_text(size = 3.5, color = "black", margin = margin(b = 2, t = 2)),
        axis.text.y = element_text(size = 3.5, color = "black"),
        axis.text.x =element_text(size = 3, color = "black") ,
        axis.title.x = element_text(size = 4, color = "black") ,
        axis.title.y = element_text(size = 4, color = "black"),
        legend.text = element_text(size = 3, margin = margin(l = 2)),
        legend.title = element_text(size = 4, margin = margin(b = 2) ),
        legend.box.spacing = unit(1, "mm"),
        legend.key.spacing.y = unit(0, "mm"),
        legend.key.size = unit(2,"mm"),
        legend.spacing.x = unit(0.1, "mm"),
        panel.spacing=unit(0.2, "mm"))


# ggsave("Figures/CovariatesAllScenarios_NatPamRB_Order_54min_v2.jpeg",
#        width = unit(4, "inches"),height = unit(2.5, "inches"), dpi = 600)

## Marginal probs ----

ggplot(full.covscen.mar, 
       aes(x = log(n.sites), y = mu.p.bias,
           col = Species, group = interaction(Scenario, Species)))+
  geom_rect(aes(xmin = 2.5, xmax = 8.5, ymin = -0.05, ymax = 0.05), fill = "grey90",
            inherit.aes = F, alpha = 0.7)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "black")+
  geom_line(data = full.covscen.gen,
            aes(x = log(n.sites), y =mu.p.bias, group = Gen.Par), col = "grey65", lwd = 0.2)+
  geom_line(lwd = 0.3)+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  facet_grid(Scenario~Lik, labeller = labeller(Lik = lik_labels))+
  labs(x = "Number of Sites", y = "Mean Relative Bias (RB)", col = "Species")+
  scale_colour_manual(values= c(safe_colorblind_palette[-1]))+
  coord_cartesian(ylim =  c(-0.2, 0.2),  xlim = c(4, 8))+
  theme_bw()+guides(alpha = "none")+
  scale_x_continuous(breaks = log(nsites), labels = nsites)+
  scale_y_continuous(labels = scales::percent_format() )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 3.5, color = "black",  margin = margin(r = 2, l = 2)),
        strip.text.x = element_text(size = 3.5, color = "black", margin = margin(b = 2, t = 2)),
        axis.text.y = element_text(size = 3.5, color = "black"),
        axis.text.x =element_text(size = 3, color = "black") ,
        axis.title.x = element_text(size = 4, color = "black") ,
        axis.title.y = element_text(size = 4, color = "black"),
        legend.text = element_text(size = 3, margin = margin(l = 2)),
        legend.title = element_text(size = 4, margin = margin(b = 2) ),
        legend.box.spacing = unit(1, "mm"),
        legend.key.spacing.y = unit(0, "mm"),
        legend.key.size = unit(2,"mm"),
        legend.spacing.x = unit(0.1, "mm"),
        panel.spacing=unit(0.2, "mm"))

# ggsave("Figures/CovariatesAllScenarios_MargProbRB_54min_v2.jpeg",
#        width = unit(4, "inches"),height = unit(2.5, "inches"), dpi = 600)


# Let's make a fancy table ------
library(kableExtra)

colfunc<-colorRampPalette(c("red","pink"))
Col.pallete <- c(colfunc(5), "white") # n = number of steps (sequentially change)

table.color <- function(x){ case_when(
  abs(x) <= 5~ Col.pallete[6],
  abs(x) > 5 & abs(x) <=15 ~ Col.pallete[5],
  abs(x) > 15 & abs(x) <=35 ~ Col.pallete[4],
  abs(x) > 35 & abs(x) <=75 ~ Col.pallete[3],
  abs(x) > 75 & abs(x) <100 ~ Col.pallete[2],
  abs(x) >= 100 ~ Col.pallete[1]
)}

## Natural parameters ----
### one slope  -----
f1sl.nat.wide <- full.covscen.nat %>%
  filter(Scenario == "1 Slope" & n.sites >50) %>%
  select(n.sites, Parameter,  Lik, mu.p.bias) %>%
  mutate(mu.p.bias = round(mu.p.bias*100, 2),
         # Param.Type = factor(Param.Type,levels = c("Natural", "General")),
         Parameter = factor(Parameter, levels = c("f1.0",
                                                  "f2.0", 
                                                  "f12.0", "f12.1"))) %>%
  arrange(n.sites,Parameter, Lik) %>%
  pivot_wider(names_from = c(Parameter), values_from = mu.p.bias)%>%
  arrange(Lik)


kable.1sl <- kableExtra::kbl(f1sl.nat.wide, booktabs = T, 
                             col.names = c("N","Likelihood",
                                           colnames(f1sl.nat.wide)[-c(1:2)]),
                             digits = 1,align = "c", 
                             linesep = "", format = "latex") %>%
  add_header_above(c(" ", "", "1st Order" = 2, "2nd Order" = 2)) %>%
  # kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 2,latex_hline = "major", 
                row_group_label_position = "first") 



# Loop over columns for column_spec
for (col_num in 3:ncol(f1sl.nat.wide)) {
  
  kable.1sl <- kableExtra::column_spec(kable.1sl,
                                       col_num,
                                       background =  unlist(lapply(f1sl.nat.wide[[col_num]],function(x) table.color(x) )))
}

kable.1sl

### Three slopes  -----
f3sl.nat.wide <- full.covscen.nat %>%
  filter(Scenario == "3 Slopes" & n.sites >50) %>%
  select(n.sites, Parameter,  Lik, mu.p.bias) %>%
  mutate(mu.p.bias = round(mu.p.bias*100, 2),
         # Param.Type = factor(Param.Type,levels = c("Natural", "General")),
         Parameter = factor(Parameter, levels = c("f1.0", "f1.1",
                                                  "f2.0", "f2.1", 
                                                  "f12.0", "f12.1"))) %>%
  arrange(n.sites,Parameter, Lik) %>%
  pivot_wider(names_from = c(Parameter), values_from = mu.p.bias)%>%
  arrange(Lik)


kable.3sl <- kableExtra::kbl(f3sl.nat.wide, booktabs = T, 
                             col.names = c("N","Likelihood",
                                           colnames(f3sl.nat.wide)[-c(1:2)]),
                             digits = 1,align = "c", 
                             linesep = "", format = "latex") %>%
  add_header_above(c(" ", "", "1st Order" = 4, "2nd Order" = 2)) %>%
  # kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 2,latex_hline = "major", 
                row_group_label_position = "first") 



# Loop over columns for column_spec
for (col_num in 3:ncol(f3sl.nat.wide)) {
  
  kable.3sl <- kableExtra::column_spec(kable.3sl,
                                       col_num,
                                       background =  unlist(lapply(f3sl.nat.wide[[col_num]],function(x) table.color(x) )))
}

kable.3sl


### Four slopes  -----
f4sl.nat.wide <- full.covscen.nat %>%
  filter(Scenario == "4 Slopes" & n.sites >50) %>%
  select(n.sites, Parameter,  Lik, mu.p.bias) %>%
  mutate(mu.p.bias = round(mu.p.bias*100, 2),
         # Param.Type = factor(Param.Type,levels = c("Natural", "General")),
         Parameter = factor(Parameter, levels = c("f1.0", "f1.1",
                                                  "f2.0", "f2.1", 
                                                  "f12.0", "f12.1", "f12.2"))) %>%
  arrange(n.sites,Parameter, Lik) %>%
  pivot_wider(names_from = c(Parameter), values_from = mu.p.bias)%>%
  arrange(Lik)


kable.4sl <- kableExtra::kbl(f4sl.nat.wide, booktabs = T, 
                             col.names = c("N","Likelihood",
                                           colnames(f4sl.nat.wide)[-c(1:2)]),
                             digits = 1,align = "c", 
                             linesep = "", format = "latex") %>%
  add_header_above(c(" ", "", "1st Order" = 4, "2nd Order" = 3)) %>%
  # kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 2,latex_hline = "major", 
                row_group_label_position = "first") 



# Loop over columns for column_spec
for (col_num in 3:ncol(f4sl.nat.wide)) {
  
  kable.4sl <- kableExtra::column_spec(kable.4sl,
                                       col_num,
                                       background =  unlist(lapply(f4sl.nat.wide[[col_num]],function(x) table.color(x) )))
}

kable.4sl


### Five slopes  -----
f5sl.nat.wide <- full.covscen.nat %>%
  filter(Scenario == "5 Slopes" & n.sites >50) %>%
  select(n.sites, Parameter,  Lik, mu.p.bias) %>%
  mutate(mu.p.bias = round(mu.p.bias*100, 2),
         # Param.Type = factor(Param.Type,levels = c("Natural", "General")),
         Parameter = factor(Parameter, levels = c("f1.0", "f1.1", "f1.2",
                                                  "f2.0", "f2.1", 
                                                  "f12.0", "f12.1", "f12.2"))) %>%
  arrange(n.sites,Parameter, Lik) %>%
  pivot_wider(names_from = c(Parameter), values_from = mu.p.bias)%>%
  arrange(Lik)


kable.5sl <- kableExtra::kbl(f5sl.nat.wide, booktabs = T, 
                             col.names = c("N","Likelihood",
                                           colnames(f5sl.nat.wide)[-c(1:2)]),
                             digits = 1,align = "c", 
                             linesep = "", format = "latex") %>%
  add_header_above(c(" ", "", "1st Order" = 4, "2nd Order" = 4)) %>%
  # kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 2,latex_hline = "major", 
                row_group_label_position = "first") 



# Loop over columns for column_spec
for (col_num in 3:ncol(f5sl.nat.wide)) {
  
  kable.5sl <- kableExtra::column_spec(kable.5sl,
                                       col_num,
                                       background =  unlist(lapply(f5sl.nat.wide[[col_num]],function(x) table.color(x) )))
}

kable.5sl



# Convergence rates ----

full.scen.nat <- read.csv("Results/NullAllScenNat.csv")
full.sp.nat   <- read.csv("Results/NetworkAllScenNat.csv")
shared.vars   <- intersect(colnames(full.scen.nat), colnames(full.covscen.nat))

all.scen.nat <- rbind(full.scen.nat[shared.vars], full.sp.nat[shared.vars], full.covscen.nat[shared.vars])


conv.sum <- all.scen.nat %>%
  mutate(Scenario = factor(Scenario, 
                           levels = c("Scenario1",
                                      "Scenario2",
                                      "Scenario3",
                                      "Scenario4",
                                      "Scenario5",
                                      "Scenario6",
                                      "3 Species",
                                      "4 Species",
                                      "5 Species",
                                      "1 Slope",
                                      "3 Slopes",
                                      "4 Slopes",
                                      "5 Slopes"))) %>%
  select(Scenario, Lik,n.sites, ndatasets) %>%
  group_by(n.sites, Scenario, Lik) %>%
  summarise(Conv.rate = unique(ndatasets)/100) %>%
  ungroup()%>%
  arrange(Lik, n.sites) %>%
  pivot_wider(values_from = Conv.rate, names_from = c(Scenario), values_fill = 1) %>%
  rename(N = n.sites, Likelihood = Lik)

conv.sum[-c(1,2)] <- apply(conv.sum[-c(1,2)], 2 ,function(a) format(round(a, digits=2), nsmall = 2) )


kableExtra::kbl(conv.sum, booktabs = T, format = "latex", 
                linesep = "", align = "c", digits = 2,
                col.names = c("N", "Likelihood", 
                              1:6, paste0(3:5 , "Spp."),paste0(c(1, 3:5), "Slopes") ))%>%
  add_header_above(c(" ", "", "Null models" = 6, "Network models" = 3, "Covariate models" = 4)) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 2,latex_hline = "major", 
                row_group_label_position = "first") 

