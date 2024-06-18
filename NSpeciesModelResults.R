##############################################
##### 3+Species model results processing #####
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

# function for retreiving general params
gen.param.fun <- function(beta, nspecies = 2, 
                          occ_formulas = NULL, occ_covs = NULL, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  if(is.null(occ_formulas)){
    occ_formulas <- rep("~1", nspecies + nspecies*(nspecies-1)/2)
    nocc_covs = 0} else{
    nocc_covs <- max(as.numeric(gsub(pattern = "[a-z]", "",
                                     str_extract(occ_formulas, "cov[0-9]"))), na.rm = T)
    
    if(is.infinite(nocc_covs)) nocc_covs = 0
  }
  
  
  if(is.null(occ_covs) & nocc_covs !=0) {
    occ_covs <- data.frame(matrix(rep(seq(-1, 1, 
                                          length.out  = 100), nocc_covs),
                                  ncol = nocc_covs))
    names(occ_covs) <- paste('occ_cov',1:nocc_covs,sep='')
  } 
  
  if(nocc_covs !=0) names(occ_covs) <- paste('occ_cov',1:ncol(occ_covs),sep='')
  n.nat.params <- nspecies + nspecies*(nspecies-1)/2
  nat.params <- paste0("f", c(1:nspecies, 
                              apply(combn(1:nspecies, 2),
                                    2,paste0, collapse = "" )))
  
  
  if(is.null(occ_covs))  occ_covs = data.frame(X = 1)
  
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
      N <- ifelse(nocc_covs != 0, nrow(occ_covs),1)
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



# 3 Species Scenario  ----

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


## Read object
scen.3sp <- read_rds("Results/Scenario1_3spNullseed1337_simResults_v2.rds")


## Normal likelihood 
State_sim.3sp <- scen.3sp$State.params %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)) |any(is.na(SE)) , 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam))
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

scen.3sp.conv <-State_sim.3sp %>%
  group_by(n.sites, Parameter) %>%
  summarise(Conv.rate = mean(conv)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim.3sp.sum <- State_sim.3sp %>%
  group_by(n.sites, Parameter) %>%
  filter(conv == 1) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias), mu.CV = mean(CV, na.rm = T),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered, na.rm = T),                    # Coverage rate 
            PWR = mean(below.alpha, na.rm = T),               # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate (should be 100%)
            ndatasets = n(),                                  # Number datasets in summary
            Lik = "LL") %>% 
  ungroup()


#### General parameters (psi11, psi10, psi01, psi00)

GenParam_sim3sp <- State_sim.3sp %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = as.list(Estimate), 
                                                nspecies = nspecies,
                                                occ_formulas = occ_formulas))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               nspecies = nspecies,
                                               occ_formulas = occ_formulas))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             nspecies = nspecies,
                                             occ_formulas = occ_formulas))) %>%
  ungroup


GenParam_sim3sp <-  GenParam_sim3sp %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim3sp.sum <- GenParam_sim3sp %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()


#### Derived Parameters 


Marg.prob.3sp <-  GenParam_sim3sp %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2", "Sp3"),
            Marg.est = c(sum(psi.est[substr(Gen.Par,1,1) ==1]),
                         sum(psi.est[substr(Gen.Par,2,2) ==1]),
                         sum(psi.est[substr(Gen.Par,3,3) ==1])),
            Marg.og = c(sum(og.psi[substr(Gen.Par,1,1) ==1]),
                        sum(og.psi[substr(Gen.Par,2,2) ==1]),
                        sum(og.psi[substr(Gen.Par,3,3) ==1])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()




State_sim.3sp$Parameter == "f1.0"


Cond.prob.3sp <-  State_sim.3sp %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1",
                          "Sp1=1|Sp3=0", "Sp1=1|Sp3=1",
                          "Sp2=1|Sp1=0", "Sp2=1|Sp1=1",
                          "Sp2=1|Sp3=0", "Sp2=1|Sp3=1",
                          "Sp3=1|Sp1=0", "Sp3=1|Sp1=1",
                          "Sp3=1|Sp2=0", "Sp3=1|Sp2=1"),
            # Estimated
            Cond.est = c(plogis(Estimate[Parameter == "f1.0"]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(Estimate[Parameter %in% c("f1.0")]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f23.0")]))),
            # Original
            Cond.og =  c(plogis(og_param.val[Parameter == "f1.0"]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(og_param.val[Parameter %in% c("f1.0")]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f23.0")]))),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()


## Summaries of derived

Marg.prob3sp.sum <- Marg.prob.3sp %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,  ndatasets = n(),               # Upper bias CI
            Lik = "LL") %>% 
  ungroup()


Cond.prob3sp.sum <- Cond.prob.3sp %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,    ndatasets = n(),             # Upper bias CI
            Lik = "LL") %>% 
  ungroup()



## Penalised likelihood  -----
State_sim.3sp.pl <- scen.3sp$State.params.pl %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)) |any(is.na(SE)) , 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam))
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

scen.3sp.pl.conv <-State_sim.3sp.pl %>%
  group_by(n.sites, Parameter) %>%
  summarise(Conv.rate = mean(conv)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim.3sp.pl.sum <- State_sim.3sp.pl %>%
  group_by(n.sites, Parameter) %>%
  filter(conv == 1) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias), mu.CV = mean(CV, na.rm = T),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered, na.rm = T),                    # Coverage rate 
            PWR = mean(below.alpha, na.rm = T),               # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate (should be 100%)
            ndatasets = n(),                                  # Number datasets in summary
            Lik = "PL") %>% 
  ungroup()


#### General parameters (psi11, psi10, psi01, psi00)

GenParam_sim3sp.pl <- State_sim.3sp.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = as.list(Estimate), 
                                                nspecies = nspecies,
                                                occ_formulas = occ_formulas))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               nspecies = nspecies,
                                               occ_formulas = occ_formulas))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             nspecies = nspecies,
                                             occ_formulas = occ_formulas))) %>%
  ungroup


GenParam_sim3sp.pl <-  GenParam_sim3sp.pl %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim3sp.pl.sum <- GenParam_sim3sp.pl %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()


#### Derived Parameters 


Marg.prob.3sp.pl <-  GenParam_sim3sp.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2", "Sp3"),
            Marg.est = c(sum(psi.est[substr(Gen.Par,1,1) ==1]),
                         sum(psi.est[substr(Gen.Par,2,2) ==1]),
                         sum(psi.est[substr(Gen.Par,3,3) ==1])),
            Marg.og = c(sum(og.psi[substr(Gen.Par,1,1) ==1]),
                        sum(og.psi[substr(Gen.Par,2,2) ==1]),
                        sum(og.psi[substr(Gen.Par,3,3) ==1])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()




State_sim.3sp.pl$Parameter == "f1.0"


Cond.prob.3sp.pl <-  State_sim.3sp.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1",
                          "Sp1=1|Sp3=0", "Sp1=1|Sp3=1",
                          "Sp2=1|Sp1=0", "Sp2=1|Sp1=1",
                          "Sp2=1|Sp3=0", "Sp2=1|Sp3=1",
                          "Sp3=1|Sp1=0", "Sp3=1|Sp1=1",
                          "Sp3=1|Sp2=0", "Sp3=1|Sp2=1"),
            # Estimated
            Cond.est = c(plogis(Estimate[Parameter == "f1.0"]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(Estimate[Parameter %in% c("f1.0")]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f23.0")]))),
            # Original
            Cond.og =  c(plogis(og_param.val[Parameter == "f1.0"]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(og_param.val[Parameter %in% c("f1.0")]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f23.0")]))),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()


## Summaries of derived

Marg.prob3sp.pl.sum <- Marg.prob.3sp.pl %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,  ndatasets = n(),               # Upper bias CI
            Lik = "PL") %>% 
  ungroup()


Cond.prob3sp.pl.sum <- Cond.prob.3sp.pl %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,    ndatasets = n(),             # Upper bias CI
            Lik = "PL") %>% 
  ungroup()


###### Combining LL and PL into one for plotting

## Natural parameters
full.scen3sp.nat <- rbind(State_sim.3sp.sum, State_sim.3sp.pl.sum)

## General parameters

full.scen3sp.gen <- rbind(GenParam_sim3sp.sum, GenParam_sim3sp.pl.sum)

## Marginal probabilities
full.scen3sp.mar <- rbind(Marg.prob3sp.sum, Marg.prob3sp.pl.sum)

## Conditional probabilities
full.scen3sp.con <- rbind(Cond.prob3sp.sum, Cond.prob3sp.pl.sum)







# 4 Species Scenario  ----


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


## Read object
scen.4sp <- read_rds("Results/Scenario1_4spNullseed1337_simResults_v2.rds")


## Normal likelihood 
State_sim.4sp <- scen.4sp$State.params %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)) |any(is.na(SE)) , 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam))
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

scen.4sp.conv <-State_sim.4sp %>%
  group_by(n.sites, Parameter) %>%
  summarise(Conv.rate = mean(conv)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim.4sp.sum <- State_sim.4sp %>%
  group_by(n.sites, Parameter) %>%
  filter(conv == 1) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias), mu.CV = mean(CV, na.rm = T),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered, na.rm = T),                    # Coverage rate 
            PWR = mean(below.alpha, na.rm = T),               # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate (should be 100%)
            ndatasets = n(),                                  # Number datasets in summary
            Lik = "LL") %>% 
  ungroup()


#### General parameters (psi11, psi10, psi01, psi00)

GenParam_sim4sp <- State_sim.4sp %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = as.list(Estimate), 
                                                nspecies = nspecies,
                                                occ_formulas = occ_formulas))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               nspecies = nspecies,
                                               occ_formulas = occ_formulas))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             nspecies = nspecies,
                                             occ_formulas = occ_formulas))) %>%
  ungroup


GenParam_sim4sp <-  GenParam_sim4sp %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim4sp.sum <- GenParam_sim4sp %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()


#### Derived Parameters 




Marg.prob.4sp <-  GenParam_sim4sp %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2", "Sp3", "Sp4"),
            Marg.est = c(sum(psi.est[substr(Gen.Par,1,1) ==1]),
                         sum(psi.est[substr(Gen.Par,2,2) ==1]),
                         sum(psi.est[substr(Gen.Par,3,3) ==1]),
                         sum(psi.est[substr(Gen.Par,4,4) ==1])),
            Marg.og = c(sum(og.psi[substr(Gen.Par,1,1) ==1]),
                        sum(og.psi[substr(Gen.Par,2,2) ==1]),
                        sum(og.psi[substr(Gen.Par,3,3) ==1]),
                        sum(og.psi[substr(Gen.Par,4,4) ==1])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()


State_sim.4sp %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(ok =  c(plogis(Estimate[Parameter == "f1.0"]),
                    plogis(sum(Estimate[Parameter %in% c("f1.0", "f12.0")])),
                    plogis(Estimate[Parameter %in% c("f1.0")]),
                    plogis(sum(Estimate[Parameter %in% c("f1.0", "f13.0")])),
                    plogis(Estimate[Parameter %in% c("f1.0")]),
                    plogis(sum(Estimate[Parameter %in% c("f1.0", "f14.0")]))),
            )


Cond.prob.4sp <-  State_sim.4sp %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1",
                          "Sp1=1|Sp3=0", "Sp1=1|Sp3=1",
                          "Sp1=1|Sp4=0", "Sp1=1|Sp4=1",
                          "Sp2=1|Sp1=0", "Sp2=1|Sp1=1",
                          "Sp2=1|Sp3=0", "Sp2=1|Sp3=1",
                          "Sp2=1|Sp4=0", "Sp2=1|Sp4=1",
                          "Sp3=1|Sp1=0", "Sp3=1|Sp1=1",
                          "Sp3=1|Sp2=0", "Sp3=1|Sp2=1",
                          "Sp3=1|Sp4=0", "Sp3=1|Sp4=1",
                          "Sp4=1|Sp1=0", "Sp4=1|Sp1=1",
                          "Sp4=1|Sp2=0", "Sp4=1|Sp2=1",
                          "Sp4=1|Sp3=0", "Sp4=1|Sp3=1"),
            # Estimated
            Cond.est = c(plogis(Estimate[Parameter == "f1.0"]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(Estimate[Parameter %in% c("f1.0")]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(Estimate[Parameter %in% c("f1.0")]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f14.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f24.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f23.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f14.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f14.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f24.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f34.0")]))),
            # Original
            Cond.og =  c(plogis(og_param.val[Parameter == "f1.0"]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(og_param.val[Parameter %in% c("f1.0")]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(og_param.val[Parameter %in% c("f1.0")]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f14.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f24.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f23.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f34.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f14.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f24.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f34.0")]))
            ),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()


## Summaries of derived

Marg.prob4sp.sum <- Marg.prob.4sp %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,  ndatasets = n(),               # Upper bias CI
            Lik = "LL") %>% 
  ungroup()


Cond.prob4sp.sum <- Cond.prob.4sp %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,    ndatasets = n(),             # Upper bias CI
            Lik = "LL") %>% 
  ungroup()



## Penalised likelihood  -----
State_sim.4sp.pl <- scen.4sp$State.params.pl %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)) |any(is.na(SE)) , 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam))
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

scen.4sp.pl.conv <-State_sim.4sp.pl %>%
  group_by(n.sites) %>%
  summarise(Conv.rate = mean(conv)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim.4sp.pl.sum <- State_sim.4sp.pl %>%
  group_by(n.sites, Parameter) %>%
  filter(conv == 1) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias), mu.CV = mean(CV, na.rm = T),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered, na.rm = T),                    # Coverage rate 
            PWR = mean(below.alpha, na.rm = T),               # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate (should be 100%)
            ndatasets = n(),                                  # Number datasets in summary
            Lik = "PL") %>% 
  ungroup()


#### General parameters (psi11, psi10, psi01, psi00)

GenParam_sim4sp.pl <- State_sim.4sp.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = as.list(Estimate), 
                                                nspecies = nspecies,
                                                occ_formulas = occ_formulas))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               nspecies = nspecies,
                                               occ_formulas = occ_formulas))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             nspecies = nspecies,
                                             occ_formulas = occ_formulas))) %>%
  ungroup


GenParam_sim4sp.pl <-  GenParam_sim4sp.pl %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim4sp.pl.sum <- GenParam_sim4sp.pl %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()


ggplot(GenParam_sim4sp.pl.sum)+
  geom_hline(yintercept = 0.05, col = "red", linetype = "dashed")+
  geom_hline(yintercept = -0.05, col = "red", linetype = "dashed")+
  geom_line(aes(x = log(n.sites), y = mu.p.bias, col = Gen.Par))+theme_bw()

#### Derived Parameters 




Marg.prob.4sp.pl <-  GenParam_sim4sp.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2", "Sp3", "Sp4"),
            Marg.est = c(sum(psi.est[substr(Gen.Par,1,1) ==1]),
                         sum(psi.est[substr(Gen.Par,2,2) ==1]),
                         sum(psi.est[substr(Gen.Par,3,3) ==1]),
                         sum(psi.est[substr(Gen.Par,4,4) ==1])),
            Marg.og = c(sum(og.psi[substr(Gen.Par,1,1) ==1]),
                        sum(og.psi[substr(Gen.Par,2,2) ==1]),
                        sum(og.psi[substr(Gen.Par,3,3) ==1]),
                        sum(og.psi[substr(Gen.Par,4,4) ==1])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()





Cond.prob.4sp.pl <-  State_sim.4sp.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1",
                          "Sp1=1|Sp3=0", "Sp1=1|Sp3=1",
                          "Sp1=1|Sp4=0", "Sp1=1|Sp4=1",
                          "Sp2=1|Sp1=0", "Sp2=1|Sp1=1",
                          "Sp2=1|Sp3=0", "Sp2=1|Sp3=1",
                          "Sp2=1|Sp4=0", "Sp2=1|Sp4=1",
                          "Sp3=1|Sp1=0", "Sp3=1|Sp1=1",
                          "Sp3=1|Sp2=0", "Sp3=1|Sp2=1",
                          "Sp3=1|Sp4=0", "Sp3=1|Sp4=1",
                          "Sp4=1|Sp1=0", "Sp4=1|Sp1=1",
                          "Sp4=1|Sp2=0", "Sp4=1|Sp2=1",
                          "Sp4=1|Sp3=0", "Sp4=1|Sp3=1"),
            # Estimated
            Cond.est = c(plogis(Estimate[Parameter == "f1.0"]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(Estimate[Parameter %in% c("f1.0")]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(Estimate[Parameter %in% c("f1.0")]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f14.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f24.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f23.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f14.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f14.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f24.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f34.0")]))),
            # Original
            Cond.og =  c(plogis(og_param.val[Parameter == "f1.0"]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(og_param.val[Parameter %in% c("f1.0")]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(og_param.val[Parameter %in% c("f1.0")]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f14.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f24.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f23.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f14.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f14.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f24.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f34.0")]))
            ),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()


## Summaries of derived

Marg.prob4sp.pl.sum <- Marg.prob.4sp.pl %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,  ndatasets = n(),               # Upper bias CI
            Lik = "PL") %>% 
  ungroup()


ggplot(Marg.prob4sp.pl.sum)+
  geom_hline(yintercept = 0.05, col = "red", linetype = "dashed")+
  geom_hline(yintercept = -0.05, col = "red", linetype = "dashed")+
  geom_line(aes(x = log(n.sites), y = mu.p.bias, col = Species))+theme_bw()



Cond.prob4sp.pl.sum <- Cond.prob.4sp.pl %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,  
            ndatasets = n(),             # Upper bias CI
            Lik = "PL") %>% 
  ungroup()




###### Combining LL and PL into one for plotting

## Natural parameters
full.scen4sp.nat <- rbind(State_sim.4sp.sum, State_sim.4sp.pl.sum)

full.scen4sp.nat <- full.scen4sp.nat %>%
  mutate(order = ifelse(nchar(Parameter)>4, "2nd", "1st"))

## General parameters

full.scen4sp.gen <- rbind(GenParam_sim4sp.sum, GenParam_sim4sp.pl.sum)

## Marginal probabilities
full.scen4sp.mar <- rbind(Marg.prob4sp.sum, Marg.prob4sp.pl.sum)

## Conditional probabilities
full.scen4sp.con <- rbind(Cond.prob4sp.sum, Cond.prob4sp.pl.sum)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

colour2 <- RColorBrewer::brewer.pal(10, "Paired")


ggplot(full.scen4sp.nat, aes(x = log(n.sites), y = mu.p.bias))+
  geom_line( aes( linetype = Lik, col = Parameter), lwd = 1)+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased"), col = Parameter))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  # facet_wrap(~Parameter)+
  geom_hline(yintercept = 0.05, col = "red", linetype = "dashed")+
  geom_hline(yintercept = -0.05, col = "red", linetype = "dashed")+
  labs(x = "log(Number of Sites)", y = "Mean Relative Bias", title = "4 Species", col = "Par")+
  scale_colour_manual(values= c(colour2))+
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")+
  coord_cartesian(ylim =  c(-2.5, 2.5))+
  guides(alpha = "none")+
  facet_grid(~order)






# 5 Species Scenario  ----


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

## Read object
scen.5sp <- read_rds("Results/Scenario1_5spNullseed1337_simResults_v2.rds")


## Normal likelihood 
State_sim.5sp <- scen.5sp$State.params %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)) |any(is.na(SE)) , 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam))
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

scen.5sp.conv <-State_sim.5sp %>%
  group_by(n.sites, Parameter) %>%
  summarise(Conv.rate = mean(conv)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim.5sp.sum <- State_sim.5sp %>%
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

GenParam_sim5sp <- State_sim.5sp %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = as.list(Estimate), 
                                                nspecies = nspecies,
                                                occ_formulas = occ_formulas))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               nspecies = nspecies,
                                               occ_formulas = occ_formulas))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             nspecies = nspecies,
                                             occ_formulas = occ_formulas))) %>%
  ungroup


GenParam_sim5sp <-  GenParam_sim5sp %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim5sp.sum <- GenParam_sim5sp %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "LL") %>% 
  ungroup()


#### Derived Parameters 




Marg.prob.5sp <-  GenParam_sim5sp %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2", "Sp3", "Sp4", "Sp5"),
            Marg.est = c(sum(psi.est[substr(Gen.Par,1,1) ==1]),
                         sum(psi.est[substr(Gen.Par,2,2) ==1]),
                         sum(psi.est[substr(Gen.Par,3,3) ==1]),
                         sum(psi.est[substr(Gen.Par,4,4) ==1]),
                         sum(psi.est[substr(Gen.Par,5,5) ==1])),
            Marg.og = c(sum(og.psi[substr(Gen.Par,1,1) ==1]),
                        sum(og.psi[substr(Gen.Par,2,2) ==1]),
                        sum(og.psi[substr(Gen.Par,3,3) ==1]),
                        sum(og.psi[substr(Gen.Par,4,4) ==1]),
                        sum(og.psi[substr(Gen.Par,5,5) ==1])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()



Cond.prob.5sp <-  State_sim.5sp %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1",
                          "Sp1=1|Sp3=0", "Sp1=1|Sp3=1",
                          "Sp1=1|Sp4=0", "Sp1=1|Sp4=1",
                          "Sp1=1|Sp5=0", "Sp1=1|Sp5=1",
                          "Sp2=1|Sp1=0", "Sp2=1|Sp1=1",
                          "Sp2=1|Sp3=0", "Sp2=1|Sp3=1",
                          "Sp2=1|Sp4=0", "Sp2=1|Sp4=1",
                          "Sp2=1|Sp5=0", "Sp2=1|Sp5=1",
                          "Sp3=1|Sp1=0", "Sp3=1|Sp1=1",
                          "Sp3=1|Sp2=0", "Sp3=1|Sp2=1",
                          "Sp3=1|Sp4=0", "Sp3=1|Sp4=1",
                          "Sp3=1|Sp5=0", "Sp3=1|Sp5=1",
                          "Sp4=1|Sp1=0", "Sp4=1|Sp1=1",
                          "Sp4=1|Sp2=0", "Sp4=1|Sp2=1",
                          "Sp4=1|Sp3=0", "Sp4=1|Sp3=1",
                          "Sp4=1|Sp5=0", "Sp4=1|Sp5=1",
                          "Sp5=1|Sp1=0", "Sp5=1|Sp1=1",
                          "Sp5=1|Sp2=0", "Sp5=1|Sp2=1",
                          "Sp5=1|Sp3=0", "Sp5=1|Sp3=1",
                          "Sp5=1|Sp4=0", "Sp5=1|Sp4=1"),
            # Estimated
            Cond.est = c(plogis(Estimate[Parameter == "f1.0"]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(Estimate[Parameter %in% c("f1.0")]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(Estimate[Parameter %in% c("f1.0")]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f14.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f1.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f15.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f24.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f25.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f23.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f34.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f35.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f14.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f24.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f34.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f45.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0", "f15.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0", "f25.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0", "f35.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0", "f45.0")]))),
            # Original
            Cond.og =  c(plogis(og_param.val[Parameter == "f1.0"]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(og_param.val[Parameter %in% c("f1.0")]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(og_param.val[Parameter %in% c("f1.0")]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f14.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f15.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f24.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f25.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f23.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f34.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f35.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f14.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f24.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f34.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f45.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0", "f15.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0", "f25.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0", "f35.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0", "f45.0")]))),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()


## Summaries of derived

Marg.prob5sp.sum <- Marg.prob.5sp %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,  ndatasets = n(),               # Upper bias CI
            Lik = "LL") %>% 
  ungroup()


Cond.prob5sp.sum <- Cond.prob.5sp %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,    ndatasets = n(),             # Upper bias CI
            Lik = "LL") %>% 
  ungroup()



# Penalised likelihood  ----

## 
State_sim.5sp.pl <- scen.5sp$State.params.pl %>%
  group_by(n.sites, niter) %>%
  mutate(conv = ifelse(any(is.nan(SE)) |any(is.na(SE)) , 0, conv)) %>%
  ungroup %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam))
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

scen.5sp.pl.conv <-State_sim.5sp.pl %>%
  group_by(n.sites, Parameter) %>%
  summarise(Conv.rate = mean(conv)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim.5sp.pl.sum <- State_sim.5sp.pl %>%
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
            Lik = "PL") %>% 
  ungroup()


#### General parameters (psi11, psi10, psi01, psi00)

GenParam_sim5sp.pl <- State_sim.5sp.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(gen.param.fun(beta = as.list(Estimate), 
                                                nspecies = nspecies,
                                                occ_formulas = occ_formulas))),
            og.psi = as.vector(t(gen.param.fun(beta = beta, 
                                               nspecies = nspecies,
                                               occ_formulas = occ_formulas))),
            Gen.Par = colnames(gen.param.fun(beta = beta, 
                                             nspecies = nspecies,
                                             occ_formulas = occ_formulas))) %>%
  ungroup


GenParam_sim5sp.pl <-  GenParam_sim5sp.pl %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim5sp.pl.sum <- GenParam_sim5sp.pl %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            ndatasets = n(),
            Lik = "PL") %>% 
  ungroup()


#### Derived Parameters 




Marg.prob.5sp.pl <-  GenParam_sim5sp.pl %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2", "Sp3", "Sp4", "Sp5"),
            Marg.est = c(sum(psi.est[substr(Gen.Par,1,1) ==1]),
                         sum(psi.est[substr(Gen.Par,2,2) ==1]),
                         sum(psi.est[substr(Gen.Par,3,3) ==1]),
                         sum(psi.est[substr(Gen.Par,4,4) ==1]),
                         sum(psi.est[substr(Gen.Par,5,5) ==1])),
            Marg.og = c(sum(og.psi[substr(Gen.Par,1,1) ==1]),
                        sum(og.psi[substr(Gen.Par,2,2) ==1]),
                        sum(og.psi[substr(Gen.Par,3,3) ==1]),
                        sum(og.psi[substr(Gen.Par,4,4) ==1]),
                        sum(og.psi[substr(Gen.Par,5,5) ==1])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()



Cond.prob.5sp.pl <-  State_sim.5sp.pl %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(Cond.prob = c("Sp1=1|Sp2=0", "Sp1=1|Sp2=1",
                          "Sp1=1|Sp3=0", "Sp1=1|Sp3=1",
                          "Sp1=1|Sp4=0", "Sp1=1|Sp4=1",
                          "Sp1=1|Sp5=0", "Sp1=1|Sp5=1",
                          "Sp2=1|Sp1=0", "Sp2=1|Sp1=1",
                          "Sp2=1|Sp3=0", "Sp2=1|Sp3=1",
                          "Sp2=1|Sp4=0", "Sp2=1|Sp4=1",
                          "Sp2=1|Sp5=0", "Sp2=1|Sp5=1",
                          "Sp3=1|Sp1=0", "Sp3=1|Sp1=1",
                          "Sp3=1|Sp2=0", "Sp3=1|Sp2=1",
                          "Sp3=1|Sp4=0", "Sp3=1|Sp4=1",
                          "Sp3=1|Sp5=0", "Sp3=1|Sp5=1",
                          "Sp4=1|Sp1=0", "Sp4=1|Sp1=1",
                          "Sp4=1|Sp2=0", "Sp4=1|Sp2=1",
                          "Sp4=1|Sp3=0", "Sp4=1|Sp3=1",
                          "Sp4=1|Sp5=0", "Sp4=1|Sp5=1",
                          "Sp5=1|Sp1=0", "Sp5=1|Sp1=1",
                          "Sp5=1|Sp2=0", "Sp5=1|Sp2=1",
                          "Sp5=1|Sp3=0", "Sp5=1|Sp3=1",
                          "Sp5=1|Sp4=0", "Sp5=1|Sp4=1"),
            # Estimated
            Cond.est = c(plogis(Estimate[Parameter == "f1.0"]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(Estimate[Parameter %in% c("f1.0")]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(Estimate[Parameter %in% c("f1.0")]),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f14.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f1.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f1.0", "f15.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f24.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f2.0", "f25.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f23.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f34.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f3.0", "f35.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f14.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f24.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f34.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f4.0", "f45.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0", "f15.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0", "f25.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0", "f35.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0")])),
                         plogis(sum(Estimate[Parameter %in% c("f5.0", "f45.0")]))),
            # Original
            Cond.og =  c(plogis(og_param.val[Parameter == "f1.0"]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f12.0")])),
                         plogis(og_param.val[Parameter %in% c("f1.0")]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f13.0")])),
                         plogis(og_param.val[Parameter %in% c("f1.0")]),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f14.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f1.0", "f15.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f12.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f23.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f24.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f2.0", "f25.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f13.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f23.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f34.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f3.0", "f35.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f14.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f24.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f34.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f4.0", "f45.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0", "f15.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0", "f25.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0", "f35.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0")])),
                         plogis(sum(og_param.val[Parameter %in% c("f5.0", "f45.0")]))),
            bias = Cond.est - Cond.og, prop.bias = (Cond.est - Cond.og)/Cond.og) %>%
  ungroup()


## Summaries of derived

Marg.prob5sp.pl.sum <- Marg.prob.5sp.pl %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,  ndatasets = n(),               # Upper bias CI
            Lik = "PL") %>% 
  ungroup()


Cond.prob5sp.pl.sum <- Cond.prob.5sp.pl %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,    ndatasets = n(),             # Upper bias CI
            Lik = "PL") %>% 
  ungroup()




###### Combining LL and PL into one for plotting

## Natural parameters
full.scen5sp.nat <- rbind(State_sim.5sp.sum, State_sim.5sp.pl.sum)

full.scen5sp.nat <- full.scen5sp.nat %>%
  mutate(order = ifelse(nchar(Parameter)>4, "2nd", "1st"))

## General parameters

full.scen5sp.gen <- rbind(GenParam_sim5sp.sum, GenParam_sim5sp.pl.sum)

## Marginal probabilities
full.scen5sp.mar <- rbind(Marg.prob5sp.sum, Marg.prob5sp.pl.sum)

ggplot(full.scen5sp.mar)+
  geom_hline(yintercept = 0.05, col = "red", linetype = "dotted")+
  geom_hline(yintercept = -0.05, col = "red", linetype = "dotted")+
  geom_line(aes(x = log(n.sites), y = mu.p.bias, col = Species, linetype = Lik), lwd = 0.8)+
  theme_bw()
## Conditional probabilities
full.scen5sp.con <- rbind(Cond.prob5sp.sum, Cond.prob5sp.pl.sum)



## Power plot for 5 species

ggplot(full.scen5sp.nat %>% filter((order == "2nd")))+
  geom_line(aes(x = log(n.sites), y = PWR, col = og.val, linetype = Lik))+
  facet_wrap(~Parameter)+
  geom_hline(yintercept = 0.95, col = "red", linetype = "dotted")
