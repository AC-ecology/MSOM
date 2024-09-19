########################################
##### Nul model results processing #####
########################################

# Workflow:
## 1. Load result object
## 2. Calculate bias, coverage rate, CV and power (for 2nd order and covariate terms)
## 3. Calculate bias, coverage rate for general parameters
## 4. Summarise results per sample size scenario
## 5.Plot with ggplot the power of 2nd order parameters and bias of all parameters

require(readr)
require(tidyverse)

### Shared features:

ln.sites <- seq(3, 8, by = 0.5)
nsites <- as.integer(exp(ln.sites))
nsim <- 100
nspecies <- 2

source("MSOM_SimFun.R")

# Scenario 1 ----

## Read object
scen1 <- read_rds("Results/scen1_seed1337_simResults2.rds")

## Normal likelihood 
State_sim1 <- scen1$State.params %>%
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

### Summary of Normal likelihood (LL)
State_sim1.sum <- State_sim1 %>%
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
GenParam_sim1 <- State_sim1 %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(psi.fun(f = cbind(t(Estimate))))),
            og.psi = as.vector(t(psi.fun(f = cbind(t(og_param.val))))),
            Gen.Par = colnames(psi.fun(f = cbind(t(Estimate))))) %>%
  ungroup()


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
            Lik = "LL") %>% 
  ungroup()


#### Derived Parameters 

##### Marginal probability

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
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()



Cond.prob1.sum <- Cond.prob1 %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()

# Scenario 2 ----

## Read object
scen2 <- read_rds("Results/scen2_seed1337_simResults2.rds")

## Normal likelihood 
State_sim2 <- scen2$State.params %>%
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

### Summary of Normal likelihood (LL)
State_sim2.sum <- State_sim2 %>%
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
GenParam_sim2 <- State_sim2 %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(psi.fun(f = cbind(t(Estimate))))),
            og.psi = as.vector(t(psi.fun(f = cbind(t(og_param.val))))),
            Gen.Par = colnames(psi.fun(f = cbind(t(Estimate))))) %>%
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
            Lik = "LL") %>% 
  ungroup()


#### Derived Parameters 
##### Marginal probability
Marg.prob2 <-  GenParam_sim2 %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()



##### Conditional probability
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

# Scenario 3 ----

## Read object
scen3 <- read_rds("Results/scen3_seed1337_simResults2.rds")

## Normal likelihood 
State_sim3 <- scen3$State.params %>%
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
### Summary of Normal likelihood (LL)
State_sim3.sum <- State_sim3 %>%
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
GenParam_sim3 <- State_sim3 %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(psi.fun(f = cbind(t(Estimate))))),
            og.psi = as.vector(t(psi.fun(f = cbind(t(og_param.val))))),
            Gen.Par = colnames(psi.fun(f = cbind(t(Estimate))))) %>%
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
            Lik = "LL") %>% 
  ungroup()


#### Derived Parameters 
##### Marginal probability
Marg.prob3 <-  GenParam_sim3 %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()

##### Coditional probability
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
            Lik = "LL") %>% 
  ungroup()


Cond.prob3.sum <- Cond.prob3 %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()

# Scenario 4 ----
## Read object
scen4 <- read_rds("Results/scen4_seed1337_simResults2.rds")

## Normal likelihood 
State_sim4 <- scen4$State.params %>%
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

### Summary of Normal likelihood (LL)
State_sim4.sum <- State_sim4 %>%
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
GenParam_sim4 <- State_sim4 %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(psi.fun(f = cbind(t(Estimate))))),
            og.psi = as.vector(t(psi.fun(f = cbind(t(og_param.val))))),
            Gen.Par = colnames(psi.fun(f = cbind(t(Estimate))))) %>%
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
            Lik = "LL") %>% 
  ungroup()

#### Derived Parameters 
##### Marginal probability
Marg.prob4 <-  GenParam_sim4 %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()

##### Conditional probability
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
            Lik = "LL") %>% 
  ungroup()

Cond.prob4.sum <- Cond.prob4 %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()

# Scenario 5 ----
## Read object
scen5 <- read_rds("Results/Scenario5_MidPosNullseed1337_simResults_v2.rds")

## Normal likelihood 
State_sim5 <- scen5$State.params %>%
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

### Summary of Normal likelihood (LL)
State_sim5.sum <- State_sim5 %>%
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
GenParam_sim5 <- State_sim5 %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(psi.fun(f = cbind(t(Estimate))))),
            og.psi = as.vector(t(psi.fun(f = cbind(t(og_param.val))))),
            Gen.Par = colnames(psi.fun(f = cbind(t(Estimate))))) %>%
  ungroup()


GenParam_sim5 <-  GenParam_sim5 %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()


GenParam_sim5.sum <- GenParam_sim5 %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()
#### Derived Parameters 
##### Marginal probability
Marg.prob5 <-  GenParam_sim5 %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()

##### Conditional probability
Cond.prob5 <-  GenParam_sim5 %>%
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

## Summaries of derived parameters
##### Marginal probability
Marg.prob5.sum <- Marg.prob5 %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()


##### Conditional probability
Cond.prob5.sum <- Cond.prob5 %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()


# Scenario 6 ----
## Read object
scen6 <- read_rds("Results/Scenario6_MidNegNullseed1337_simResults_v2.rds")

## Normal likelihood 
State_sim6 <- scen6$State.params %>%
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

### Summary of Normal likelihood (LL)
State_sim6.sum <- State_sim6 %>%
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
GenParam_sim6 <- State_sim6 %>%
  filter(conv == 1) %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(psi.fun(f = cbind(t(Estimate))))),
            og.psi = as.vector(t(psi.fun(f = cbind(t(og_param.val))))),
            Gen.Par = colnames(psi.fun(f = cbind(t(Estimate))))) %>%
  ungroup()

GenParam_sim6 <-  GenParam_sim6 %>%
  group_by(Gen.Par) %>%
  mutate(bias = (psi.est-og.psi), 
         prop.bias = (psi.est-og.psi)/og.psi) %>%
  ungroup()

GenParam_sim6.sum <- GenParam_sim6 %>%
  group_by(n.sites, Gen.Par) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()

#### Derived Parameters 
##### Marginal probability
Marg.prob6 <-  GenParam_sim6 %>%
  group_by(n.sites, niter) %>%
  summarise(Species = c("Sp1", "Sp2"),
            Marg.est = c(sum(psi.est[Gen.Par %in% c("10", "11")]),
                         sum(psi.est[Gen.Par %in% c("01", "11")])),
            Marg.og = c(sum(og.psi[Gen.Par %in% c("10", "11")]),
                        sum(og.psi[Gen.Par %in% c("01", "11")])),
            bias = Marg.est - Marg.og, prop.bias = (Marg.est - Marg.og)/Marg.og) %>%
  ungroup()


##### Conditional probability
Cond.prob6 <-  GenParam_sim6 %>%
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
Marg.prob6.sum <- Marg.prob6 %>%
  group_by(n.sites, Species) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()

Cond.prob6.sum <- Cond.prob6 %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "LL") %>% 
  ungroup()


### Only likelihood (all scenarios) -----

##### Natural parameters
full.scen.nat <- rbind(State_sim1.sum, State_sim2.sum, State_sim5.sum,
                       State_sim3.sum, State_sim4.sum, State_sim6.sum)
full.scen.nat$Scenario <- rep(paste0("Scenario", c(1, 2, 5, 3, 4, 6)), each = nrow(State_sim3.sum))
full.scen.nat$IntStrength <- rep(c("Strong", "Weak", "Moderate" , "Strong", "Weak",  "Moderate"), each = nrow(State_sim3.sum))
full.scen.nat$IntType <- rep(c("Positive interaction ", "Negative interaction"), each = nrow(State_sim3.sum)*3)
full.scen.nat$NatPam <- rep(c("f1.0", "f12.0", "f2.0"),nrow( full.scen.nat)/3)
full.scen.nat$CoeffReg <- rep(c("alpha0", "beta0", "gamma0"), nrow(full.scen.nat)/3)
full.scen.nat$Parameter <- rep(c("f1.0", "f12.0", "f2.0"),nrow( full.scen.nat)/3)

### General Parameters

full.scen.gen <- rbind(GenParam_sim1.sum, GenParam_sim2.sum, GenParam_sim5.sum,
                       GenParam_sim3.sum, GenParam_sim4.sum, GenParam_sim6.sum)

full.scen.gen$Scenario <- rep(paste0("Scenario", c(1, 2, 5, 3, 4, 6)), each = nrow(GenParam_sim5.sum))
full.scen.gen$IntStrength <- rep(c("Strong", "Weak", "Moderate" , "Strong", "Weak",  "Moderate"), each = nrow(GenParam_sim5.sum))
full.scen.gen$IntType <- rep(c("Positive interaction ", "Negative interaction"), each = nrow(GenParam_sim5.sum)*3)
full.scen.gen$Parameter <- full.scen.gen$Gen.Par



##### Marginal Probabilities
full.scen.mar <- rbind(Marg.prob1.sum, Marg.prob2.sum, Marg.prob5.sum,
                       Marg.prob3.sum, Marg.prob4.sum, Marg.prob6.sum)
full.scen.mar$Scenario <- rep(paste0("Scenario", c(1, 2, 5, 3, 4, 6)), each = nrow(Marg.prob1.sum))
full.scen.mar$IntStrength <- rep(c("Strong", "Weak", "Moderate" , "Strong", "Weak",  "Moderate"), each = nrow(Marg.prob1.sum))
full.scen.mar$IntType <- rep(c("Positive interaction ", "Negative interaction"), each = nrow(Marg.prob1.sum)*3)
full.scen.mar$Prob <- full.scen.mar$Species
full.scen.mar$ProbType <- "Marginal"
##### Conditional Probabilities
full.scen.con <- rbind(Cond.prob1.sum, Cond.prob2.sum, Cond.prob5.sum,
                       Cond.prob3.sum, Cond.prob4.sum, Cond.prob6.sum)
full.scen.con$Scenario <- rep(paste0("Scenario", c(1, 2, 5, 3, 4, 6)), 
                              each = nrow(Cond.prob5.sum))
full.scen.con$IntStrength <- rep(c("Strong", "Weak", "Moderate" , "Strong", "Weak",  "Moderate"),
                                 each = nrow(Cond.prob5.sum))
full.scen.con$IntType <- rep(c("Positive interaction ", "Negative interaction"),
                             each = nrow(Cond.prob5.sum)*3)


full.scen.con$Prob <- full.scen.con$Cond.prob
full.scen.con$ProbType <- "Conditional"

#### Generally Natural (combo) 


shared_cols <- intersect(colnames(full.scen.gen), colnames(full.scen.nat))

full.scen.natgen <- rbind(full.scen.nat[, shared_cols],full.scen.gen[,shared_cols] )

full.scen.natgen$Param.Type <- c(rep("Natural", nrow(full.scen.nat)),
                                 rep("General", nrow(full.scen.gen)))

#### Marginally conditional (combo)
shared_cols1 <- intersect(colnames(full.scen.con), colnames(full.scen.mar))
full.scen.mar.con <-  rbind(full.scen.mar[, shared_cols1],full.scen.con[,shared_cols1])
### Making Bias table  ------
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
# scales::show_col(Col.pallete)
# scales::show_col(alpha("red", seq(0.9, 0.4, length.out = 5)))

### Gen.Nat Long Table 1  --------------
natgen.long <- full.scen.natgen %>% 
  select(n.sites, IntType, IntStrength,  Parameter, Param.Type, mu.p.bias) %>%
  mutate(mu.p.bias = round(mu.p.bias*100, 2),
         IntType = ifelse(grepl("Positive",IntType), "pos", "neg"),
         IntStrength = factor(substr(IntStrength, 1, 1), levels = c("W", "M", "S"))) %>%
  pivot_wider(names_from = c(Param.Type, Parameter), values_from = mu.p.bias) %>%
  arrange(IntType, IntStrength )

order.vec1 <- c(rep(0, 3), 0, 1, -1,0, 1, -1, 0)

natgen.long1 <- natgen.long[c(1:ncol(natgen.long)+order.vec1)]

## Positive scenario
natgen.long.pos <- natgen.long1 %>% filter(IntType == "pos")

# Store table to be worked on
kable_out <- kableExtra::kbl(natgen.long1, booktabs = T, 
                             col.names = c("N", "Direction", "Strength", rep(c("f1", "f2", "f12"), 1),
                                           rep(c("00", "10", "01", "11"), 1)),
                             digits = 1,align = "c", 
                             linesep = "", format = "html") %>%
  # kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 2:3,latex_hline = "major", row_group_label_position = "first") 

# Loop over columns for column_spec
for (col_num in 4:ncol(natgen.long1)) {
  
  kable_out <- kableExtra::column_spec(kable_out,
                                       col_num,
                                       background =  unlist(lapply(natgen.long1[[col_num]],function(x) table.color(x) )))
}

kable_out


###### Without moderate and only positive ----
natgen.long.sub <- natgen.long1 %>% 
  filter(IntStrength  != "M" & IntType == "neg") %>%
  select(-IntType)
kable_out <- kableExtra::kbl(natgen.long.sub, booktabs = T, 
                             col.names = c("N","Strength", rep(c("f1", "f2", "f12"), 1),
                                           rep(c("00", "10", "01", "11"), 1)),
                             digits = 1,align = "c", 
                             linesep = "", format = "html") %>%
  add_header_above(c(" ", "", "Natural" = 3, "General" = 4)) %>%
  # kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 2:3,latex_hline = "major", row_group_label_position = "first") 

# Loop over columns for column_spec
for (col_num in 3:ncol(natgen.long.sub)) {
  
  kable_out <- kableExtra::column_spec(kable_out,
                                       col_num,
                                       background =  unlist(lapply(natgen.long.sub[[col_num]],function(x) table.color(x) )))
}

kable_out


#### Mar-Con table! ------------

marcon.long <- full.scen.mar.con %>% 
  select(n.sites, IntType, IntStrength,  Prob, mu.p.bias) %>%
  mutate(mu.p.bias = round(mu.p.bias*100, 2),
         IntType = ifelse(grepl("Positive",IntType), "pos", "neg"),
         IntStrength = factor(substr(IntStrength, 1, 1), levels = c("W", "M", "S"))) %>%
  pivot_wider(names_from = c(Prob), values_from = mu.p.bias) %>%
  arrange(IntType, IntStrength )

mar.con.sub <- marcon.long %>%
  filter(IntStrength  != "M" & IntType == "neg") %>%
  select(-IntType)

### Table 

colnames(mar.con.sub)
kable_out <- kableExtra::kbl(mar.con.sub, booktabs = T, 
                             col.names = c("N","Strength", colnames(mar.con.sub)[-c(1:2)]),
                             digits = 1,align = "c", 
                             linesep = "", format = "latex") %>%
  add_header_above(c(" ", "", "Marginal" = 2, "Conditional" = 4)) %>%
  # kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 2:3,latex_hline = "major", row_group_label_position = "first") 

# Loop over columns for column_spec
for (col_num in 3:ncol(mar.con.sub)) {
  
  kable_out <- kableExtra::column_spec(kable_out,
                                       col_num,
                                       background =  unlist(lapply(mar.con.sub[[col_num]],function(x) table.color(x) )))
}

kable_out

#### Wide Table(s)  --------------

natgen.wide <- full.scen.natgen %>% 
  select(n.sites, IntStrength, IntType, Parameter, Param.Type, mu.p.bias) %>%
  mutate(mu.p.bias = round(mu.p.bias*100, 2),
         IntType = ifelse(grepl("Positive",IntType), "pos", "neg"),
         IntStrength = substr(IntStrength, 1, 1)) %>%
  pivot_wider(names_from = c(IntType,IntStrength,Param.Type, Parameter), values_from = mu.p.bias)

order.vec <- c(0, rep(c(0, 1, -1), 6), rep(c(0, 1, -1,0), 6))
natgen.wide1 <- natgen.wide[c(1:ncol(natgen.wide)+order.vec)]


natgen.wide <- full.scen.natgen %>% 
  select(n.sites, IntStrength, IntType, Parameter, Param.Type, mu.p.bias) %>%
  mutate(mu.p.bias = round(mu.p.bias*100, 2),
         IntType = ifelse(grepl("Positive",IntType), "pos", "neg"),
         IntStrength = substr(IntStrength, 1, 1)) %>%
  pivot_wider(names_from = c(IntType,IntStrength,Param.Type, Parameter), values_from = mu.p.bias)

order.vec <- c(0, rep(c(0, 1, -1), 6), rep(c(0, 1, -1,0), 6))
natgen.wide1 <- natgen.wide[c(1:ncol(natgen.wide)+order.vec)]

kbl(natgen.wide1, booktabs = T, 
    col.names = c("N", rep(c("f1", "f2", "f12"), 6),
                  rep(c("00", "10", "01", "11"), 6)),digits = 1,align = "l", 
 linesep = "") %>%
  kable_styling(latex_options = "striped") %>%
  add_header_above(c(" ", "Strong" = 3, "Weak" = 3, "Mid" = 3,
                     "Strong" = 3, "Weak" = 3, "Mid" = 3,
                     "Strong" = 4, "Weak" = 4, "Mid" = 4,
                     "Strong" = 4, "Weak" = 4, "Mid" = 4)) %>%
  add_header_above(c(" ", "Positive" = 9, "Negative" = 9,
                     "Positive" = 12, "Negative" = 12)) %>%
  add_header_above(c(" ","Natural Parameters" = 18, "General Parameters" = 24), bold = T, italic = T) %>%
  kable_styling(latex_options = c("striped", "scale_down")) %>%
  column_spec(1, bold = T) %>%
  kable_styling(latex_options = c("striped", "hold_position")) %>%
  kable_classic()




## Test
colorines <- unlist(lapply(ok[,2],function(x) table.color(x) ))


# Store table to be worked on
kable_out <- kableExtra::kbl(natgen.wide1, booktabs = T, 
                             col.names = c("N", rep(c("f1", "f2", "f12"), 6),
                                           rep(c("00", "10", "01", "11"), 6)),
                             digits = 1,align = "c", format = "latex") %>%
  add_header_above(c(" ", "Strong" = 3, "Weak" = 3, "Mid" = 3,
                     "Strong" = 3, "Weak" = 3, "Mid" = 3,
                     "Strong" = 4, "Weak" = 4, "Mid" = 4,
                     "Strong" = 4, "Weak" = 4, "Mid" = 4)) %>%
  add_header_above(c(" ", "Positive" = 9, "Negative" = 9,
                     "Positive" = 12, "Negative" = 12)) %>%
  add_header_above(c(" ","Natural Parameters" = 18, "General Parameters" = 24), bold = T, italic = T) %>%
  kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T)

# Loop over columns for column_spec
for (col_num in 2:ncol(natgen.wide1)) {
  
  kable_out <- kableExtra::column_spec(kable_out,
                                       col_num,
                                       background =  unlist(lapply(natgen.wide1[[col_num]],function(x) table.color(x) )))
}

kable_out

##### Table without moderate scenario ----

natgen.wide <- full.scen.natgen %>% 
  filter(IntStrength != "Moderate") %>%
  select(n.sites, IntStrength, IntType, Parameter, Param.Type, mu.p.bias) %>%
  mutate(mu.p.bias = round(mu.p.bias*100, 2),
         IntType = ifelse(grepl("Positive",IntType), "pos", "neg"),
         IntStrength = substr(IntStrength, 1, 1)) %>%
  pivot_wider(names_from = c(IntType,IntStrength,Param.Type, Parameter), values_from = mu.p.bias)

order.vec <- c(0, rep(c(0, 1, -1), 4), rep(c(0, 1, -1,0), 4))
natgen.wide1 <- natgen.wide[c(1:ncol(natgen.wide)+order.vec)]


kbl(natgen.wide1, booktabs = T, 
    col.names = c("N", rep(c("f1", "f2", "f12"), 4),
                  rep(c("00", "10", "01", "11"), 4)),digits = 1,align = "l", 
    linesep = "") %>%
  kable_styling(latex_options = "striped") %>%
  add_header_above(c(" ", "Strong" = 3, "Weak" = 3, 
                     "Strong" = 3, "Weak" = 3, 
                     "Strong" = 4, "Weak" = 4, 
                     "Strong" = 4, "Weak" = 4 )) %>%
  add_header_above(c(" ", "Positive" = 6, "Negative" = 6,
                     "Positive" = 8, "Negative" = 8)) %>%
  add_header_above(c(" ","Natural Parameters" = 12, "General Parameters" = 16), bold = T, italic = T) %>%
  kable_styling(latex_options = c("striped", "scale_down")) %>%
  column_spec(1, bold = T) %>%
  kable_styling(latex_options = c("striped", "hold_position")) %>%
  kable_classic()


colfunc<-colorRampPalette(c("red","yellow","#71C93F"))

Col.pallete <- alpha(colfunc(6), 0.5) # n = number of steps (sequentially change)

## Test
colorines <- unlist(lapply(ok[,2],function(x) table.color(x) ))

kbl(ok, booktabs = T, escape = F) %>% 
  column_spec(column = c(2), background =colorines )%>%
  kable_classic() %>%
  kable_styling(latex_options = c("scale_down")) 



# Store table to be worked on
kable_out <- kableExtra::kbl(natgen.wide1, booktabs = T, 
                             col.names = c("N", rep(c("f1", "f2", "f12"), 4),
                                           rep(c("00", "10", "01", "11"), 4)),digits = 1,align = "l", 
                             linesep = "", format = "latex") %>%
  add_header_above(c(" ", "Strong" = 3, "Weak" = 3, 
                     "Strong" = 3, "Weak" = 3, 
                     "Strong" = 4, "Weak" = 4, 
                     "Strong" = 4, "Weak" = 4 )) %>%
  add_header_above(c(" ", "Positive" = 6, "Negative" = 6,
                     "Positive" = 8, "Negative" = 8)) %>%
  add_header_above(c(" ","Natural Parameters" = 12, "General Parameters" = 16), bold = T, italic = T) %>%
  # kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T)

# Loop over columns for column_spec
for (col_num in 2:ncol(natgen.wide1)) {
  
  kable_out <- kableExtra::column_spec(kable_out,
                                       col_num,
                                       background =  unlist(lapply(natgen.wide1[[col_num]],function(x) table.color(x) )))
}

kable_out


#### Table without moderate scenario and changing order ----------------

natgen.wide2 <- full.scen.natgen %>% 
  filter(IntStrength != "Moderate") %>%
  arrange(Scenario) %>%
  select(n.sites, IntStrength, IntType, Parameter, Param.Type, mu.p.bias) %>%
  mutate(mu.p.bias = round(mu.p.bias*100, 2),
         IntType = ifelse(grepl("Positive",IntType), "pos", "neg"),
         IntStrength = substr(IntStrength, 1, 1)) %>%
  pivot_wider(names_from = c(IntType,Param.Type,Parameter, IntStrength), values_from = mu.p.bias)

order.vec <- c(0, rep(c(0, 1, -1, 0, 1, -1,0), 4))

natgen.wide3 <- natgen.wide2[c(1:ncol(natgen.wide2)+order.vec)]



colfunc<-colorRampPalette(c("red","yellow","#71C93F"))

Col.pallete <- alpha(colfunc(6), 0.5) # n = number of steps (sequentially change)

## Test
colorines <- unlist(lapply(ok[,2],function(x) table.color(x) ))

kbl(ok, booktabs = T, escape = F) %>% 
  column_spec(column = c(2), background =colorines )%>%
  kable_classic() %>%
  kable_styling(latex_options = c("scale_down")) 



# Store table to be worked on
kable_out1 <- kableExtra::kbl(natgen.wide3, booktabs = T, 
                             col.names = c("N", rep(c("f_1", "f_2", "f_12","00", "10", "01", "11"), 4)),digits = 1,align = "l", 
                             linesep = "", format = "html") %>%
  add_header_above(c(" ", "Natural" = 3, "General" = 4,
                     "Natural" = 3, "General" = 4, 
                     "Natural" = 3, "General" = 4, 
                     "Natural" = 3, "General" = 4)) %>%
  add_header_above(c(" ", "Strong" = 7, "Weak" = 7,
                     "Strong" = 7, "Weak" = 7)) %>%
  add_header_above(c(" ","Positive" = 14, "Negative" = 14), bold = T, italic = T) %>%
  kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T)

# Loop over columns for column_spec
for (col_num in 2:ncol(natgen.wide3)) {
  
  kable_out1 <- kableExtra::column_spec(kable_out1,
                                       col_num,
                                       background =  unlist(lapply(natgen.wide3[[col_num]],function(x) table.color(x) )))
}

### Subset so we have 2 tables (1 for positive 1 for negative)
# Store table to be worked on
kable_out.pos <- kableExtra::kbl(natgen.wide1[,1:13], booktabs = T, 
                              col.names = c("N", rep(c("f_1", "f_2", "f_12","00", "10", "01", "11"), 2)),
                              digits = 1,align = "l", 
                              linesep = "", format = "latex") %>%
  add_header_above(c(" ", "Natural" = 3, "General" = 4,
                     "Natural" = 3, "General" = 4, 
                     "Natural" = 3, "General" = 4, 
                     "Natural" = 3, "General" = 4)) %>%
  add_header_above(c(" ", "Strong" = 7, "Weak" = 7,
                     "Strong" = 7, "Weak" = 7)) %>%
  add_header_above(c(" ","Positive" = 14, "Negative" = 14), bold = T, italic = T) %>%
  kable_styling(latex_options = c( "scale_down")) %>%
  column_spec(1, bold = T)

# Loop over columns for column_spec
for (col_num in 2:ncol(natgen.wide1)) {
  
  kable_out1 <- kableExtra::column_spec(kable_out1,
                                        col_num,
                                        background =  unlist(lapply(natgen.wide3[[col_num]],function(x) table.color(x) )))
}

## The plot parade ----

# Colourblind Palettes  

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
### All Scenarios -----

p.full.scen.v1 <- ggplot(full.scen.nat, aes(x = n.sites, y = mu.p.bias,
                                         col = IntStrength))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "black")+
  geom_line(aes(col = IntStrength))+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  facet_grid(IntType~CoeffReg)+
  labs(x = "Number of Sites", y = "Mean Relative Bias (RB)", col = "Interaction strength")+
  scale_colour_manual(values= c(cbbPalette[2], cbbPalette[3]))+
  coord_cartesian(ylim =  c(-2.5, 2.5))+
  theme_bw()

p.full.scen.v2 <- ggplot(full.scen.nat, aes(x = factor(n.sites), y = mu.p.bias,
                                            col = IntStrength, group = IntStrength))+
  geom_hline(yintercept = 0, linetype = "dashed", col = "black")+
  geom_line()+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  facet_grid(IntType~factor(NatPam, levels = c("f1", "f2", "f12")))+
  labs(x = "Number of Sites", y = "Mean Relative Bias (RB)", col = "Interaction strength")+
  scale_colour_manual(values= c(cbbPalette[2], cbbPalette[3]))+
  coord_cartesian(ylim =  c(-2.5, 2.5))+
  theme_bw()+guides(alpha = "none")+
  theme(axis.text.x = element_text(size  = 8))
p.full.scen.v2 + guides(alpha = "none")


### Power  -----
 ggplot(full.scen.nat %>%
                           filter(NatPam == "f12.0"), aes(x = factor(n.sites), y = PWR,
                                            col = IntStrength, group = IntStrength))+
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "red", lwd = 0.2)+
  geom_line(lwd = 0.5)+
  geom_point(size = 0.8)+
  scale_alpha_manual(values = c(1, 0.3))+
  facet_grid(~IntType)+
  labs(x = "Number of Sites", y = "Power", col = "Interaction strength")+
  scale_y_continuous(labels = scales::percent_format() )+
  scale_colour_manual(values= c(cbbPalette[4], cbbPalette[2], cbbPalette[3]))+
  theme_bw()+ guides(color = guide_legend(override.aes = list(size = 0.3)),
                     lwd =guide_legend(override.aes = list(lwd = 0.3)) )+
  theme(        strip.text.y = element_text(size = 5, color = "black", face = "bold"),
                strip.text.x = element_text(size = 5, color = "black", face = "bold"),
                axis.text.y = element_text(size = 4, color = "black"),
                axis.text.x =element_text(size = 3.5, color = "black") ,
                axis.title.x = element_text(size = 5, color = "black") ,
                axis.title.y = element_text(size = 5, color = "black"),
                legend.text = element_text(size = 3, margin = margin(l = 2)),
                legend.title = element_text(size = 4, margin = margin(b = 2) ),
                legend.box.spacing = unit(1, "mm"),
                legend.key.spacing.y = unit(0, "mm"),
                legend.key.size = unit(2,"mm"),
                legend.spacing.x = unit(0.1, "mm"))
 
 # ggsave("Figures/Null_AllScenarios_PowerNatPam_ModInt_v1.jpeg", 
 #        width = unit(4, "inches"),height = unit(2, "inches"), dpi = 600)


### Marginal probabilities ----

ggplot(full.scen.mar, 
       aes(x = log(n.sites), y = mu.p.bias,
                          col = Species, group = interaction(IntStrength, Species)))+
  geom_rect(aes(xmin = 2.5, xmax = 8.5, ymin = -0.05, ymax = 0.05), fill = "grey95",
            inherit.aes = F, alpha = 0.1)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "black")+
  geom_line()+
  geom_point(size = 2.5,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  facet_grid(IntType~factor(IntStrength, levels = c("Weak", "Moderate", "Strong")))+
  labs(x = "Number of Sites", y = "Mean Relative Bias (RB)", col = "Species")+
  scale_colour_manual(values= c(safe_colorblind_palette))+
  coord_cartesian(ylim =  c(-0.2, 0.2),  xlim = c(3, 8))+
  theme_bw()+guides(alpha = "none")+
  scale_x_continuous(breaks = log(nsites), labels = nsites)+
  scale_y_continuous(labels = scales::percent_format() )+
  theme(axis.text.x = element_text(size  = 8))


#### Marginal and general

ggplot(full.scen.mar, 
       aes(x = log(n.sites), y = mu.p.bias,
           col = Species, group = interaction(IntStrength, Species)))+
  geom_rect(aes(xmin = 2.5, xmax = 8.5, ymin = -0.05, ymax = 0.05), fill = "grey90",
            inherit.aes = F, alpha = 0.7)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "black")+
  geom_line(data = full.scen.gen,
            aes(x = log(n.sites), y =mu.p.bias, group = Gen.Par), col = "grey65", lwd = 0.2)+
  geom_line(lwd = 0.3)+
  geom_point(size = 0.4,
             aes(alpha = ifelse(abs(mu.p.bias)>0.05, "Biased", "Unbiased")))+
  scale_alpha_manual(values = c(1, 0.3))+labs(alpha = "RB>0.05")+
  facet_grid(IntType~factor(IntStrength, levels = c("Weak", "Moderate", "Strong")))+
  labs(x = "Number of Sites", y = "Mean Relative Bias (RB)", col = "Species")+
  scale_colour_manual(values= c(safe_colorblind_palette[-1]))+
  coord_cartesian(ylim =  c(-0.2, 0.2),  xlim = c(3, 8))+
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

# ggsave("Figures/Null_2Sp_AllScenarios_MarGen_ModInt.jpeg", 
#        width = unit(4, "inches"),height = unit(2, "inches"), dpi = 600)
# 






