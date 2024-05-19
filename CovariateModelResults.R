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

source("MSOM_SimFun.R")

# Scenario 1 ----

## Read object
scen.cov1 <- read_rds("Scenario1_2spCovariateseed1337_simResults_v2.rds")


## Normal likelihood 
State_simcov1 <- scen.cov1$State.params %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()



### Summary of Normal likelihood (LL)
State_simcov1.sum <- State_simcov1 %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias), mu.CV = mean(CV, na.rm = T),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered, na.rm = T),                               # Coverage rate
            PWR = mean(below.alpha, na.rm = T),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            Lik = "LL") %>% 
  ungroup()

ggplot(State_simcov1.sum) +
  geom_point(aes(x= n.sites, y = Conv.rate))+
  facet_wrap(~Parameter)



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim1 <- State_simcov1 %>%
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



## Penalised likelihood 

State_simcov1.pl <- scen1$State.params.pl %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_simcov1.pl.sum <- State_simcov1.pl %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            Lik = "PL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim1.pl <- State_simcov1.pl %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(psi.fun(f = cbind(t(Estimate))))),
            og.psi = as.vector(t(psi.fun(f = cbind(t(og_param.val))))),
            Gen.Par = colnames(psi.fun(f = cbind(t(Estimate))))) %>%
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
            Lik = "PL") %>% 
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
            Lik = "PL") %>% 
  ungroup()



Cond.prob1.pl.sum <- Cond.prob1.pl %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            Lik = "PL") %>% 
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

## Read object
scen2 <- read_rds("scen2_seed1337_simResults2.rds")




## Normal likelihood 
State_sim2 <- scen2$State.params %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim2.sum <- State_sim2 %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            Lik = "LL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim2 <- State_sim2 %>%
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
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim2.pl.sum <- State_sim2.pl %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            Lik = "PL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim2.pl <- State_sim2.pl %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(psi.fun(f = cbind(t(Estimate))))),
            og.psi = as.vector(t(psi.fun(f = cbind(t(og_param.val))))),
            Gen.Par = colnames(psi.fun(f = cbind(t(Estimate))))) %>%
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
            Lik = "PL") %>% 
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

## Read object
scen3 <- read_rds("scen3_seed1337_simResults2.rds")




## Normal likelihood 
State_sim3 <- scen3$State.params %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim3.sum <- State_sim3 %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            Lik = "LL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim3 <- State_sim3 %>%
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


## Penalised likelihood 

State_sim3.pl <- scen3$State.params.pl %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim3.pl.sum <- State_sim3.pl %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            Lik = "PL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim3.pl <- State_sim3.pl %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(psi.fun(f = cbind(t(Estimate))))),
            og.psi = as.vector(t(psi.fun(f = cbind(t(og_param.val))))),
            Gen.Par = colnames(psi.fun(f = cbind(t(Estimate))))) %>%
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
            Lik = "PL") %>% 
  ungroup()



Cond.prob3.pl.sum <- Cond.prob3.pl %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
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


## Read object
scen4 <- read_rds("scen4_seed1337_simResults2.rds")




## Normal likelihood 
State_sim4 <- scen4$State.params %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim4.sum <- State_sim4 %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            Lik = "LL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim4 <- State_sim4 %>%
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

## Penalised likelihood 

State_sim4.pl <- scen4$State.params.pl %>%
  group_by(Parameter) %>%
  mutate(bias = (Estimate-og_param.val), prop.bias = (Estimate-og_param.val)/og_param.val,
         #prop.bias = abs((Estimate-og_pam)/og_pam),
         CV = SE/abs(Estimate), 
         Covered = ifelse(og_param.val >= Estimate - 1.96*SE & og_param.val <=Estimate+1.96*SE, 1, 0),
         below.alpha = ifelse(`P(>|z|)` < 0.05, 1, 0)) %>%
  ungroup()

### Summary of Normal likelihood (LL)
State_sim4.pl.sum <- State_sim4.pl %>%
  group_by(n.sites, Parameter) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
            CR = mean(Covered),                               # Coverage rate
            PWR = mean(below.alpha),                          # Mean type I error rate
            Conv.rate = mean(conv),                           # Mean convergence rate
            Lik = "PL") %>% 
  ungroup()



#### General parameters (psi11, psi10, psi01, psi00)


GenParam_sim4.pl <- State_sim4.pl %>%
  group_by(n.sites, niter) %>%
  summarise(psi.est = as.vector(t(psi.fun(f = cbind(t(Estimate))))),
            og.psi = as.vector(t(psi.fun(f = cbind(t(og_param.val))))),
            Gen.Par = colnames(psi.fun(f = cbind(t(Estimate))))) %>%
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
            Lik = "PL") %>% 
  ungroup()



Cond.prob4.pl.sum <- Cond.prob4.pl %>%
  group_by(n.sites, Cond.prob) %>%
  summarise(mu.bias = mean(bias), sd.bias = sd(bias),            # Mean and St.D of bias
            mu.p.bias = mean(prop.bias),
            bias.lci = mu.bias - 1.96*sd.bias,                # Lower bias CI
            bias.uci = mu.bias + 1.96*sd.bias,                # Upper bias CI
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
## The colour-blind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
### Scenario 1 -----
#### Natural parameters

##### Prop.bias
g1.nat <- ggplot(full.scen1.nat, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 1 (Str +ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-2.5, 2.5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")

##### Power

g1.pwr <- ggplot(full.scen1.nat%>%
                   filter(Parameter == "beta12.0"), aes(x = log(n.sites), y = PWR, col = Lik))+
  geom_point(size = 2.5, pch = 18)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  #facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Power", title = "Scenario 1 (Str +ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(0, 1)+
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "red")

#### General parameters
g1.gen <- ggplot(full.scen1.gen, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5, pch = 15)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Gen.Par)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 1 (Str +ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


#### Derived parameters

##### Marginal probs
g1.mar <-  ggplot(full.scen1.mar, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Species)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 1 (Str +ve)")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")

##### conditional probs
g1.con <- ggplot(full.scen1.con, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5, pch = 15)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Cond.prob)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 1 (Str +ve)")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")



### Scenario 2 ------
#### Natural parameters

##### Prop.bias
g2.nat <- ggplot(full.scen2.nat, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 2 (Weak +ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-2.5, 2.5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


##### Power

g2.pwr <- ggplot(full.scen2.nat%>%
                   filter(Parameter == "beta12.0"), aes(x = log(n.sites), y = PWR, col = Lik))+
  geom_point(size = 2.5, pch = 18)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  #facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Power", title = "Scenario 2 (Weak +ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(0, 1)+
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "red")

#### General parameters
g2.gen <- ggplot(full.scen2.gen, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5, pch = 15)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Gen.Par)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 2 (Weak +ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


#### Derived parameters

##### Marginal probs
g2.mar <-  ggplot(full.scen2.mar, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Species)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 2 (Weak +ve)")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")

##### conditional probs
g2.con <- ggplot(full.scen2.con, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5, pch = 15)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Cond.prob)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 2 (Weak +ve)")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


### Scenario 3 -----
#### Natural parameters

##### Prop.bias
g3.nat <- ggplot(full.scen3.nat, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 3 (Str -ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-5.5, 5.5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


##### Power

g3.pwr <- ggplot(full.scen3.nat%>%
                   filter(Parameter == "beta12.0"), aes(x = log(n.sites), y = PWR, col = Lik))+
  geom_point(size = 2.5, pch = 18)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  #facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Power", title = "Scenario 3 (Str -ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(0, 1)+
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "red")

#### General parameters
g3.gen <- ggplot(full.scen3.gen, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5, pch = 15)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Gen.Par)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 3 (Str -ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-.8, .8)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


#### Derived parameters

##### Marginal probs
g3.mar <-  ggplot(full.scen3.mar, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Species)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 3 (Str -ve)")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")

##### conditional probs
g3.con <- ggplot(full.scen3.con, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5, pch = 15)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Cond.prob)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 3 (Str -ve)")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")



### Scenario 4 -----
#### Natural parameters

##### Prop.bias
g4.nat <- ggplot(full.scen4.nat, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 4 (Weak -ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-5.5, 5.5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


##### Power

g4.pwr <- ggplot(full.scen4.nat%>%
                   filter(Parameter == "beta12.0"), aes(x = log(n.sites), y = PWR, col = Lik))+
  geom_point(size = 2.5, pch = 18)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  #facet_wrap(~Parameter)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Power", title = "Scenario 4 (Weak -ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(0, 1)+
  geom_hline(yintercept = 0.95, linetype = "dashed", col = "red")

#### General parameters
g4.gen <- ggplot(full.scen4.gen, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5, pch = 15)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Gen.Par)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 4 (Weak -ve)")+
  scale_colour_manual(values= c(cbbPalette[1], cbbPalette[2]))+ylim(-.8, .8)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")


#### Derived parameters

##### Marginal probs
g4.mar <-  ggplot(full.scen4.mar, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Species)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 4 (Weak -ve)")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")

##### conditional probs
g4.con <- ggplot(full.scen4.con, aes(x = log(n.sites), y = mu.p.bias, col = Lik))+
  geom_point(size = 2.5, pch = 15)+
  #geom_errorbar(aes(ymin = bias.lci, ymax = bias.uci ))+
  facet_wrap(~Cond.prob)+
  theme_bw()+
  labs(x = "log(Number of Sites)", y = "Mean Prop bias", title = "Scenario 4 (Weak -ve)")+
  scale_colour_manual(values= c(cbbPalette[3], cbbPalette[7]))+ylim(-.5, .5)+
  geom_hline(yintercept = 0, linetype = "dashed", col = "green")



### Cowplot parade! ----

## Bias plots ----

### Natural parameters 
#### Positive scenarios
cowplot::plot_grid(g1.nat+theme(axis.title.x = element_blank(),
                                axis.text.x = element_blank()),
                   g2.nat, nrow = 2)

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



