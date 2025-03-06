#######################################
#### MSOM Simulation function     #####
#######################################

# Hello! In this script you will find three functions used in our study to simulate, 
# fit and extract the data and parameter values for the Rota et al. (2016) co-occurrence model



# psi.fun ----
# Description: function to derive general parameters (psi) from natural parameters (f's) 

## INPUT
# f: matrix with ncol = 2^nspecies-1 ordered by paramter order-species (e.g. with 2 species: f= cbind(f1,f2,f12)) 
#   and nrow = number of observations arising from linear predictor. If no covariates, nrow(f) can be 1
# nspecies: number of species (must be >=1)

## OUTPUT:
# psi: matrix with 2^nspecies and nrow = nrow(f) columns indicating the probabilities 
#  of each of the 2^nspecies possible states of occupancy

psi.fun <- function(f, nspecies = 2) {
  if(ncol(f)<nspecies){
    stop(print(paste0("f must contain at leat first-order parameters for the ", nspecies, " species" )))
  }
  
  if(ncol(f) < 2^nspecies-1) {
    f <- cbind(f, matrix(0, nrow = f, ncol = (2^nspecies-1 -ncol(f))))
  }
  
  if(ncol(f) > 2^nspecies-1) {
    stop(print(paste0("ncol(f) must be smaller than maximum number of natural parameters (2^nspecies-1): ",2^nspecies-1 )))
  }
  
  z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
  colnames(z) <- paste('sp',1:nspecies,sep='')
  dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)
  
  psi <- exp(f %*% t(dm))
  psi <- psi/rowSums(psi)
  
  colnames(psi) <- apply(z, 1, function(x) paste0(x, collapse = ""))
  return(psi )
}



# crv.Clipp21 -----
# Description: Function adapted from Clipp et al. (2021) to get the likelihood under the chosen
# penalty value through K-cross-validation (CV)


## INPUT
# l: penalty value for which to compute the cross validation
# data: an 'unmarkedFrameOccuMulti' object containing the dataset
# M: the number of partitions (i.e. folds) in which to split the data

## Output
# sum(cv): sum of the log-likelihood of all cross-validation folds

crv.Clipp21 <- function(l, data, M = 5){
  seed = seed
  # number of sites
  n <- nrow(data@y)
  
  # index of where to start and stop vector of site indices
  s <- seq(1, by = n / M, length.out = M)  # starting index
  e <- seq(n / M, by = n / M, length.out = M)  # ending index
  
  # initializing vector of cross validation values
  cv <- numeric(M)
  
  for(i in 1:M){  # looping through all partitions
    
    # fitting to withheld data
    # error handling so simulations don't stop with error
    err_tmp <- class(try(
      fit <- occuMulti(detformulas = det_formulas,
                       stateformulas = occ_formulas,
                       data = data[(1:n)[-(s[i]:e[i])]],
                       maxOrder = 2,
                       penalty = l, se = T, boot = 100))) == 'try-error'
    
    
    # likelihood of withheld data
    if(err_tmp){
      cv[i] <- NA
    } else{
      cv[i] <- sum(unmarked:::occuMultiLogLik(fit, data[s[i]:e[i]]))
    }
  }
  
  return(sum(cv))
  
}

# MSOM_simfit.fun.v2 ------
# Description: Function to simulate data and fit a 2-species 
# multi-species occupancy model (MSOM) from Rota et al. (2016)

## INPUT
# beta: list with length nspecies+nspecies*(nspecies-1)/2 conaining the vectors
#       of the coefficient regressions for each natural parameter. Pairwise
#       independence can be set by setting 2nd order beta to 0. Order of natural
#       parameters as follows: 
#       paste0("f", c(1:nspecies, apply(combn(1:nspecies, 2), 2,paste0, collapse = "" )))
# nsites: value or vector of the number of sites to simulate data for
# nspecies: number of species to simulate data for (currently restricted to 2 species)
# J: number of visits in the detection process (default = 3 as per Pollock's robust method)
# nocc_covs: number of occupancy covaraites to generate
# p_true: true detection probabilities (in the real scale) for each species (length = nspecies)
# ndet_covs: number of detection covariates to simulate (without any use as per 12/11/2023)
# occ_formulas: vector with the formula for each linear predictor for the natural parameters
#               following the formula notation of "~covariate". Null modell specified with "~1"
# det_formulas: vector with the formula for each linear predictor for the observation process
#               for each species following notation of "~covariate". Null modell specified with "~1".
#               Length must equal nspecies. 
# nsim: number of datasets to simulate and fit for each sample size scenario
# seed: seed to ensure replicability of value generatio and model fits

## Output:
## Object containting:
# Detection.Params: Data.frame containing the detection coefficient estimates, 
#                    st.error and p-value for the nsim*length(nsites) models fit
#                    for the detection parameters (logit-scale). Also contains indicator columns
#                    of model convergence (1 if converged else 0)
# Detection.Params.pl: Data.frame containing the detection coefficient estimates,
#                       st.error and p-value for the nsim*length(nsites) models fit
#                       for the detection parameters (logit-scale) and the shrinkage coefficient from the penalised likelihood.
#                       Also contains indicator columns of model convergence (1 if converged else 0)
# State.params: Data.frame containing the state coefficient estimates, st.error and p-value
#                 for the nsim*length(nsites) models fit for the occupancy parameters. Also 
#                  contains indicator columns of model convergence (1 if converged else 0)
# State.params.pl: Data.frame containing the state coefficient estimates, st.error and p-value
#                   for the nsim*length(nsites) models fit for the occupancy parameters. Also 
#                   contains indicator columns of model convergence (1 if converged else 0) and
#                   estimated shirnkage parameter value for model fit with penalised likelihood. 
# time.ellapsed: diff.time object that represents how long the scenario took to simulate

MSOM_simfit.fun.v2 <- function(beta, nsites, nspecies = 2, J = 3, nocc_covs = 3,
                               p_true= c(0.5, 0.5) ,
                               occ_formulas = NULL,
                               det_formulas = NULL, nsim =100,
                               seed = NULL, store.data = T, sim.only = F) {
  
  if(!store.data & sim.only) stop("Error: Cannot simulate data only and not store it")
  if(!is.list(beta) | length(beta) != (nspecies+nspecies*(nspecies-1)/2)){ 
    stop("Error: beta should be a list with length = nspecies+nspecies*(nspecies-1)/2")
  }
  
  ndet_covs = nspecies
  if(!is.null(seed)) set.seed(seed)
  init.time <- Sys.time()
  require(unmarked, quietly = TRUE)
  require(tidyverse, quietly = TRUE)
  
  if(is.null(occ_formulas))  occ_formulas = rep("~1", nspecies+ncol(combn(1:nspecies, 2)))
  if(is.null(det_formulas))  det_formulas = rep("~1", nspecies)
  #### Debugging settings
   # nsites = 50; nsim = 2; ndet_covs = nocc_covs = nspecies; sim.only = F; store.data = F  #only to use for debugging
   
   #######~~~~~~~~~~~~~~~~~~~~~
  for(N in nsites){
    occ_covs <- as.data.frame(matrix(rnorm(N * nocc_covs),ncol=nocc_covs))
    names(occ_covs) <- paste('occ_cov',1:nocc_covs,sep='')
    
    det_covs <- list()
    for (i in 1:ndet_covs){
      det_covs[[i]] <- matrix(rnorm(N*J),nrow=N)
    }
    
    names(det_covs) <- paste('det_cov',1:ndet_covs,sep='')
    
    
    if(length(det_formulas) != nspecies) stop("Error: Number of linear predictors for detection should be equal to number of species") 
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
    
    s = 1
    while(s <= nsim){
      #True state
      ztruth <- matrix(NA,nrow=N,ncol=nspecies)
      for (i in 1:N){
        if(nrow(psi) ==1 ){ 
          ztruth[i,] <- as.matrix(z[sample(2^nspecies,1,prob=psi),])
        } else{     
          ztruth[i,] <- as.matrix(z[sample(2^nspecies,1,prob=psi[i,]),])
        }
      }
      
      if(any(colSums(ztruth) == 0)) next
      
      
      # fake y data
      y <- list()
      
      for (i in 1:nspecies){
        y[[i]] <- matrix(NA,N,J)
        for (j in 1:N){
          for (k in 1:J){
            y[[i]][j,k] <- rbinom(1,1,ztruth[j,i]*p_true[i])
          }
        }
      }
      
      # Check if all species have been detected at least once
      if(any(lapply(y, sum) == 0)) next
      names(y) <- paste0('sp', 1:nspecies)
      data <- unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)
      
      if(!sim.only){
      
        cat(paste(paste0("NSites: ", N),
                    paste0("SIMULATION ITERATION: ", s), sep = '\n'))
        # Without penalty
        fit <-try(occuMulti(det_formulas,occ_formulas,data, maxOrder = 2))
        if(class(fit) %in% "try-error") next
        fit_sum <-summary(fit)
        
        mod.conv <- ifelse(fit@opt$convergence==0 & !any(abs(c(fit_sum$state$Estimate, 
                                                              fit_sum$det$Estimate)) > 9.21) , 1, 0)
                           
        # with penalty (if boundary estimates run CV code from Clipp 2021)
        fit_pl <- try(optimizePenalty(fit))
        
        if("try-error" %in% class(fit_pl)){
          # Same estimates as Clipp et al 2021
          lam <- c(0.01, 0.02, 0.04, 0.08, 0.16, 
                   0.32, 0.64, 1.28, 2.56, 5.12)
          
          fit.ll <- numeric(length(lam))
          
          for(l in 1:length(lam)){
            # Try fitting doing the cross-validation lambda selection
            fit.ll[l] <- tryCatch(crv.Clipp21(l = lam[l], data = data, M = 5), 
                                  error = function(e) rep(NA, 1))
            if("try-error" %in% class(fit.ll[l])) next
          }
          
          if("try-error" %in% class(fit.ll[l])) next
          
          if(!all(is.na(fit.ll))){
            pen.val <- lam[which(fit.ll == max(fit.ll, na.rm = T))]
          } else{next}
          fit_pl <- occuMulti(det_formulas,occ_formulas,data, 
                              maxOrder = 2, penalty = pen.val, boot = 100)
          
          pen.est <- "CrossV"
        } else{
          
          pen.val <- as.numeric(as.character(fit_pl@call['penalty']))
          if(pen.val %in% c(2^-4, 2^4)){
            lam <- matrix(c(seq(0.001, 0.05, length = 8),
                          seq(18, 50, length = 8)),
                          ncol = 8, byrow = T)[which(pen.val == c(2^-4, 2^4)),]
            fit.ll <- numeric(length(lam))

            for(l in 1:length(lam)){
              fit.ll[l] <- tryCatch(crv.Clipp21(l = lam[l], data = data, M = 5),
                                    error = function(e) rep(NA, 1))
              if("try-error" %in% class(fit.ll[l])) next
            }

            if("try-error" %in% class(fit.ll[l])) next

            if(!all(is.na(fit.ll))){
              pen.val <- lam[which(fit.ll == max(fit.ll, na.rm = T))]
            } else{next}
            fit_pl <- occuMulti(det_formulas,occ_formulas,data,
                                maxOrder = 2, penalty = pen.val, boot = 100)
          }
          pen.est <- "optPen"
        }
        fit_pl.sum <- summary(fit_pl)
        
        mod.conv_pl <- ifelse(fit_pl@opt$convergence==0 & 
                                !any(abs(c(fit_pl.sum$state$Estimate,
                                           fit_pl.sum$det$Estimate)) > 9.21),1, 0)
        
      } # sim.only?
      if(s == 1){
        if(store.data){
          y.dat <- data
          z.dat <- cbind(ztruth, nsim = s, n.sites = N)
        }
      if(sim.only) {
        s <- s + 1
      next
      }

        #  Without penalisation
        ## Detection parameters and data frames
        det_df <- fit_sum$det
        det_df$Species <- names(y)
        det_df$niter <- rep(s, nrow(det_df))
        det_df$conv <- mod.conv
        det_df$og_param.val <- qlogis(p_true)
        
        ## State (occ) parameters and data frames
        state_df <- fit_sum$state
        state_df$Parameter <- params
        state_df$niter <- rep(s, nrow(state_df))
        state_df$conv <- mod.conv
        state_df$og_param.val <- og.param.val
        
        #  With penalisation
        ## Detection parameters and data frames
        det_df.pl <- fit_pl.sum$det
        det_df.pl$Species <- names(y)
        det_df.pl$niter <- rep(s, nrow(det_df.pl))
        det_df.pl$conv <- mod.conv_pl
        det_df.pl$og_param.val <- qlogis(p_true)
        det_df.pl$lambda <- pen.val
        det_df.pl$pen.est <- pen.est
        
        ## State (occ) parameters and data frames
        state_df.pl <- fit_pl.sum$state
        state_df.pl$Parameter <- params
        state_df.pl$niter <- rep(s, nrow(state_df.pl))
        state_df.pl$conv <- mod.conv_pl
        state_df.pl$og_param.val <- og.param.val
        state_df.pl$lambda <- pen.val
        state_df.pl$pen.est <- pen.est
        
      } else {
        if(store.data){
          y.dat <-c(y.dat,  data)
          z.dat <- rbind(z.dat, cbind(ztruth, nsim = s, n.sites = N))
        }
        if(sim.only)  {
          s <- s + 1
          next
        }
        
        
        #  Without penalisation
        ## Detection parameters and data frames
        det_df1 <- fit_sum$det
        det_df1$Species <- names(y)
        det_df1$niter <- rep(s, nrow(det_df1))
        det_df1$og_param.val <- qlogis(p_true)
        det_df1$conv <- mod.conv
       
        
        ## State (occ) parameters and data frames
        state_df1 <- fit_sum$state
        state_df1$Parameter <- params
        state_df1$conv <- mod.conv
        state_df1$niter <- rep(s, nrow(state_df1))
        state_df1$og_param.val <- og.param.val
       
        
        det_df <- rbind(det_df, det_df1)
        state_df <- rbind(state_df, state_df1)
        
        #  With penalisation
        ## Detection parameters and data frames
        det_df.pl1 <- fit_pl.sum$det
        det_df.pl1$Species <- names(y)
        det_df.pl1$niter <- rep(s, nrow(det_df.pl1))
        det_df.pl1$conv <- mod.conv_pl
        det_df.pl1$og_param.val <- qlogis(p_true)
        det_df.pl1$lambda <- pen.val
        det_df.pl1$pen.est <- pen.est
        ## State (occ) parameters and data frames
        state_df.pl1 <- fit_pl.sum$state
        state_df.pl1$Parameter <-params
        state_df.pl1$niter <- rep(s, nrow(state_df.pl1))
        state_df.pl1$conv <- mod.conv_pl
        state_df.pl1$og_param.val <- og.param.val
        state_df.pl1$lambda <- pen.val
        state_df.pl1$pen.est <- pen.est
        
        det_df.pl <- rbind(det_df.pl, det_df.pl1)
        state_df.pl <- rbind(state_df.pl, state_df.pl1)
        
      }
      s <- s + 1
    } #s
    det_df$n.sites <- N
    state_df$n.sites <- N
    
    det_df.pl$n.sites <- N
    state_df.pl$n.sites <- N
    
    if(N == nsites[1]){
      if(store.data){
        y.dat.full <- y.dat
        z.dat.full <- z.dat
      }
      if(sim.only) next
      
      # Without penalisation
      N_det_df <- det_df
      N_state_df <- state_df
      
      # With penalisation
      N_det_df.pl <- det_df.pl
      N_state_df.pl <- state_df.pl
      
    } else{
      if(store.data){
        y.dat.full <- c(y.dat.full, y.dat)
        z.dat.full <- rbind(z.dat.full, z.dat)
      }
      if(sim.only) next
      
      N_det_df <- rbind(N_det_df, det_df)
      N_state_df <- rbind(N_state_df, state_df)
      
      #With penalisation
      N_det_df.pl <- rbind(N_det_df.pl, det_df.pl)
      N_state_df.pl <- rbind(N_state_df.pl, state_df.pl)
    }
  }# N?
  diff_time <- Sys.time()-init.time
  
  
  if(sim.only){ return(out = list(z.truth = z.dat.full, sim.dat = y.dat.full,
                                  time.ellapsed = diff_time))}
  
  if(store.data){
    out = list(Detection.Params = N_det_df, State.params = N_state_df, 
             Detection.Params.PL = N_det_df.pl, State.params.pl =N_state_df.pl,
             z.truth = z.dat.full, sim.dat = y.dat.full,
             time.ellapsed = diff_time)
  } else{
    out = list(Detection.Params = N_det_df, State.params = N_state_df, 
               Detection.Params.PL = N_det_df.pl, State.params.pl =N_state_df.pl,
               time.ellapsed = diff_time)
  }
  return(out = out)
}


##########################
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
