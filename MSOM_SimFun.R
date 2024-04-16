#######################################
#### MSOM Simulation funciton     #####
#######################################


# function to derive general parameters (psi) from natural parameters (f's)
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

# psi.fun(f = cbind(beta1, beta2, beta12))


# Function to simulate data and fit a 2-species 
# multi-species occupancy model (MSOM) from Rota et al. (2016)

## INPUT
# beta1: vector with regression coefficient(s) for linear predictor (in logit-scale)
#       for first-order natural parameter f1 (i.e. species 1 occurrence)
# beta2: vector with regression coefficient(s) for linear predictor (in logit-scale)
#       for first-order natural parameter f2 (i.e. species 1 occurrence)
# beta12: vector with regression coefficient(s) for linear predictor (in logit-scale)
#       for second-order natural parameter f12 (i.e. species 1 & 2 c0-occurrence)
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


MSOM_SimandFit.fun <- function(beta1, beta2, beta12, nsites, nspecies = 2, J = 3, nocc_covs = 3,
                               p_true= c(0.5, 0.5) , ndet_covs = 2,
                               occ_formulas = rep("~1", 2^2-1),
                               det_formulas = rep("~1", 2), nsim =100, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  init.time <- Sys.time()
  require(unmarked, quietly = TRUE)
  require(tidyverse, quietly = TRUE)
  
  for(N in nsites){
    occ_covs <- as.data.frame(matrix(rnorm(N * nocc_covs),ncol=nocc_covs))
    names(occ_covs) <- paste('occ_cov',1:nocc_covs,sep='')
    
    det_covs <- list()
    for (i in 1:ndet_covs){
      det_covs[[i]] <- matrix(rnorm(N*J),nrow=N)
    }
    
    names(det_covs) <- paste('det_cov',1:ndet_covs,sep='')
    
    
    if(length(det_formulas) > nspecies) stop("Error: Number of linear predictors for detection should be equal to number of species") 
    ## True parameter values
    ### Single species occupancies
    n.nat.params <- 2^nspecies-1
    
    occ_dms <- lapply(X = as.list(occ_formulas), 
                      function(x)model.matrix(as.formula(x), data = occ_covs))
    
    # CHANGE THIS TO BE ADAPTIVE FOR 2+ SPP :
    if(ncol(occ_dms[[1]]) != length(beta1)){
      stop("Error: Number of regression coefficients for linear predictor f1 must equal number of covariates in formula + 1")
      }
      f1 <- occ_dms[[1]]%*%beta1
      
    if(ncol(occ_dms[[2]]) != length(beta2)){
      stop("Error: Number of regression coefficients for linear predictor f2 must equal number of covariates in formula + 1")
      }
      f2 <-occ_dms[[2]]%*%beta2
      if(length(occ_formulas) < n.nat.params ) {
        f12 <- 0
      } else {
        if(ncol(occ_dms[[3]]) != length(beta12)){
          stop("Error: Number of regression coefficients for linear predictor f12 must equal number of covariates in formula + 1")
        }
        f12 <-occ_dms[[3]]%*%beta12
      }
    
    f <- cbind(f1,f2,f12)
    z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
    colnames(z) <- paste('sp',1:nspecies,sep='')
    dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)
    
    psi <- exp(f %*% t(dm))
    psi <- psi/rowSums(psi)
    
    # psi <- psi.fun(f = f, nspecies = nspecies) # not worth within loop
    
    
    for(s in 1:nsim){
      #True state
      ztruth <- matrix(NA,nrow=N,ncol=nspecies)
      for (i in 1:N){
        if(nrow(psi) ==1 ){ 
          ztruth[i,] <- as.matrix(z[sample(4,1,prob=psi),])
        } else{     
          ztruth[i,] <- as.matrix(z[sample(4,1,prob=psi[i,]),])
        }
      }
      
      
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
      names(y) <- c('sp1','sp2')
      data <- unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)
      
      # Without penalty
      fit <- occuMulti(det_formulas,occ_formulas,data)
      fit_sum <-summary(fit)
      mod.conv <- ifelse(fit@opt$convergence==0, 1, 0)
      
      # with penalty
      fit_pl <- optimizePenalty(fit)
      pen.val <- as.numeric(as.character(fit_pl@call['penalty']))
      fit_pl.sum <- summary(fit_pl)
      mod.conv_pl <- ifelse(fit_pl@opt$convergence==0, 1, 0)
      
      
      if(s == 1){
        #  Without penalisation
        ## Detection parameters and data frames
        det_df <- fit_sum$det
        det_df$Species <- names(y)
        det_df$niter <- rep(s, nrow(det_df))
        det_df$conv <- mod.conv
        det_df$og_param.val <- qlogis(p_true)
        
        ## State (occ) parameters and data frames
        state_df <- fit_sum$state
        state_df$Parameter <- c(paste0("beta1.", seq(0, length(beta1)-1)),
                                paste0("beta2.", seq(0, length(beta2)-1)),
                                paste0("beta12.", seq(0, length(beta12)-1)))
        state_df$niter <- rep(s, nrow(state_df))
        state_df$conv <- mod.conv
        state_df$og_param.val <- c(beta1, beta2, beta12)
        
        #  With penalisation
        ## Detection parameters and data frames
        det_df.pl <- fit_pl.sum$det
        det_df.pl$Species <- names(y)
        det_df.pl$niter <- rep(s, nrow(det_df.pl))
        det_df.pl$conv <- mod.conv_pl
        det_df.pl$og_param.val <- qlogis(p_true)
        det_df.pl$lambda <- pen.val
        
        ## State (occ) parameters and data frames
        state_df.pl <- fit_pl.sum$state
        state_df.pl$Parameter <- c(paste0("beta1.", seq(0, length(beta1)-1)),
                                   paste0("beta2.", seq(0, length(beta2)-1)),
                                   paste0("beta12.", seq(0, length(beta12)-1)))
        state_df.pl$niter <- rep(s, nrow(state_df.pl))
        state_df.pl$conv <- mod.conv_pl
        state_df.pl$og_param.val <- c(beta1, beta2, beta12)
        state_df.pl$lambda <- pen.val
        
      } else {
        #  Without penalisation
        ## Detection parameters and data frames
        det_df1 <- fit_sum$det
        det_df1$Species <- names(y)
        det_df1$niter <- rep(s, nrow(det_df1))
        det_df1$og_param.val <- qlogis(p_true)
        det_df1$conv <- mod.conv
       
        
        ## State (occ) parameters and data frames
        state_df1 <- fit_sum$state
        state_df1$Parameter <- c(paste0("beta1.", seq(0, length(beta1)-1)),
                                 paste0("beta2.", seq(0, length(beta2)-1)),
                                 paste0("beta12.", seq(0, length(beta12)-1)))
        state_df1$niter <- rep(s, nrow(state_df1))
        state_df1$conv <- mod.conv
        state_df1$og_param.val <- c(beta1, beta2, beta12)
       
        
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
        
        ## State (occ) parameters and data frames
        state_df.pl1 <- fit_pl.sum$state
        state_df.pl1$Parameter <- c(paste0("beta1.", seq(0, length(beta1)-1)),
                                    paste0("beta2.", seq(0, length(beta2)-1)),
                                    paste0("beta12.", seq(0, length(beta12)-1)))
        state_df.pl1$niter <- rep(s, nrow(state_df.pl1))
        state_df.pl1$conv <- mod.conv_pl
        state_df.pl1$og_param.val <- c(beta1, beta2, beta12)
        state_df.pl1$lambda <- pen.val
        
        det_df.pl <- rbind(det_df.pl, det_df.pl1)
        state_df.pl <- rbind(state_df.pl, state_df.pl1)
        
      }
    } #s
    det_df$n.sites <- N
    state_df$n.sites <- N
    
    det_df.pl$n.sites <- N
    state_df.pl$n.sites <- N
    if(N == nsites[1]){
      # Without penalisation
      N_det_df <- det_df
      N_state_df <- state_df
      
      # With penalisation
      N_det_df.pl <- det_df.pl
      N_state_df.pl <- state_df.pl
      
    } else{
      N_det_df <- rbind(N_det_df, det_df)
      N_state_df <- rbind(N_state_df, state_df)
      
      #With penalisation
      N_det_df.pl <- rbind(N_det_df.pl, det_df.pl)
      N_state_df.pl <- rbind(N_state_df.pl, state_df.pl)
    }
  }# N?
  
  diff_time <- Sys.time()-init.time
  
  return(out = list(Detection.Params = N_det_df, State.params = N_state_df, 
                    Detection.Params.PL = N_det_df.pl, State.params.pl =N_state_df.pl,
                    time.ellapsed = diff_time))
}





