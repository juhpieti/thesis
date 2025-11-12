####################################################################################
### THIS SCRIPT INCLUDES A LOTS OF FUNCTIONS THAT ARE UTILIZED IN OTHERS SCRIPTS ###
### FUNCTIONS ARE LOADED WITH load(path/helpers.R) IN OTHER SCRIPTS              ###
####################################################################################

# load in libraries
library(rdist) #for calculating distance matrices
library(cubature) #for calculating integrals

exp_covariance <- function(s1,s2,sigma2,l) {
  ### Calculates Covariance matrix using exponential covariance function
  # x1 is Nx2 matrix of observations
  # x2 is Mx2 matrix of observations
  # sigma2 = magnitude of the cov. function
  # l = lengthscale parameter of the cov. function
  # RETURNS NxN covariance matrix
  
  X.dist <- cdist(s1,s2) # pairwise distances 
  return(sigma2*exp(-X.dist/l))
}

scale_covariates <- function(X,X_new = NULL) {
  ### scales the covariates to zero mean unit variance
  ### mean(x) and sd(x) are calculated from the training data X
  ### if new covariates X_new are introduced, they are scaled with the parameters
  ### learned from the training data
  # X: matrix for covariates at sampling locations
  # X_new: matrix for covariates at new locations
  # RETURNS: scaled X if X_new not given (=NULL), scaled X_new otherwise
  
  # calculate column mean and standard deviaton
  means <- apply(X,2,mean)
  sds <- apply(X,2,sd)
  
  if (!is.null(X_new)) { # if X_new is given, return that instead of X
    X <- X_new
  }
  
  # scale to 0 mean 1 variance
  for(i in 1:ncol(X)) {
    X[,i] <- (X[,i] - means[i]) / sds[i]
  }
  return(X)
}

add_second_order_terms <- function(X,name_list) {
  ### adds second order terms to data matrix X
  # X: a data matrix (nxp)
  # name_list: vector of column names to add second order terms
  # RETURNS a new data matrix (nx(p+m)) where m is the number of second order terms
  for (name in name_list) {
    col_name <- paste0(name,"^2")
    X[,col_name] <- X[,name]^2
  }
  return(X)
}

add_interactions <- function(X,name1,name2) {
  ### adds interaction term between variables name1 and name2 in the design matrix.
  # X: a data matrix (nxp)
  # name1: name for the first variable
  # name2: name for the second variable
  col_name <- paste0(name1,"*",name2)
  X[,col_name] <- X[,name1]*X[,name2]
  return(X)
}

calc_loo <- function(stan_fit) {
  ### calculates leave-one-out log-score 
  ### requires that the stan_fit object includes output called log_lik (log-likelihoods are calculated in generate_quantities block)
  # stan_fit: stan fit object
  
  fit.loo <- loo(stan_fit)
  return(fit.loo$estimates[1,1]) # returns the sum of pointwise log-scores
}

inv_logit <- function(lin.pred) {
  ### inverse logit function 1/1+exp(-x)
  # lin.pred: linear predictor (x)
  return(1/(1+exp(-lin.pred)))
}


### CALCULATING EXPECTATIONS FOR Y WHEN Y ~ LC-BETA(MU,RHO)

integrate_LC_beta <- function(mu,rho,a,b=0.5,right_censored=FALSE) {
  ### calculates the expectation of y when y~LC-Beta(mu,rho)
  ### formulas from Pietilä thesis 2025
  # mu: mean of latent beta
  # rho: precision/sample size of latent beta
  # a: left-censoring constant
  # b: right-censoring constant
  # right_censored: Boolean to tell whether right-censoring is used
  
  ### integral crashes in certain cases, try to fix that beforehand
  ### this is probably not a good way to approach, think how to avoid this better?
  rho = max(rho,0.1) # integration crashes with very small rho ("very" U-shaped distribution)
  if (mu >= 0.9975) { # integration crashes if lots of mass near 1
    return(0.9975) # just return a high value
  }
  if ((mu > 0.99) & (rho < 0.25)) { # again large mu combained with small rho fails integral (explodes in the right end)
    return(0.995) # again, return a high value
  }
  
  # safe_integrate <- function(f,lower,upper) {
  #   tryCatch(integrate(f,lower,upper)$value,
  #            error = function (e) {
  #              message("Integration failed with mu=",mu,", rho=",rho, "| error: ", conditionMessage(e))
  #       return(NA)
  #     }
  #   )
  # }
  
  ### calculate the expectation
  if (!right_censored) {
    # calculate the integral by transforming the density of latent beta (again, see Pietilä Thesis 2025 for details)
    calc_density <- function(x) (x*dbeta((x+a)/(a+1),mu*rho,(1-mu)*rho)/(a+1)) 
    # integrate over (0,1) interval
    return(integrate(calc_density,0,1)$value)
    #return(safe_integrate(calc_density,0,1))
  } else {
    # if right-censoring is used, one needs to add 1xP(y=1), which is achieved when W > 1 <=> V > (1+a)/(1+a+b)
    calc_density <- function(x) (x*dbeta((x+a)/(a+b+1),mu*rho,(1-mu)*rho)/(a+b+1))
    return(1-pbeta((a+1)/(a+b+1),mu*rho,(1-mu)*rho)+integrate(calc_density,0,1)$value)
    #return(safe_integrate(calc_density,0,1))
  }
}

integrate_ZI_LC_beta <- function(mu,pi,rho,a,b=0.5,right_censored=FALSE) {
  ### calculates the expectation of y when y~LC-Beta(mu,rho) with zero-inflation
  ### formulas from Pietilä thesis 2025
  # mu: mean of latent beta
  # pi: probability of suitability
  # rho: precision/sample size of latent beta
  # a: left-censoring constant
  # b: right-censoring constant
  # right_censored: Boolean to tell whether right-censoring is used
  
  ### integral crashes in certain cases, try to fix that beforehand
  ### this is probably not a good way to approach, think how to avoid this better?
  rho = max(rho,0.1) # integration crashes with very small rho
  if (mu >= 0.995) { # integration crashes if lots of mass near 1
    return(pi*0.995) # just return a high value
  }
  if  ((mu > 0.99) & (rho < 0.25)) { #integration crashes if mean is high and rho is small (explodes in the right end)
    return(pi*0.995) # return a high value
  }
  
  ### calculate the expectation
  if (!right_censored) {
    # calculate the integral by transforming the density of latent beta (again, see Pietilä Thesis 2025 for details)
    calc_density <- function(x) (x*pi*dbeta((x+a)/(a+1),mu*rho,(1-mu)*rho)/(a+1))
    return(integrate(calc_density,0,1)$value)
  } else {
    # if right-censoring is used, one needs to add 1xP(y=1), which is achieved when W > 1 <=> V > (1+a)/(1+a+b)
    calc_density <- function(x) (x*pi*dbeta((x+a)/(a+b+1),mu*rho,(1-mu)*rho)/(a+b+1))
    return(pi*(1-pbeta((a+1)/(a+b+1),mu*rho,(1-mu)*rho)) + integrate(calc_density,0,1)$value)
  }
}

calculate_probzero_LC_beta <- function(mu,rho,a,b=0.5,right_censored=FALSE) {
  ### calculates the probability of zero when y~LC-Beta(mu,rho)
  ### this is just P(y=0) = P(W<0) = P(V<a/(a+1)), check the details from Pietilä thesis (2025)
  # mu: mean of latent beta
  # rho: precision/sample size of latent beta
  # a: left-censoring constant
  # b: right-censoring constant
  # right_censored: Boolean to tell whether right-censoring is used
  
  if (!right_censored) {
    return(pbeta(a/(a+1),mu*rho,(1-mu)*rho)) 
  } else {
    return(pbeta(a/(a+b+1),mu*rho,(1-mu)*rho))
  }
}

calculate_probzero_ZI_LC_beta <- function(mu,pi,rho,a,b=0.5,right_censored=FALSE) {
  ### calculates the probability of zero when y~LC-Beta(mu,rho) with zero-inflation
  ### this is just (1-pi) + pi*P(y=0) = 1-pi + pi*(V<a/(a+1)), check the details from Pietilä thesis (2025)
  # mu: mean of latent beta
  # rho: precision/sample size of latent beta
  # a: left-censoring constant
  # b: right-censoring constant
  # right_censored: Boolean to tell whether right-censoring is used
  if (!right_censored) {
    return(1-pi + pi*pbeta(a/(a+1),mu*rho,(1-mu)*rho)) 
  } else {
    return(1-pi + pi*pbeta(a/(a+b+1),mu*rho,(1-mu)*rho))
  }
}

### PREDICTING WITH LC-BETA MODELS

predict_beta_regression <- function(stan_fit, X.pred, X.orig, thinning = 10, a = 1, rho_modeled = FALSE, C = 1000, b = 0.5, right_censored = FALSE, min_rho=0) {
  ### function to make predictions with left-censored beta regression
  # X.pred: prediction matrix (mxp)
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # thinning: use every thinning:th posterior sample in predictions
  # a: left-censoring constant
  # rho_modeled: TRUE: rho modeled with covariates, FALSE: common rho
  # C: if rho modeled, what is it's upper limit? rho = C*inv_logit(a+XB)
  # b: right-censoring constant
  # right_censored: Boolean to tell whether right-censoring is used
  # min_rho: minimum possible value for rho(x) (if rho is modeled with covariates)
  # RETURNS list of rep x m matrix of samples from posterior predictive for different quantities
  ### 1) latent linear predictors f 
  ### 2) predicted Ys
  ### 3) expected Ys
  ### 4) probabilies of zero
  ### 5) rhos
  
  # load the posterior samples
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta_1","beta_2"))
  
  # prepare rho
  if (rho_modeled) {
    # coefficient related to rho(x)
    alpha_rho_sam <- as.matrix(stan_fit, pars = c("alpha_rho"))
    beta_rho_sam <- as.matrix(stan_fit, pars = c("beta_rho_1","beta_rho_2"))
  } else {
    # posterior sample of common rho
    rho_sam <- as.matrix(stan_fit, pars = c("rho"))
  }
  
  # prepare prediction matrix (scale, add second order terms)
  X.pred.scaled <- scale_covariates(X.orig,X.pred)
  Xpred <- add_second_order_terms(X.pred.scaled, colnames(X.pred))
  Xpred <- as.matrix(Xpred) # df to matrix for matrix calculations
  
  # initialize matrices for output
  f_sam <- c()
  y_sam <- c()
  EY_sam <- c()
  probzero_sam <- c()
  rho_sample <- c()
  
  # loop over posterior samples (take every thinning:th sample)
  for (i in seq(1,nrow(beta_sam),thinning)) {
    
    # corresponding coefficients
    alpha_i <- alpha_sam[i,]
    beta_i <- beta_sam[i,]
    
    # prepare rho
    if (rho_modeled){
      # rho(x) modeled with covariates
      rho_i <- min_rho + (C-min_rho)*inv_logit(as.vector(alpha_rho_sam[i,]  + Xpred %*% beta_rho_sam[i,]))
    } else {
      # every location has common rho
      rho_i <- rep(rho_sam[i,], nrow(Xpred))
    }
    
    # save the rhos for output
    rho_sample <- rbind(rho_sample,rho_i)

    # latent f (linear predictor a + xb)
    f_i <- as.vector(alpha_i + Xpred %*% beta_i)
    f_sam <- rbind(f_sam,f_i)
    
    # mean for beta distribution using logit-link
    mu_i <- inv_logit(f_i)
    
    ### calculate expectations and probability of zeroes
    ### NOTE: mapply iterates over pairs of (mu,rho) with common a
    EYi <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a))
    probzero_i <- mapply(calculate_probzero_LC_beta, mu=mu_i, rho=rho_i, MoreArgs = list(a=a))

    # save the expectations and prob. of zeros
    EY_sam <- rbind(EY_sam,EYi)
    probzero_sam <- rbind(probzero_sam, probzero_i)
    
    ### also predict Ys
    # sample latent Vs
    V_i <- rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
    
    # scale & censor
    if (!right_censored) {
      # scale to (-a,1), everything < 0 are treated as 0
      y_i <- sapply(V_i, function(v) (max(0,(a+1)*v - a)))
    } else {
      # if right-censoring is used, scale to (-a,1+b), everything < 0 are treated as 0, everything > 1 are treated as 1
      y_i <- sapply(V_i, function(v) (max(0, min(1,(a+b+1)*v - a))))
    }
    # save the sample
    y_sam <- rbind(y_sam, y_i)
  }
  
  # return everything
  return(list(f_sam = f_sam,
              y_sam = y_sam,
              EY_sam = EY_sam,
              probzero_sam = probzero_sam,
              rho_sample = rho_sample))
}

predict_ZI_beta_regression <- function(stan_fit, X.pred, X.orig, thinning = 10, a = 1, rho_modeled = FALSE, C = 1000, b = 0.5, right_censored = FALSE, min_rho = 0) {
  ### function to make predictions with zero-inflated left-censored beta regression
  # X.pred: prediction matrix (mxp)
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # thinning: use every thinning:th posterior sample in predictions
  # a: left-censoring constant
  # rho_modeled: TRUE: rho modeled with covariates, FALSE: common rho
  # C: if rho modeled, what is it's upper limit? rho = C*inv_logit(a+XB)
  # b: right-censoring constant
  # right_censored: Boolean to tell whether right-censoring is used
  # min_rho: minimum possible value for rho(x) (if rho is modeled with covariates)
  # RETURNS list of rep x m matrix of samples from posterior predictive for different quantities
  ### 1) latent linear predictors f 
  ### 2) predicted Ys
  ### 3) expected Ys
  ### 4) probabilies of zero
  ### 5) probability of suitability 
  ### 6) expectations given suitability
  ### 7) rhos
  
  # load the posterior samples
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  alpha_pi_sam <- as.matrix(stan_fit, pars = c("alpha_pi"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta_1","beta_2"))
  beta_pi_sam <- as.matrix(stan_fit, pars = c("beta_pi_1","beta_pi_2"))
  
  # prepare rho
  if (rho_modeled) {
    # rho(x) modeled with covariates
    alpha_rho_sam <- as.matrix(stan_fit, pars = c("alpha_rho"))
    beta_rho_sam <- as.matrix(stan_fit, pars = c("beta_rho_1","beta_rho_2"))
  } else {
    # every location has common rho
    rho_sam <- as.matrix(stan_fit, pars = c("rho"))
  }

  # prepare prediction matrix (scale, add second order terms)
  X.pred.scaled <- scale_covariates(X.orig,X.pred)
  Xpred <- add_second_order_terms(X.pred.scaled, colnames(X.pred))
  Xpred <- as.matrix(Xpred) # from df to matrix for matrix calculation
  
  # initialize matrices for output
  f_sam <- c()
  y_sam <- c()
  EY_sam <- c()
  probzero_sam <- c()
  prob_suit_sam <- c()
  EY_if_suitable_sam <- c()
  rho_sample <- c()
  
  # loop over posterior samples (take every thinning:th sample)
  for (i in seq(1,nrow(beta_sam),thinning)) {
    # corresponding coefficients
    alpha_i <- alpha_sam[i,]
    alpha_pi_i <- alpha_pi_sam[i,]
    beta_i <- beta_sam[i,]
    beta_pi_i <- beta_pi_sam[i,]
    
    # calculate rho
    if (rho_modeled){
      rho_i <- min_rho + (C-min_rho)*inv_logit(as.vector(alpha_rho_sam[i,]  + Xpred %*% beta_rho_sam[i,]))
    } else {
      rho_i <- rep(rho_sam[i,], nrow(Xpred)) # every location has common rho
    }
    
    # save rho for output
    rho_sample <- rbind(rho_sample, rho_i)

    # latent f (linear predictor) for mean of beta
    f_i <- as.vector(alpha_i + Xpred %*% beta_i)
    f_sam <- rbind(f_sam,f_i)
    
    # calculate the mean of beta using logit-link
    mu_i <- inv_logit(f_i)
    
    # latent f (linear predictor) for probability of suitability
    f_pi_i <- as.vector(alpha_pi_i + Xpred %*% beta_pi_i)
    
    # calculate probability of suitaility using logit-link
    pi_i <- inv_logit(f_pi_i)
    prob_suit_sam <- rbind(prob_suit_sam,pi_i)
    
    ### calculate expectations and probability of zeroes
    ### NOTE: mapply iterates over pairs of (mu,pi,rho) with common a
    EYi <- mapply(integrate_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
    EYi_if_suitable <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a)) #NOTE: this does not care about prob. of suitability
    probzero_i <- mapply(calculate_probzero_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
    
    # save for output
    EY_sam <- rbind(EY_sam,EYi)
    EY_if_suitable_sam <- rbind(EY_if_suitable_sam, EYi_if_suitable)
    probzero_sam <- rbind(probzero_sam, probzero_i)
    
    ### also predict Ys
    # sample location suitabilities
    Z_i <- rbinom(nrow(Xpred),1,pi_i) #0/1 vector to tell whether the location is suitable (1 => Y from beta) or unsuitable (0 => Y=0)
    # sample latent beta variables
    V_i <- rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
    
    # scale & censor
    if (!right_censored) {
      # scale to (-a,1), everything < 0 are treated as 0
      y_i <- sapply(V_i, function(v) (max(0,(a+1)*v - a)))
    } else {
      # if right-censoring is used, scale to (-a,1+b), everything < 0 are treated as 0, everything > 1 are treated as 1
      y_i <- sapply(V_i, function(v) (max(0,min(1,(a+b+1)*v - a))))
    }
    # if the location was unsuitable, automatically y=0
    y_i <- Z_i*y_i
    y_sam <- rbind(y_sam, y_i)
  }
  
  # return everything
  return(list(f_sam = f_sam,
              y_sam = y_sam,
              EY_sam = EY_sam,
              probzero_sam = probzero_sam,
              prob_suit_sam = prob_suit_sam,
              EY_if_suitable_sam = EY_if_suitable_sam,
              rho_sample = rho_sample))
}

prepare_P_matrix <- function(pred.locations, grid_center.locations) {
  ### prepares a sparse (n_locations X n_spatial_cells) matrix that indicates the spatial random effect block each point belongs to
  ### each rows sums to 1, 1 indicating the random effect block
  ### Given P matrix and vector z of n_spatial_cells spatial grid cells, Pz gives a vector of spatial random effects for each prediction location!
  # pred.locations: includes prediction locations as Nx2-matrix (coordinates in TM35FIN system)
  # grid_center.locations: includes the spatial grid centers as matrix (coordinates in TM35FIN system)
  ### RETURNS matrix P
  
  # turn the matrixes into terra vectors
  # this requires library(terra)
  grid_center.vect <- vect(grid_center.locations, geom = c("x","y"), crs = "EPSG:3067")
  pred_locs.vect <- vect(pred.locations, geom = c("x","y"), crs = "EPSG:3067")
  
  # find nearest spatial grid cell for each prediction location
  nearest_grid_cells.df <- as.data.frame(nearest(pred_locs.vect, grid_center.vect))
  # vector telling the ID of the closest grid cell
  nearest_id <- nearest_grid_cells.df$to_id 
  
  # initialize P as zeros
  P <- matrix(0, nrow = nrow(pred.locations), ncol = nrow(grid_center.locations))
  # loop over prediction locations, add 1 to the column corresponding to the closest spatial grid cell
  for (i in 1:nrow(pred.locations)) {
    P[i,nearest_id[i]] <- 1
  }
  
  # name rows and columns & return
  colnames(P) <- 1:nrow(grid_center.locations)
  rownames(P) <- 1:nrow(pred.locations)
  return(P)
}


predict_spatial_beta_regression <- function(stan_fit, X.pred, X.orig, pred.locations, S.pred, S.obs, thinning = 10, a = 1, rho_modeled = FALSE, C = 1000) {
  ### function to make predictions with spatial left-censored beta regression
  # X.pred: prediction matrix (mxp)
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # pred.locations: locations for each prediction cell (in meters)
  # S.pred: locations of coarse spatial grid cells (in kilometers)
  # S.obs: locations of observed coarse spatial grid cells (in kilometers, WATCH THAT THIS IS THE SAME AS WHEN FITTING A MODEL)
  # thinning: use every thinning:th posterior sample in predictions
  # a: left-censoring constant
  # rho_modeled: TRUE: rho modeled with covariates, FALSE: common rho
  # C: if rho modeled, what is it's upper limit? rho = C*inv_logit(a+XB)
  # RETURNS list of rep x m matrix of samples from posterior predictive for different quantities
  ### 1) latent linear predictors f 
  ### 2) predicted Ys
  ### 3) expected Ys
  ### 4) probabilies of zero
  ### 5) predicted spatial random effects
  ### 6) rhos
  
  # load in posterior samples
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta_1","beta_2"))
  
  # prepare rho
  if (rho_modeled) {
    # rho(x) modeled with covariates
    alpha_rho_sam <- as.matrix(stan_fit, pars = c("alpha_rho"))
    beta_rho_sam <- as.matrix(stan_fit, pars = c("beta_rho_1","beta_rho_2"))
  } else {
    # every location has common rho
    rho_sam <- as.matrix(stan_fit, pars = c("rho"))
  }
  
  # spatial random effects
  phi_sam <- as.matrix(stan_fit, pars = c("phi"))
  # covariance function parameters
  l_sam <- as.matrix(stan_fit, pars = c("l"))
  s2_sam <- as.matrix(stan_fit, pars = c("s2_cf"))
  
  # prepare prediction matrix (scale, add second order terms)
  X.pred.scaled <- scale_covariates(X.orig,X.pred)
  Xpred <- add_second_order_terms(X.pred.scaled, colnames(X.pred))
  Xpred <- as.matrix(Xpred) # from df to matrix for calculation
  
  # prepare P matrix (mxS) that tells in which spatial random effect cell each prediction point belongs to
  P <- prepare_P_matrix(pred.locations,S.pred*1000) # turn grid locations to meters
  
  # initialize matrices for output
  f_sam <- c()
  phi_pred_sam <- c()
  y_sam <- c()
  EY_sam <- c()
  probzero_sam <- c()
  rho_sample <- c()
  
  # loop over posterior samples (take every thinning:th sample)
  for (i in seq(1,nrow(beta_sam),thinning)) {
    # corresponding coefficients
    alpha_i <- alpha_sam[i,]
    beta_i <- beta_sam[i,]
    
    # calculate rho
    if (rho_modeled){
      rho_i <- C*inv_logit(as.vector(alpha_rho_sam[i,]  + Xpred %*% beta_rho_sam[i,]))
    } else {
      rho_i <- rep(rho_sam[i,], nrow(Xpred)) # every location has common rho
    }
    
    # save for output
    rho_sample <- rbind(rho_sample, rho_i)

    # latent f (linear predictor)
    f_i <- as.vector(alpha_i + Xpred %*% beta_i)
    
    # sample spatial random effect
    phi_obs_i <- phi_sam[i,] #these are spatial random effects that are sampled (included in the training data), correspond to locations S.obs
    l_i <- l_sam[i,]
    s2_i <- s2_sam[i,]
    
    # calculate covariance block-matrices
    K_pred_obs <- exp_covariance(S.pred,S.obs,s2_i,l_i)
    K_pred <- exp_covariance(S.pred,S.pred,s2_i,l_i)
    K_obs <- exp_covariance(S.obs,S.obs,s2_i,l_i)
    
    # mean and covariance for predicting phi in new locations given set of parameters and observed phi
    # justification for formulas can be read from Pietilä thesis (2025)
    phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_obs_i)
    phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
    
    ### sample spatial random effects over the prediction locations
    # add jitter on diagonal for computational stability
    phi_pred_Cov <- phi_pred_Cov + 1e-08*diag(length(phi_pred_m))
    
    # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
    L <- t(chol(phi_pred_Cov))
    phi_pred_i <- phi_pred_m + L %*% rnorm(length(phi_pred_m))
    
    # save spatial random effects
    phi_pred_sam <- rbind(phi_pred_sam, as.vector(P %*% phi_pred_i))
    
    # add the spatial random effect to linear predictor
    f_i <- f_i + as.vector(P %*% phi_pred_i) # P matrix picks the correct random effects (from the grid cell that the prediction point belongs to)
    f_sam <- rbind(f_sam,f_i)
    
    # mean of the beta distribution using logit-link
    mu_i <- inv_logit(f_i)
    
    ### calculate expectations and probability of zeroes
    ### NOTE: mapply iterates over pairs of (mu,rho) with common a
    EYi <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a))
    probzero_i <- mapply(calculate_probzero_LC_beta, mu=mu_i, rho=rho_i, MoreArgs = list(a=a))
    
    # save the expectations and prob. of zeros
    EY_sam <- rbind(EY_sam,EYi)
    probzero_sam <- rbind(probzero_sam, probzero_i)
    
    ### also predict Ys
    # latent beta variables
    V_i <- rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
    # scale & censor
    y_i <- sapply(V_i, function(v) (max(0,(a+1)*v - a)))
    y_sam <- rbind(y_sam, y_i)
  }
  
  # return everything
  return(list(f_sam = f_sam,
              y_sam = y_sam,
              EY_sam = EY_sam,
              probzero_sam = probzero_sam,
              phi_mu_pred_sam = phi_pred_sam,
              rho_sample = rho_sample))
}

predict_spatial_ZI_beta_regression <- function(stan_fit, X.pred, X.orig, pred.locations, S.pred, S.obs, thinning = 10, a = 1, rho_modeled = FALSE, C = 1000) {
  ### function to make predictions with zero-inflated spatial left-censored beta regression
  # X.pred: prediction matrix (mxp)
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # pred.locations: locations for each prediction cell (in meters)
  # S.pred: locations of coarse spatial grid cells (in kilometers)
  # S.obs: locations of observed coarse spatial grid cells (in kilometers, WATCH THAT THIS IS THE SAME AS WHEN FITTING A MODEL)
  # thinning: use every thinning:th posterior sample in predictions
  # a: left-censoring constant
  # rho_modeled: TRUE: rho modeled with covariates, FALSE: common rho
  # C: if rho modeled, what is it's upper limit? rho = C*inv_logit(a+XB)
  # RETURNS list of rep x m matrix of samples from posterior predictive for different quantities
  ### 1) latent linear predictors f 
  ### 2) predicted Ys
  ### 3) expected Ys
  ### 4) probabilies of zero
  ### 5) probability of suitability
  ### 6) predicted spatial random effects (for mean of beta distribution)
  ### 7) predicted spatial random effects (for probability of suitability)
  ### 8) rhos
  
  # load in posterior samples
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  alpha_pi_sam <- as.matrix(stan_fit, pars = c("alpha_pi"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta_1","beta_2"))
  beta_pi_sam <- as.matrix(stan_fit, pars = c("beta_pi_1","beta_pi_2"))
  
  # prepare rho
  if (rho_modeled) {
    # rho(x) modeled with covariates
    alpha_rho_sam <- as.matrix(stan_fit, pars = c("alpha_rho"))
    beta_rho_sam <- as.matrix(stan_fit, pars = c("beta_rho_1","beta_rho_2"))
  } else {
    # every location has common rho
    rho_sam <- as.matrix(stan_fit, pars = c("rho"))
  }
  
  # spatial random effects (mean of beta & prob. of suitability)
  phi_sam <- as.matrix(stan_fit, pars = c("phi_mu"))
  phi_pi_sam <- as.matrix(stan_fit, pars = c("phi_pi"))
  
  # covariance function parameters (for mean of beta)
  l_sam <- as.matrix(stan_fit, pars = c("l_mu"))
  s2_sam <- as.matrix(stan_fit, pars = c("s2_mu"))
  
  # covariance function parameters (for prob. of suitability)
  l_pi_sam <- as.matrix(stan_fit, pars = c("l_pi"))
  s2_pi_sam <- as.matrix(stan_fit, pars = c("s2_pi"))
  
  # prepare prediction matrix (scale, add second order terms)
  X.pred.scaled <- scale_covariates(X.orig,X.pred)
  Xpred <- add_second_order_terms(X.pred.scaled, colnames(X.pred))
  Xpred <- as.matrix(Xpred) # from df to matrix for calculation
  
  # prepare P matrix (mxS) that tells in which spatial random effect cell each prediction point belongs to
  P <- prepare_P_matrix(pred.locations,S.pred*1000) # turn grid locations to meters
  
  # initialize matrices for output
  f_sam <- c()
  y_sam <- c()
  EY_sam <- c()
  probzero_sam <- c()
  prob_suit_sam <- c()
  EY_if_suitable_sam <- c()
  phi_mu_pred_sam <- c()
  phi_pi_pred_sam <- c()
  rho_sample <- c()
  
  # loop over posterior samples (take every thinning:th sample)
  for (i in seq(1,nrow(beta_sam),thinning)) {
    # corresponding coefficients
    alpha_i <- alpha_sam[i,]
    alpha_pi_i <- alpha_pi_sam[i,]
    beta_i <- beta_sam[i,]
    beta_pi_i <- beta_pi_sam[i,]
    
    # calculate rho
    if (rho_modeled){
      rho_i <- C*inv_logit(as.vector(alpha_rho_sam[i,]  + Xpred %*% beta_rho_sam[i,]))
    } else {
      rho_i <- rep(rho_sam[i,], nrow(Xpred)) # every location has common rho
    }
    
    # save for output
    rho_sample <- rbind(rho_sample, rho_i)
    
    ### 1) PREDICT SPATIAL RANDOM EFFECTS FOR MEAN OF BETA DISTRIBUTION
    # latent f (linear predictor)
    f_i <- as.vector(alpha_i + Xpred %*% beta_i)
    
    phi_obs_i <- phi_sam[i,] #these are spatial random effects that are sampled (included in the training data), correspond to locations S.obs
    l_i <- l_sam[i,]
    s2_i <- s2_sam[i,]
    
    # calculate covariance block-matrices
    K_pred_obs <- exp_covariance(S.pred,S.obs,s2_i,l_i)
    K_pred <- exp_covariance(S.pred,S.pred,s2_i,l_i)
    K_obs <- exp_covariance(S.obs,S.obs,s2_i,l_i)
    
    # mean and covariance for predicting phi in new locations given set of parameters and observed phi
    # justification for formulas can be read from Pietilä thesis (2025)
    phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_obs_i)
    phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
    
    ### add jitter for computational stability
    phi_pred_Cov <- phi_pred_Cov + 1e-08*diag(length(phi_pred_m))
    
    # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
    L <- t(chol(phi_pred_Cov))
    phi_pred_i <- phi_pred_m + L %*% rnorm(length(phi_pred_m))
    
    # save spatial effects
    phi_mu_pred_sam <- rbind(phi_mu_pred_sam, as.vector(P %*% phi_pred_i))
    
    # add the spatial random effect to linear predictor
    f_i <- f_i + as.vector(P %*% phi_pred_i) # P matrix picks the correct random effects (from the grid cell that the prediction point belongs to)
    f_sam <- rbind(f_sam,f_i)
    
    # mean of the beta distribution using logit-link
    mu_i <- inv_logit(f_i)
    
    ### 2) PREDICT SPATIAL RANDOM EFFECTS FOR PROBABILITY OF ZERO
    # latent f (linear predictor)
    f_pi_i <- as.vector(alpha_pi_i + Xpred %*% beta_pi_i)
    
    # sample spatial random effects
    phi_obs_i <- phi_pi_sam[i,] #these are spatial random effects that are sampled (included in the training data), correspond to locations S.obs
    l_i <- l_pi_sam[i,]
    s2_i <- s2_pi_sam[i,]
    
    # calculate covariance block-matrices
    K_pred_obs <- exp_covariance(S.pred,S.obs,s2_i,l_i)
    K_pred <- exp_covariance(S.pred,S.pred,s2_i,l_i)
    K_obs <- exp_covariance(S.obs,S.obs,s2_i,l_i)
    
    # mean and covariance for predicting phi in new locations given set of parameters and observed phi
    # justification for formulas can be read from Pietilä thesis (2025)
    phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_obs_i)
    phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
    
    ### add jitter for computational stability
    phi_pred_Cov <- phi_pred_Cov + 1e-08*diag(length(phi_pred_m))
    
    # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
    L <- t(chol(phi_pred_Cov))
    phi_pi_pred_i <- phi_pred_m + L %*% rnorm(length(phi_pred_m))
    
    # save spatial effects
    phi_pi_pred_sam <- rbind(phi_pi_pred_sam, as.vector(P %*% phi_pi_pred_i))
    
    # add the spatial random effect to linear predictor
    f_pi_i <- f_pi_i + as.vector(P %*% phi_pi_pred_i) # P matrix picks the correct random effects (from the grid cell that the prediction point belongs to)

    ### probability of suitability using logit-link
    pi_i <- inv_logit(f_pi_i)
    prob_suit_sam <- rbind(prob_suit_sam, pi_i)
    
    ### calculate expectations and probability of zeroes
    ### NOTE: mapply iterates over pairs of (mu,pi,rho) with common a
    EYi <- mapply(integrate_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
    EYi_if_suitable <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a)) #NOTE: this does not care about prob. of suitability
    probzero_i <- mapply(calculate_probzero_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
    
    # save for output
    EY_sam <- rbind(EY_sam,EYi)
    EY_if_suitable_sam <- rbind(EY_if_suitable_sam, EYi_if_suitable)
    probzero_sam <- rbind(probzero_sam, probzero_i)
  
    ### also predict Ys
    # sample suitabilities of locations
    Z_i <- rbinom(nrow(Xpred),1,pi_i) #0/1 vector to tell whether the location is suitable (1 => Y from beta) or unsuitable (0 => Y=0)
    # sample latent beta variables
    V_i <- rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
    # scale & censor
    y_i <- sapply(V_i, function(v) (max(0,(a+1)*v - a)))
    # if the location was unsuitable, automatically y = 0
    y_i <- Z_i*y_i 
    y_sam <- rbind(y_sam, y_i)
  }
  
  # return everything
  return(list(f_sam = f_sam,
              y_sam = y_sam,
              EY_sam = EY_sam,
              probzero_sam = probzero_sam,
              prob_suit_sam = prob_suit_sam,
              EY_if_suitable_sam = EY_if_suitable_sam,
              phi_mu_pred_sam = phi_mu_pred_sam,
              phi_pi_pred_sam = phi_pi_pred_sam,
              rho_sample = rho_sample))
}
  
####### UNEXAMINED FROM THIS ONWARDS #################


### IS THIS THE SAME THAN USED IN ANALYSE_RESULTS?
plot_and_save_responses <- function(mod_list,X,grid_length = 200,im_width,im_height,thinning=200) {
  
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  post.mean_list <- list() #list with n_model components, each are grid_length x p matrices
  prob.zero_list <- list() #same for probability of zero
  
  mod_idx <- 1
  subfolder <- paste0("n_",nrow(X))
  for(mod in mod_list) {
    # save the posterior means
    post.mean_mat <- c() #grid_length x p matrix 
    prob.zero_mat <- c()
    
    # save all the predictions for current model
    all_vars_preds <- list() #list with p components, each n_rep x grid_length
    
    ### plots for separate curves if ZI model
    if (mod_idx %in% c(2,4)) {
      png(paste0("plots/final_results/",subfolder,"/M",mod_idx,"/separate_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
      par(mfrow = c(3,3),
          mar = c(5,4,2,2))
    }
    
    for (i in 1:ncol(X)) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      if (mod_idx %in% c(1,3)) {
        res <- predict_beta_regression(mod,predX,X,thinning,1)
      } else {
        res <- predict_ZI_beta_regression(mod,predX,X,thinning,1)
      }
      res.EY <- res$EY_sam # n_rep x grid_length
      all_vars_preds[[colnames(X)[i]]] <- res.EY
      post.means <- colMeans(res.EY)
      
      #plot the probability of presences and expectation given presence
      if (mod_idx %in% c(2,4)) {
        
        pi <- colMeans(res$prob_suit_sam)
        EY_if_suitable <- colMeans(res$EY_if_suitable_sam)
        
        plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), main = colnames(X)[i],type="l")
        lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
        if (i == 1){
          legend("bottomright", legend = c("Pr(suitable)","E[Y|suitable]"), lty = 1, col = c("blue","red"), lwd = 1, cex = 1)
        }
        # add legend
        #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
        #legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
      }
      
      # save posterior means
      post.mean_mat <- cbind(post.mean_mat, post.means)
      prob.zero_mat <- cbind(prob.zero_mat, colMeans(res$probzero_sam))
    }
    
    if (mod_idx %in% c(2,4)) {
      dev.off()
    }
    
    # now we want to normalize such that largest posterior men gets value 1!
    max_post_mean <- max(post.mean_mat)
    post.mean_mat <- post.mean_mat / max_post_mean #now max value is 1
    
    # save scaled posterior means
    post.mean_list[[mod_idx]] <- post.mean_mat
    prob.zero_list[[mod_idx]] <- prob.zero_mat
    
    # scale also all the other samples
    all_vars_preds <- lapply(all_vars_preds, function(x) (x/max_post_mean))
    
    
    # draw plots
    png(paste0("plots/final_results/",subfolder,"/M",mod_idx,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    
    par(mfrow = c(3,3),
        mar = c(5,4,2,2))
    
    for (i in 1:ncol(X)) { # index for variable
      plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
           xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i])
      
      n_rep <- nrow(res.EY)
      for (j in 1:n_rep) { # index for sample
        lines(grid_matrix[,i],all_vars_preds[[i]][j,], col = "lightgrey", lwd = 0.8)
      }
      lines(grid_matrix[,i], post.mean_mat[,i], col = "black", lwd = 1)
    }
    
    dev.off()
    
    mod_idx <- mod_idx + 1
  }
  
  ### then also plot the colmeans on the same plot
  cols <- rainbow(length(post.mean_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  
  png(paste0("plots/final_results/",subfolder,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i], lwd = 1)
    for (j in 1:length(post.mean_list)) {
      lines(grid_matrix[,i],post.mean_list[[j]][,i], col = cols[j], lwd = 1)
    }
    if (i == 1){
      legend("topright", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
  #legend("center", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1.5, bty = "n")
  dev.off()
  
  ### same for probability of zero?
  png(paste0("plots/final_results/",subfolder,"/prob_zero_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  
  cols <- rainbow(length(prob.zero_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = "value", ylab = "prob. of zero", main = colnames(X)[i], lwd = 1)
    for (j in 1:length(prob.zero_list)) {
      lines(grid_matrix[,i],prob.zero_list[[j]][,i], col = cols[j], lwd = 1)
    }
    if (i == 1){
      legend("topright", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
  #legend("center", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1.5, bty = "n")
  dev.off()
}

plot_responses <- function(mod,X,grid_length = 200,mod_idx,rho_modeled=FALSE,C=1000,thinning=200) {
  ### mod_idx = 1,2,3 or 4 (base,ZI,RE,RE+ZI)
  
  ### create a matrix that includes a grid from min to max observed values for each covariate
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  ### create a matrix that just repeats the mean value of each covariate
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  # save the posterior means
  post.mean_mat <- c() #grid_length x p matrix 
  prob.zero_mat <- c()
  
  # save all the predictions for current model
  all_vars_preds <- list() #list with p components, each n_rep x grid_length
  
  ### plots for separate curves if ZI model
  if (mod_idx %in% c(2,4)) {
    par(mfrow = c(3,3),
        mar = c(5,4,2,2))
  }
  
  for (i in 1:ncol(X)) {
    predX <- colmeans_matrix # take all the variables to their mean in the training data
    predX[,i] <- grid_matrix[,i] # replace one covariate by corresponding grid
    predX <- as.data.frame(predX)
    
    if (mod_idx %in% c(1,3)) {
      res <- predict_beta_regression(mod,predX,X,thinning,a=1,rho_modeled=rho_modeled,C=C)
    } else {
      res <- predict_ZI_beta_regression(mod,predX,X,thinning,a=1,rho_modeled=rho_modeled,C=C)
    }
    
    res.EY <- res$EY_sam # n_rep x grid_length
    all_vars_preds[[colnames(X)[i]]] <- res.EY
    post.means <- colMeans(res.EY)
    
    #plot the probability of presences and expectation given presence
    if (mod_idx %in% c(2,4)) {
      
      pi <- colMeans(res$prob_suit_sam)
      EY_if_suitable <- colMeans(res$EY_if_suitable_sam)
      
      plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), main = colnames(X)[i],type="l")
      lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
      if (i == 1){
        legend("bottomright", legend = c("Pr(suitable)","E[Y|suitable]"), lty = 1, col = c("blue","red"), lwd = 1, cex = 1)
      }
      # add legend
      #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
      #legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
    }
    
    # save posterior means
    post.mean_mat <- cbind(post.mean_mat, post.means)
    prob.zero_mat <- cbind(prob.zero_mat, colMeans(res$probzero_sam))
  }
  
  # # now we want to normalize such that largest posterior men gets value 1!
  # max_post_mean <- max(post.mean_mat)
  # post.mean_mat <- post.mean_mat / max_post_mean #now max value is 1
  # 
  # # save scaled posterior means
  # #post.mean_list[[mod_idx]] <- post.mean_mat
  # #prob.zero_list[[mod_idx]] <- prob.zero_mat
  # 
  # # scale also all the other samples
  # all_vars_preds <- lapply(all_vars_preds, function(x) (x/max_post_mean))
  # 
  # 
  # # draw plots
  # 
  # par(mfrow = c(3,3),
  #     mar = c(5,4,2,2))
  # 
  # for (i in 1:ncol(X)) { # index for variable
  #   plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
  #        xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i])
  #   
  #   n_rep <- nrow(res.EY)
  #   for (j in 1:n_rep) { # index for sample
  #     lines(grid_matrix[,i],all_vars_preds[[i]][j,], col = "lightgrey", lwd = 0.8)
  #   }
  #   lines(grid_matrix[,i], post.mean_mat[,i], col = "black", lwd = 1)
  # }
}


plot_responses_one_model <- function(mod,X,grid_length = 200,mod_idx,rho_modeled=FALSE,C=1000,thinning=200) {
  ### create a matrix that includes a grid from min to max observed values for each covariate
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  ### create a matrix that just repeats the mean value of each covariate
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  
}
  
  
plot_responses_beta_regression <- function(stan_fit, X.train, X.orig, is_standardized = TRUE, second_order = TRUE, xmin=-3, xmax=3, ymin=0, ymax=1,
                                           plot_nx=3, plot_ny=3,a=1,type="expectation") {
  ###
  # xmin: min depth
  # xmax: max depth
  # type: either "expectation" or "probzero"
  
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta_1","beta_2"))
  rho_sam <- as.matrix(stan_fit, pars = c("rho"))
  
  ### response curves
  par(mfrow = c(plot_nx,plot_ny))
  
  ### remove the second order terms since they will be recreated for prediction values
  if (second_order) {
    X.train <- X.train[,-grep("\\^2", colnames(X.train))]
  }
  
  ### prepare a grid matrix
  #depth_idx <- which(colnames(X.train) == "depth")
  x_grid <- c()
  for (i in 1:ncol(X.train)) {
    #x_grid <- cbind(x_grid, seq(0,100,length=100))
    x_grid <- cbind(x_grid, seq(min(X.orig[,i]),max(X.orig[,i]),length=1000))
  }
  #x_grid[,depth_idx] <- seq(xmin,xmax,length=100)
  x_grid_orig <- x_grid # original scale grid
  
  ### scale if the X.train also standardized
  if (is_standardized) {
    x_grid <- scale_covariates(X.orig,x_grid)
  }
  
  colnames(x_grid) <- colnames(X.train)
  x_grid <- as.data.frame(x_grid)
  
  ### add second order if they are in the model
  if (second_order) {
    x_grid <- add_second_order_terms(x_grid, colnames(x_grid))
  }
  
  ### go through the variables
  for (i in 1:ncol(X.train)) {
    var_name <- colnames(X.train)[i]
    f_i <- mean(alpha_sam) + x_grid[,i] * mean(beta_sam[,i])
    if (second_order) {
      f_i <- f_i + x_grid[,i+ncol(X.train)] * mean(beta_sam[,i+ncol(X.train)])
    }
    
    mu_i <- inv_logit(f_i)
    rho_mean <- mean(rho_sam)
    
    ### calculate the expected values for every grid point
    EYi <- c()
    grid_01 <- seq(0.01,1-0.01,0.01)
    
    for(j in 1:length(x_grid[,i])) {
      #EYi <- c(EYi, sum(0.01*dbeta((grid_01+a)/(a+1),mu_i[j]*phi_mean,(1-mu_i[j])*phi_mean)/(a+1)))
      if (type == "expectation") {
        calc_density <- function(x) (x*dbeta((x+a)/(a+1),mu_i[j]*rho_mean,(1-mu_i[j])*rho_mean)/(a+1))
        EYi <- c(EYi, integrate(calc_density,0,1)$value)
      } else {
        EYi <- c(EYi, pbeta(a/(a+1),mu_i[j]*rho_mean,(1-mu_i[j])*rho_mean)) #probability of zero
      }
    }
    
    y_str <- "Expected cover"
    if(type != "expectation") {
      y_str <- "P(Y=0)"
    }
    plot(x_grid_orig[,i],EYi,type="l",
         ylab=y_str,xlab="value",
         main=colnames(X.train)[i],ylim = c(ymin,ymax))
  }
}



plot_responses_ZI_beta_regression <- function(stan_fit, X.train, X.orig, is_standardized = TRUE, second_order = TRUE, xmin=-3, xmax=3, ymin=0, ymax=1,
                                              plot_nx=3, plot_ny=3,a=1,type="expectation") {
  ###
  # xmin: min depth
  # xmax: max depth
  # type: either "expectation" or "probzero"
  
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  alpha_pi_sam <- as.matrix(stan_fit, pars = c("alpha_pi"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta"))
  beta_pi_sam <- as.matrix(stan_fit, pars = c("beta_pi"))
  rho_sam <- as.matrix(stan_fit, pars = c("rho"))
  
  ### response curves
  par(mfrow = c(plot_nx,plot_ny))
  
  ### remove the second order terms since they will be recreated for prediction values
  if (second_order) {
    X.train <- X.train[,-grep("\\^2", colnames(X.train))]
  }
  
  ### prepare a grid matrix
  depth_idx <- which(colnames(X.train) == "depth")
  x_grid <- c()
  for (i in 1:ncol(X.train)) {
    #x_grid <- cbind(x_grid, seq(0,100,length=100))
    x_grid <- cbind(x_grid, seq(min(X.orig[,i]),max(X.orig[,i]),length=100))
  }
  #x_grid[,depth_idx] <- seq(xmin,xmax,length=100)
  x_grid_orig <- x_grid # original scale grid
  
  ### scale if the X.train also standardized
  if (is_standardized) {
    x_grid <- scale_covariates(X.orig,x_grid)
  }
  
  colnames(x_grid) <- colnames(X.train)
  x_grid <- as.data.frame(x_grid)
  
  ### add second order if they are in the model
  if (second_order) {
    x_grid <- add_second_order_terms(x_grid, colnames(x_grid))
  }
  
  ### go through the variables
  for (i in 1:ncol(X.train)) {
    var_name <- colnames(X.train)[i]
    
    ### latent variables
    f_mu <- mean(alpha_sam) + x_grid[,i] * mean(beta_sam[,i])
    f_pi <- mean(alpha_pi_sam) + x_grid[,i] * mean(beta_pi_sam[,i])
    
    if (second_order) {
      f_mu <- f_mu + x_grid[,i+ncol(X.train)] * mean(beta_sam[,i+ncol(X.train)])
      f_pi <- f_pi + x_grid[,i+ncol(X.train)] * mean(beta_pi_sam[,i+ncol(X.train)])
    }
    
    mu <- inv_logit(f_mu)
    prob_suitability <- inv_logit(f_pi)
    rho_mean <- mean(rho_sam)
    
    ### calculate the expected values for every grid point
    EYi <- c()
    #grid_01 <- seq(0.01,1-0.01,0.01)
    
    for(j in 1:length(x_grid[,i])) {
      #EYi <- c(EYi, sum(0.01*dbeta((grid_01+a)/(a+1),mu_i[j]*phi_mean,(1-mu_i[j])*phi_mean)/(a+1)))
      if (type == "expectation") {
        calc_density <- function(x) (x*prob_suitability[j]*dbeta((x+a)/(a+1),mu[j]*rho_mean,(1-mu[j])*rho_mean)/(a+1))
        EYi <- c(EYi, integrate(calc_density,0.001,1-0.001)$value)
      } else {
        EYi <- c(EYi, 1-prob_suitability[j] + prob_suitability[j]*pbeta(a/(a+1),mu[j]*rho_mean,(1-mu[j])*rho_mean))
      }
    }
    
    y_str <- "Expected cover"
    if (type != "expectation") {
      y_str <- "Pr(Y=0)"
    }
    plot(x_grid_orig[,i],EYi,type="l",
         ylab=y_str,xlab="value",
         main=colnames(X.train)[i],ylim = c(ymin,ymax))
  }
}

plot_separate_responses_ZI_beta_regression <- function(stan_fit, X.train, X.orig, is_standardized = TRUE, second_order = TRUE, xmin=-3, xmax=3, ymin=0, ymax=1, plot_nx=3, plot_ny=3,a=1) {
  ###
  # xmin: min depth
  # xmax: max depth
  
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  alpha_pi_sam <- as.matrix(stan_fit, pars = c("alpha_pi"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta"))
  beta_pi_sam <- as.matrix(stan_fit, pars = c("beta_pi"))
  rho_sam <- as.matrix(stan_fit, pars = c("rho"))
  
  ### response curves
  par(mfrow = c(plot_nx,plot_ny))
  
  ### remove the second order terms since they will be recreated for prediction values
  if (second_order) {
    X.train <- X.train[,-grep("\\^2", colnames(X.train))]
  }
  
  ### prepare a grid matrix
  depth_idx <- which(colnames(X.train) == "depth")
  x_grid <- c()
  for (i in 1:ncol(X.train)) {
    #x_grid <- cbind(x_grid, seq(0,100,length=100))
    x_grid <- cbind(x_grid, seq(min(X.orig[,i]),max(X.orig[,i]),length=100))
  }
  #x_grid[,depth_idx] <- seq(xmin,xmax,length=100)
  x_grid_orig <- x_grid # original scale grid
  
  ### scale if the X.train also standardized
  if (is_standardized) {
    x_grid <- scale_covariates(X.orig,x_grid)
  }
  
  colnames(x_grid) <- colnames(X.train)
  x_grid <- as.data.frame(x_grid)
  
  ### add second order if they are in the model
  if (second_order) {
    x_grid <- add_second_order_terms(x_grid, colnames(x_grid))
  }
  
  ### go through the variables
  for (i in 1:ncol(X.train)) {
    var_name <- colnames(X.train)[i]
    
    ### latent variables
    f_mu <- mean(alpha_sam) + x_grid[,i] * mean(beta_sam[,i])
    f_pi <- mean(alpha_pi_sam) + x_grid[,i] * mean(beta_pi_sam[,i])
    
    if (second_order) {
      f_mu <- f_mu + x_grid[,i+ncol(X.train)] * mean(beta_sam[,i+ncol(X.train)])
      f_pi <- f_pi + x_grid[,i+ncol(X.train)] * mean(beta_pi_sam[,i+ncol(X.train)])
    }
    
    mu <- inv_logit(f_mu)
    prob_suitability <- inv_logit(f_pi)
    rho_mean <- mean(rho_sam)
    
    ### calculate the expected values for every grid point
    EYi <- c()
    grid_01 <- seq(0.01,1-0.01,0.01)
    
    for(j in 1:length(x_grid[,i])) {
      #EYi <- c(EYi, sum(0.01*dbeta((grid_01+a)/(a+1),mu_i[j]*phi_mean,(1-mu_i[j])*phi_mean)/(a+1)))
      calc_density <- function(x) (x*dbeta((x+a)/(a+1),mu[j]*rho_mean,(1-mu[j])*rho_mean)/(a+1))
      EYi <- c(EYi, integrate(calc_density,0.001,1-0.001)$value)
    }
    
    plot(NULL, xlim = c(min(x_grid_orig[,i]),max(x_grid_orig[,i])), ylim = 0:1,
         ylab="",xlab="value",
         main=colnames(X.train)[i])
    
    lines(x_grid_orig[,i],EYi,col="red")
    lines(x_grid_orig[,i],prob_suitability,col="blue")
  }
  
  ### create a legend
  plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  legend("center", legend = c("E[Coverage|Suitable]","Pr(Suitable)"), col = c("red","blue"), lty = 1, bty = "n")
}


average_distribution_beta_regression <- function(mod,a=1,sp_name) {
  par(mfrow = c(1,1),
      mar = c(4,4,2,0))
  
  x_grid <- seq(-a+0.01,1-0.01,length=500)
  alpha_sam <- as.matrix(mod, pars = c("alpha"))
  alpha_mean <- mean(alpha_sam)
  mu_mean <- inv_logit(alpha_mean)
  rho_sam <- as.matrix(mod, pars = c("rho"))
  rho_mean <- mean(rho_sam)
  
  plot(x_grid, dbeta((x_grid+a)/(a+1),mu_mean*rho_mean, (1-mu_mean)*rho_mean)/(a+1),
       type = "l", ylab = "p(y)", xlab = "y",
       main = paste0("Beta distribution of coverage in avg. conditions (", sp_name,")"))
  
  ### end point for polygon based on the trend 
  end_x <- 0
  end_y <- 0
  if (dbeta(0.01,mu_mean*rho_mean,(1-mu_mean)*rho_mean) > dbeta(0.5,mu_mean*rho_mean,(1-mu_mean)*rho_mean)) {
    end_x <- c(0,-1)
    end_y <- c(0,0)
  }
  
  polygon(c(x_grid[x_grid <= 0], end_x),
          c(dbeta((a+x_grid[x_grid <= 0])/(a+1), mu_mean*rho_mean, (1-mu_mean)*rho_mean)/(a+1), end_y), col = "lightblue")
  
  prob_y0 <- pbeta(a/(a+1),mu_mean*rho_mean,(1-mu_mean)*rho_mean)
  
  text(-0.5,1,paste0("Pr(Y=0) = ",round(prob_y0,2)))
}

average_distribution_ZI_beta_regression <- function(mod,a=1,sp_name) {
  par(mfrow = c(1,1),
      mar = c(4,4,2,0))
  
  alpha_sam <- as.matrix(mod, pars = c("alpha"))
  alpha_mean <- mean(alpha_sam)
  mu_mean <- inv_logit(alpha_mean)
  
  alpha_pi_sam <- as.matrix(mod, pars = c("alpha_pi"))
  alpha_pi_mean <- mean(alpha_pi_sam)
  pi_mean <- inv_logit(alpha_pi_mean)
  
  rho_sam <- as.matrix(mod, pars = c("rho"))
  rho_mean <- mean(rho_sam)
  
  x_grid <- seq(-a+0.01,1-0.01,length=500)
  
  plot(x_grid, dbeta((x_grid+a)/(a+1),mu_mean*rho_mean, (1-mu_mean)*rho_mean)/(a+1),
       type = "l", ylab = "p(y|suitable)", xlab = "y",
       main = paste0("Beta distribution of coverage in avg. conditions (", sp_name,")"))
  
  ### end point for polygon based on the trend 
  end_x <- 0
  end_y <- 0
  if (dbeta(0.01,mu_mean*rho_mean,(1-mu_mean)*rho_mean) > dbeta(0.5,mu_mean*rho_mean,(1-mu_mean)*rho_mean)) {
    end_x <- c(0,-1)
    end_y <- c(0,0)
  }

  polygon(c(x_grid[x_grid <= 0], end_x),
          c(dbeta((a+x_grid[x_grid <= 0])/(a+1), mu_mean*rho_mean, (1-mu_mean)*rho_mean)/(a+1), end_y), col = "lightblue")
  
  prob_y0 <- pbeta(a/(a+1),mu_mean*rho_mean,(1-mu_mean)*rho_mean)
  
  text(-0.5,1,paste0("Pr(Y=0 | Suitable) = ",round(prob_y0,2)))
  text(-0.5,1.5,paste0("Pr(Suitable) = ", round(pi_mean,2)))
}

average_distribution_binomial_regression <- function(mod,sp_name) {
  par(mfrow = c(1,1),
      mar = c(4,4,2,0))
  
  y_grid <- seq(0,100,1)
  alpha_sam <- as.matrix(mod, pars = c("alpha"))
  plot(y_grid, dbinom(y_grid,100,inv_logit(mean(alpha_sam))), type = "h", lwd = 3,
       xlab = "y", ylab = "Pr(y)",
       main = paste0("Binomial distribution of coverage in average conditions (",sp_name,")"))
}

plot_spatial_effects_beta <- function(mod,s_obs,s_pred,grid.vect,obs.vect,shoreline.vect) {
  ###
  ###
  ###
  ###
  
  alpha.sam <- as.matrix(mod, pars = c("alpha"))
  beta.sam <- as.matrix(mod, pars = c("beta"))
  rho.sam <- as.matrix(mod, pars = c("rho"))
  phi.sam <- as.matrix(mod, pars = c("phi"))
  
  # take in the random effects from observed grid cells
  phi_obs <- colMeans(phi.sam)
  
  # means of the parameters
  s2_cf.sam <- as.matrix(mod_amphi.beta_spat, pars = "s2_cf")
  l.sam <- as.matrix(mod_amphi.beta_spat, pars = "l")
  s2_cf_mean <- mean(s2_cf.sam)
  l_mean <- mean(l.sam)
  
  # covariance matrixes with set of parameters
  K_pred_obs <- exp_covariance(s_pred,s_obs,s2_cf_mean,l_mean)
  K_pred <- exp_covariance(s_pred,s_pred,s2_cf_mean,l_mean)
  K_obs <- exp_covariance(s_obs,s_obs,s2_cf_mean,l_mean)
  
  # mean and covariance for predicting phi in new locations given set of parameters and observed phi
  phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_obs)
  phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
  
  # predictive mean
  grid.vect$spat_eff <- phi_pred_m
  
  grid.rast <- rast(ext = ext(grid.vect), crs = "EPSG:3067")
  grid.rast <- rasterize(grid.vect,grid.rast,field = "spat_eff")
  
  par(mfrow = c(1,1))
  plot(grid.rast, main = "prediction grid")
  plot(obs.vect, add = TRUE, col = "red")
  plot(shoreline.vect, add = TRUE, col = "lightgrey")
} 

plot_spatial_effects_ZIBeta <- function(mod,s_obs,s_pred,grid.vect,obs.vect,shoreline.vect) {
  ###
  ###
  ###
  ###
  
  # gather the parameters
  alpha.sam <- as.matrix(mod, pars = c("alpha"))
  alpha_pi.sam <- as.matrix(mod, pars = c("alpha_pi"))
  beta.sam <- as.matrix(mod, pars = c("beta"))
  beta_pi.sam <- as.matrix(mod, pars = c("beta_pi"))
  rho.sam <- as.matrix(mod, pars = c("rho"))
  phi_mu.sam <- as.matrix(mod, pars = c("phi_mu"))
  phi_pi.sam <- as.matrix(mod, pars = c("phi_pi"))
  
  ### 1) random effect for mu parameter
  # take in the random effects from observed grid cells
  phi_mu.obs <- colMeans(phi_mu.sam)
  
  # post. means of the parameters
  s2_mu.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = "s2_mu")
  l_mu.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = "l_mu")
  s2_mu_mean <- mean(s2_mu.sam)
  l_mu_mean <- mean(l_mu.sam)
  
  # covariance matrixes with set of parameters
  K_pred_obs <- exp_covariance(s_pred,s_obs,s2_mu_mean,l_mu_mean)
  K_pred <- exp_covariance(s_pred,s_pred,s2_mu_mean,l_mu_mean)
  K_obs <- exp_covariance(s_obs,s_obs,s2_mu_mean,l_mu_mean)
  
  # mean and covariance for predicting phi in new locations given set of parameters and observed phi
  phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_mu.obs)
  phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
  
  grid.vect$spat_eff <- phi_pred_m
  
  grid.rast <- rast(ext = ext(grid.vect), crs = "EPSG:3067")
  grid.rast <- rasterize(grid.vect,grid.rast,field = "spat_eff")
  
  par(mfrow = c(1,1))
  plot(grid.rast, main = "spatial effects for mu")
  plot(obs.vect, add = TRUE, col = "red")
  plot(shoreline.vect, add = TRUE, col = "lightgrey")
  
  ### 2) random effect of the prob. of suitability
  # take in the random effects from observed grid cells
  phi_pi.obs <- colMeans(phi_pi.sam)
  
  # posterior means
  s2_pi.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = "s2_pi")
  l_pi.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = "l_pi")
  s2_pi_mean <- mean(s2_pi.sam)
  l_pi_mean <- mean(l_pi.sam)
  
  # covariance matrixes with set of parameters
  K_pred_obs <- exp_covariance(s_pred,s_obs,s2_pi_mean,l_pi_mean)
  K_pred <- exp_covariance(s_pred,s_pred,s2_pi_mean,l_pi_mean)
  K_obs <- exp_covariance(s_obs,s_obs,s2_pi_mean,l_pi_mean)
  
  # mean and covariance for predicting phi in new locations given set of parameters and observed phi
  phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_pi.obs)
  phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
  
  grid.vect$spat_eff <- phi_pred_m
  
  grid.rast <- rast(ext = ext(grid.vect), crs = "EPSG:3067")
  grid.rast <- rasterize(grid.vect,grid.rast,field = "spat_eff")
  
  par(mfrow = c(1,1))
  plot(grid.rast, main = "spatial effects for pi")
  plot(obs.vect, add = TRUE, col = "red")
  plot(shoreline.vect, add = TRUE, col = "lightgrey")
} 

### PLOTTING FUNCTIONS (THESE HAVE TO BE CLEANED!!!)

plot_map <- function(pred_list, locs, pred.grid.vect, type, title_chr = "", hotspot_proportion = 0.7) {
  ### type: "cover", "hotspot", 
  
  vals <- colMeans(pred_list$EY_sam)
  vars <- apply(pred_list$EY_sam,2,var)
  test_df <- as.data.frame(cbind(locs, vals, vars))
  
  ### accumulative
  test_df_ordered <- test_df[order(-test_df$vals), ]
  test_df_ordered$hotspot <- 1
  
  # find when 70% of the expected value is reached
  tot_sum <- sum(test_df_ordered$vals)
  idx <- which(cumsum(test_df_ordered$vals) > hotspot_proportion*tot_sum)[1]
  
  # set all outside to be non-hotspot
  test_df_ordered[idx:nrow(test_df_ordered),"hotspot"] <- 0
  
  test_vect <- vect(test_df_ordered, geom=c("x","y"), crs = "EPSG:3067")
  test_rast <- rast(ext = ext(pred.grid.vect), res = 1000, crs = "EPSG:3067")
  r <- rasterize(test_vect, test_rast, field = "vals")
  r.vars <- rasterize(test_vect, test_rast, field = "vars")
  r.hotspot <- rasterize(test_vect, test_rast, field = "hotspot")
  
  if (type == "cover") {
    plot(r, colNA = "lightgrey", main = title_chr, ,plg = list(title = "E[coverage]"),
         xlab = "Easting (m)", ylab = "Northing (m)")
  } else if (type == "hotspot") {
    plot(r.hotspot, colNA = "lightgrey", main = title_chr, col = c("red","blue"),
         #plg = list(title = "hotspot", legend = c("no","yes")),
         legend = FALSE,
         xlab = "Easting (m)", ylab = "Northing (m)")
    #legend(583000,6362000, legend = c("yes","no"), title = "hotspot", col = c("blue","red"), pch = 15,cex = 1)
    legend("bottomright", legend = c("yes","no"), title = "hotspot", col = c("blue","red"), pch = 15,cex = 1, inset = c(0.001,0.04))
  }
}

plot_map_differences <- function(pred_list1,pred_list2, locs, pred.grid.vect, type, title_chr = "", hotspot_proportion = 0.8, return = FALSE) {
  ### type: "cover", "hotspot", 
  
  vals1 <- colMeans(pred_list1$EY_sam)
  vals2 <- colMeans(pred_list2$EY_sam)
  diff_EY <- vals1-vals2
  mean_square_diff <- mean(diff_EY^2)
  
  ### produce hotspots for model 1
  test_df1 <- as.data.frame(cbind(locs, vals1))
  
  ### accumulative
  test_df1_ordered <- test_df1[order(-test_df1$vals1), ]
  test_df1_ordered$hotspot <- 1
  
  # find when 70% of the expected value is reached
  tot_sum <- sum(test_df1_ordered$vals1)
  idx <- which(cumsum(test_df1_ordered$vals1) > hotspot_proportion*tot_sum)[1]
  
  # set all outside to be non-hotspot
  test_df1_ordered[idx:nrow(test_df1_ordered),"hotspot"] <- 0
  test_df1$hotspot <- test_df1_ordered$hotspot[rank(-test_df1$vals1)] #find the hotspot values
  
  ### produce hotspots for model 2
  test_df2 <- as.data.frame(cbind(locs, vals2))
  
  ### accumulative
  test_df2_ordered <- test_df2[order(-test_df2$vals2), ]
  test_df2_ordered$hotspot <- 1
  
  # find when 70% of the expected value is reached
  tot_sum <- sum(test_df2_ordered$vals2)
  idx <- which(cumsum(test_df2_ordered$vals2) > hotspot_proportion*tot_sum)[1]
  
  # set all outside to be non-hotspot
  test_df2_ordered[idx:nrow(test_df2_ordered),"hotspot"] <- 0
  test_df2$hotspot <- test_df2_ordered$hotspot[rank(-test_df2$vals2)] #find the correct hotspot values
  
  ### calculate Jaccard Index
  n_both_hotspots <- sum((test_df1$hotspot == 1) & (test_df2$hotspot == 1))
  n_mod1_hotspots <- sum(test_df1$hotspot == 1)
  n_mod2_hotspots <- sum(test_df2$hotspot == 1)
  jaccard_idx <- n_both_hotspots / (n_mod1_hotspots + n_mod2_hotspots - n_both_hotspots)
  
  shared_hotspot <- test_df1$hotspot == test_df2$hotspot
  prop_agreement <- mean(shared_hotspot)
  
  test_df <- as.data.frame(cbind(locs,diff_EY,shared_hotspot))
  
  test_vect <- vect(test_df, geom=c("x","y"), crs = "EPSG:3067")
  test_rast <- rast(ext = ext(pred.grid.vect), res = 1000, crs = "EPSG:3067")
  r <- rasterize(test_vect, test_rast, field = "diff_EY")
  r.hotspot <- rasterize(test_vect, test_rast, field = "shared_hotspot")
  
  if (type == "cover") {
    plot(r, colNA = "lightgrey", main = title_chr ,plg = list(title = "difference"),
         xlab = "Easting (m)", ylab = "Northing (m)")
  } else if (type == "hotspot") {
    plot(r.hotspot, colNA = "lightgrey", main = paste0(title_chr," (J=",signif(jaccard_idx,3),")"), col = c("red","forestgreen"),
         #plg = list(title = "agreement", legend = c("no","yes"), cex = 1),
         legend = FALSE,
         xlab = "Easting (m)", ylab = "Northing (m)")
    #legend(572000,6362000, legend = c("yes","no"), title = "agreement", col = c("forestgreen","red"), pch = 15,cex = 1)
    legend("bottomright", legend = c("yes","no"), title = "agreement", col = c("forestgreen","red"), pch = 15,cex = 1, inset = c(0.001,0.04))
  }
  
  if(return == TRUE) {
    return(jaccard_idx)
  }
}

plot_random_effects <- function(pred_list, locs, pred.grid.vect, type, title_chr = "") {
  ### type: "mu" or "pi"
  
  if (type == "mu") {
    vals <- colMeans(pred_list$phi_mu_pred_sam)
    vars <- apply(pred_list$phi_mu_pred_sam,2,var)
  } else if (type == "pi") {
    vals <- colMeans(pred_list$phi_pi_pred_sam)
    vars <- apply(pred_list$phi_pi_pred_sam,2,var)
  }
  test_df <- as.data.frame(cbind(locs, vals, vars))
  
  test_vect <- vect(test_df, geom=c("x","y"), crs = "EPSG:3067")
  test_rast <- rast(ext = ext(pred.grid.vect), res = 1000, crs = "EPSG:3067")
  r <- rasterize(test_vect, test_rast, field = "vals")
  r.vars <- rasterize(test_vect, test_rast, field = "vars")
  
  if (type == "mu") {
    title_str <- "spatial RE"
  }
  
  plot(r, colNA = "lightgrey", main = title_chr, ,plg = list(title = paste0("spatial RE")),
       xlab = "Easting (m)", ylab = "Northing (m)")
}


plot_and_save_responses <- function(mod_list,X,grid_length = 200,im_width,im_height, thinning = 20) {
  
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  post.mean_list <- list() #list with n_model components, each are grid_length x p matrices
  prob.zero_list <- list() #same for probability of zero
  
  mod_idx <- 1
  subfolder <- paste0("n_",nrow(X))
  for(mod in mod_list) {
    # save the posterior means
    post.mean_mat <- c() #grid_length x p matrix 
    prob.zero_mat <- c()
    
    # save all the predictions for current model
    all_vars_preds <- list() #list with p components, each n_rep x grid_length
    
    ### plots for separate curves if ZI model
    if (mod_idx %in% c(2,4)) {
      png(paste0("plots/final_results/",subfolder,"/M",mod_idx,"/separate_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
      par(mfrow = c(3,3),
          mar = c(5,4,2,2))
    }
    
    for (i in 1:ncol(X)) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      if (mod_idx %in% c(1,3)) {
        res <- predict_beta_regression(mod,predX,X,thinning,1)
      } else {
        res <- predict_ZI_beta_regression(mod,predX,X,thinning,1)
      }
      res.EY <- res$EY_sam # n_rep x grid_length
      all_vars_preds[[colnames(X)[i]]] <- res.EY
      post.means <- colMeans(res.EY)
      
      #plot the probability of presences and expectation given presence
      if (mod_idx %in% c(2,4)) {
        
        pi <- colMeans(res$prob_suit_sam)
        EY_if_suitable <- colMeans(res$EY_if_suitable_sam)
        
        plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), main = colnames(X)[i],type="l")
        lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
        if (i == 1){
          legend("bottomright", legend = c("Pr(suitable)","E[Y|suitable]"), lty = 1, col = c("blue","red"), lwd = 1, cex = 1)
        }
        # add legend
        #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
        #legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
      }
      
      # save posterior means
      post.mean_mat <- cbind(post.mean_mat, post.means)
      prob.zero_mat <- cbind(prob.zero_mat, colMeans(res$probzero_sam))
    }
    
    if (mod_idx %in% c(2,4)) {
      dev.off()
    }
    
    # now we want to normalize such that largest posterior men gets value 1!
    max_post_mean <- max(post.mean_mat)
    post.mean_mat <- post.mean_mat / max_post_mean #now max value is 1
    
    # save scaled posterior means
    post.mean_list[[mod_idx]] <- post.mean_mat
    prob.zero_list[[mod_idx]] <- prob.zero_mat
    
    # scale also all the other samples
    all_vars_preds <- lapply(all_vars_preds, function(x) (x/max_post_mean))
    
    
    # draw plots
    png(paste0("plots/final_results/",subfolder,"/M",mod_idx,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    
    par(mfrow = c(3,3),
        mar = c(5,4,2,2))
    
    for (i in 1:ncol(X)) { # index for variable
      plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
           xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i])
      
      n_rep <- nrow(res.EY)
      for (j in 1:n_rep) { # index for sample
        lines(grid_matrix[,i],all_vars_preds[[i]][j,], col = "lightgrey", lwd = 0.8)
      }
      lines(grid_matrix[,i], post.mean_mat[,i], col = "black", lwd = 1)
    }
    
    dev.off()
    
    mod_idx <- mod_idx + 1
  }
  
  ### then also plot the colmeans on the same plot
  cols <- rainbow(length(post.mean_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  
  png(paste0("plots/final_results/",subfolder,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i], lwd = 1)
    for (j in 1:length(post.mean_list)) {
      lines(grid_matrix[,i],post.mean_list[[j]][,i], col = cols[j], lwd = 1)
    }
    if (i == 1){
      legend("topright", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
  #legend("center", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1.5, bty = "n")
  dev.off()
  
  ### same for probability of zero?
  png(paste0("plots/final_results/",subfolder,"/prob_zero_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  
  cols <- rainbow(length(prob.zero_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = "value", ylab = "prob. of zero", main = colnames(X)[i], lwd = 1)
    for (j in 1:length(prob.zero_list)) {
      lines(grid_matrix[,i],prob.zero_list[[j]][,i], col = cols[j], lwd = 1)
    }
    if (i == 1){
      legend("topright", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
  #legend("center", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1.5, bty = "n")
  dev.off()
}


plot_and_save_responses_v2 <- function(mod_list,mod_rho_list,X,grid_length = 200,im_width,im_height, thinning = 20, use_median = FALSE, sp_name) {
  # modify species name suitable for folder/file name
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  # create a matrix with each column capturing grid of covariate values from min to max observed values
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  # create a matrix that repeats grid_length times the mean covariate value
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  post.mean_list <- list() #list with n_model components, each are grid_length x p matrices
  post.mean_rho_list <- list()
  prob.zero_list <- list() #same for probability of zero
  prob.zero_rho_list <- list()
  
  mod_idx <- 1
  subfolder <- paste0("n_",nrow(X))
  for(mod_i in 1:length(mod_list)) {
    mod <- mod_list[[mod_i]]
    mod_rho <- mod_rho_list[[mod_i]]
    
    # save the posterior means
    post.mean_mat <- c() #grid_length x p matrix 
    post.mean_rho_mat <- c()
    prob.zero_mat <- c()
    prob.zero_rho_mat <- c()
    
    # save all the predictions for current model
    all_vars_preds <- list() #list with p components, each n_rep x grid_length
    all_vars_preds_rho <- list()
    
    ### plots for separate curves if ZI model
    ### first model without rho modeled
    if (mod_idx %in% c(2,4)) {
      png(paste0("plots/final_results/",subfolder,"/M",mod_idx,"/separate_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
      par(mfrow = c(3,3),
          mar = c(5,4,2,2))
    }
    
    for (i in 1:ncol(X)) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      if (mod_idx %in% c(1,3)) {
        res <- predict_beta_regression(mod,predX,X,thinning,1,rho_modeled=FALSE)
      } else {
        res <- predict_ZI_beta_regression(mod,predX,X,thinning,1,rho_modeled=FALSE)
      }
      res.EY <- res$EY_sam # n_rep x grid_length
      all_vars_preds[[colnames(X)[i]]] <- res.EY
      
      if (use_median) {
        post.means <- apply(res.EY,2,median)
      } else {
        post.means <- colMeans(res.EY)
      }
      
      #plot the probability of presences and expectation given presence
      if (mod_idx %in% c(2,4)) {
        
        if (use_median) {
          pi <- apply(res$prob_suit_sam,2,median)
          EY_if_suitable <- apply(res$EY_if_suitable_sam,2,median)
        } else {
          pi <- colMeans(res$prob_suit_sam)
          EY_if_suitable <- colMeans(res$EY_if_suitable_sam)        
        }
        
        plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), ylab = "value", xlab = colnames(X)[i], main = "",type="l")
        lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
        if (i == 1){
          legend("topright", legend = c("Pr(suitable)","E[Y|suitable]"), lty = 1, col = c("blue","red"), lwd = 1, cex = 1)
        }
        # add legend
        #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
        #legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
      }
      
      # save posterior means
      post.mean_mat <- cbind(post.mean_mat, post.means)
      prob.zero_mat <- cbind(prob.zero_mat, colMeans(res$probzero_sam))
    }
    
    if (mod_idx %in% c(2,4)) {
      dev.off()
    }
    
    
    ### then the same for log_model
    if (mod_idx %in% c(2,4)) {
      png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/M",mod_idx,"/separate_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
      par(mfrow = c(3,3),
          mar = c(5,4,2,2))
    }
    
    for (i in 1:ncol(X)) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      if (mod_idx %in% c(1,3)) {
        res <- predict_beta_regression(mod_rho,predX,X,thinning,1,rho_modeled=TRUE,1000)
      } else {
        res <- predict_ZI_beta_regression(mod_rho,predX,X,thinning,1,rho_modeled=TRUE,1000)
      }
      res.EY <- res$EY_sam # n_rep x grid_length
      all_vars_preds_rho[[colnames(X)[i]]] <- res.EY
      
      if (use_median) {
        post.means <- apply(res.EY,2,median)
      } else {
        post.means <- colMeans(res.EY)
      }
      #plot the probability of presences and expectation given presence
      if (mod_idx %in% c(2,4)) {
        
        if (use_median) {
          pi <- apply(res$prob_suit_sam,2,median)
          EY_if_suitable <- apply(res$EY_if_suitable_sam,2,median)
        } else {
          pi <- colMeans(res$prob_suit_sam)
          EY_if_suitable <- colMeans(res$EY_if_suitable_sam)          
        }
        
        plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), main = colnames(X)[i],type="l")
        lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
        if (i == 1){
          legend("topright", legend = c("Pr(suitable)","E[Y|suitable]"), lty = 1, col = c("blue","red"), lwd = 1, cex = 1)
        }
        # add legend
        #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
        #legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
      }
      
      # save posterior means
      post.mean_rho_mat <- cbind(post.mean_rho_mat, post.means)
      prob.zero_rho_mat <- cbind(prob.zero_rho_mat, colMeans(res$probzero_sam))
    }
    
    if (mod_idx %in% c(2,4)) {
      dev.off()
    }
    
    
    # now we want to normalize such that largest posterior men gets value 1!
    # rho not modeled
    max_post_mean <- max(post.mean_mat)
    post.mean_mat <- post.mean_mat / max_post_mean #now max value is 1
    
    # rho modeled
    max_post_mean_rho <- max(post.mean_rho_mat)
    post.mean_rho_mat <- post.mean_rho_mat / max_post_mean_rho
    
    # save scaled posterior means
    post.mean_list[[mod_idx]] <- post.mean_mat
    prob.zero_list[[mod_idx]] <- prob.zero_mat
    
    post.mean_rho_list[[mod_idx]] <- post.mean_rho_mat
    prob.zero_rho_list[[mod_idx]] <- prob.zero_rho_mat
    
    # scale also all the other samples
    all_vars_preds <- lapply(all_vars_preds, function(x) (x/max_post_mean))
    all_vars_preds_rho <- lapply(all_vars_preds_rho, function(x) (x/max_post_mean_rho))
    
    
    # draw plots
    # first for common rho model
    png(paste0("plots/final_results/",subfolder,"/M",mod_idx,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    
    par(mfrow = c(3,3),
        mar = c(5,4,2,2))
    
    for (i in 1:ncol(X)) { # index for variable
      plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
           xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "")
      
      n_rep <- nrow(res.EY)
      for (j in 1:n_rep) { # index for sample
        lines(grid_matrix[,i],all_vars_preds[[i]][j,], col = "lightgrey", lwd = 0.8)
      }
      lines(grid_matrix[,i], post.mean_mat[,i], col = "black", lwd = 1)
    }
    
    dev.off()
    
    # then for rho modeled version
    png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/M",mod_idx,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    
    par(mfrow = c(3,3),
        mar = c(5,4,2,2))
    
    for (i in 1:ncol(X)) { # index for variable
      plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
           xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "")
      
      n_rep <- nrow(res.EY)
      for (j in 1:n_rep) { # index for sample
        lines(grid_matrix[,i],all_vars_preds_rho[[i]][j,], col = "lightgrey", lwd = 0.8)
      }
      lines(grid_matrix[,i], post.mean_rho_mat[,i], col = "black", lwd = 1)
    }
    
    dev.off()
    
    mod_idx <- mod_idx + 1
  }
  
  ### then also plot the colmeans on the same plot
  cols <- rainbow(length(post.mean_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  
  col_names <- c("depth","nitrate","oxygen","phosphate","temperature","salinity","current","chlorophyll","light level")
  
  ### first the base models
  png(paste0("plots/final_results/",subfolder,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  for (i in 1:ncol(X)) {
    # plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
    #      xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "", lwd = 1)
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = col_names[i], ylab = "relative exp. coverage", main = "", lwd = 1)
    for (j in 1:length(post.mean_list)) {
      lines(grid_matrix[,i],post.mean_list[[j]][,i], col = cols[j], lwd = 1)
      #lines(grid_matrix[,i],post.mean_rho_list[[j]][,i], col = cols[j], lty = 2)
    }
    if (i == 1){
      #legend("topright", legend = c("base","ZI","RE","ZI+RE","common rho","rho modeled"), lty = c(1,1,1,1,1,2), col = c(cols,"black","black"), lwd = 1, cex = 1)
      legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  dev.off()
  
  ### then the rho_modeled models
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  for (i in 1:ncol(X)) {
    # plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
    #      xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "", lwd = 1)
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = col_names[i], ylab = "relative exp. coverage", main = "", lwd = 1)
    for (j in 1:length(post.mean_list)) {
      lines(grid_matrix[,i],post.mean_rho_list[[j]][,i], col = cols[j], lwd = 1)
    }
    if (i == 1){
      legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  dev.off()
  
  ### same for probability of zero?
  ### first base models
  png(paste0("plots/final_results/",subfolder,"/prob_zero_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  
  cols <- rainbow(length(prob.zero_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = colnames(X)[i], ylab = "prob. of zero", main = "", lwd = 1)
    for (j in 1:length(prob.zero_list)) {
      lines(grid_matrix[,i],prob.zero_list[[j]][,i], col = cols[j], lwd = 1)
      #lines(grid_matrix[,i],prob.zero_rho_list[[j]][,i], col = cols[j], lty = 2) # rho modeled version
    }
    if (i == 1){
      #legend("topright", legend = c("base","ZI","RE","ZI+RE","common rho","rho modeled"), lty = c(1,1,1,1,1,2), col = c(cols,"black","black"), lwd = 1, cex = 1)
      legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  dev.off()
  
  ### then models with rho modeled
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/prob_zero_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  
  cols <- rainbow(length(prob.zero_list))
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = colnames(X)[i], ylab = "prob. of zero", main = "", lwd = 1)
    for (j in 1:length(prob.zero_list)) {
      lines(grid_matrix[,i],prob.zero_rho_list[[j]][,i], col = cols[j], lwd=1) # rho modeled version
    }
    if (i == 1){
      legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  dev.off()
}

plot_and_save_responses_rho_min2 <- function(mod_rho_list,X,grid_length = 200,im_width,im_height, thinning = 20, use_median = FALSE, sp_name) {
  # modify species name suitable for folder/file name
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  # create a matrix with each column capturing grid of covariate values from min to max observed values
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  # create a matrix that repeats grid_length times the mean covariate value
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  post.mean_list <- list() #list with n_model components, each are grid_length x p matrices
  post.mean_rho_list <- list()
  prob.zero_list <- list() #same for probability of zero
  prob.zero_rho_list <- list()
  
  mod_idx <- 1
  subfolder <- paste0("n_",nrow(X))
  for(mod_i in 1:length(mod_rho_list)) {
    #mod <- mod_list[[mod_i]]
    mod_rho <- mod_rho_list[[mod_i]]
    
    # save the posterior means
    #post.mean_mat <- c() #grid_length x p matrix 
    post.mean_rho_mat <- c()
    #prob.zero_mat <- c()
    prob.zero_rho_mat <- c()
    
    # save all the predictions for current model
    #all_vars_preds <- list() #list with p components, each n_rep x grid_length
    all_vars_preds_rho <- list()
    
    ### plots for separate curves if ZI model
    ### first model without rho modeled
    # if (mod_idx %in% c(2,4)) {
    #   png(paste0("plots/final_results/left_right_censored/",subfolder,"/M",mod_idx,"/separate_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    #   par(mfrow = c(3,3),
    #       mar = c(5,4,2,2))
    # }
    
    # for (i in 1:ncol(X)) {
    #   predX <- colmeans_matrix # take all the variables to their mean in the training data
    #   predX[,i] <- grid_matrix[,i] # replace one covariate by grid
    #   predX <- as.data.frame(predX)
    #   if (mod_idx %in% c(1,3)) {
    #     if (mod_idx == 3) {
    #       res <- predict_beta_regression(mod,predX,X,thinning,1,rho_modeled=FALSE,1000,0.5,FALSE) 
    #     } else {
    #       res <- predict_beta_regression(mod,predX,X,thinning,1,rho_modeled=FALSE,1000,0.5,TRUE)
    #     }
    #   } else {
    #     if (mod_idx == 4) {
    #       res <- predict_ZI_beta_regression(mod,predX,X,thinning,1,rho_modeled=FALSE,1000,0.5,FALSE)
    #     } else {
    #       res <- predict_ZI_beta_regression(mod,predX,X,thinning,1,rho_modeled=FALSE,1000,0.5,TRUE)
    #     }
    #   }
    #   res.EY <- res$EY_sam # n_rep x grid_length
    #   all_vars_preds[[colnames(X)[i]]] <- res.EY
    #   
    #   if (use_median) {
    #     post.means <- apply(res.EY,2,median)
    #   } else {
    #     post.means <- colMeans(res.EY)
    #   }
    #   
    #   #plot the probability of presences and expectation given presence
    #   if (mod_idx %in% c(2,4)) {
    #     
    #     if (use_median) {
    #       pi <- apply(res$prob_suit_sam,2,median)
    #       EY_if_suitable <- apply(res$EY_if_suitable_sam,2,median)
    #     } else {
    #       pi <- colMeans(res$prob_suit_sam)
    #       EY_if_suitable <- colMeans(res$EY_if_suitable_sam)        
    #     }
    #     
    #     plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), ylab = "value", xlab = colnames(X)[i], main = "",type="l")
    #     lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
    #     if (i == 1){
    #       legend("topright", legend = c("Pr(suitable)","E[Y|suitable]"), lty = 1, col = c("blue","red"), lwd = 1, cex = 1)
    #     }
    #     # add legend
    #     #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
    #     #legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
    #   }
    #   
    #   # save posterior means
    #   post.mean_mat <- cbind(post.mean_mat, post.means)
    #   prob.zero_mat <- cbind(prob.zero_mat, colMeans(res$probzero_sam))
    # }
    # 
    # if (mod_idx %in% c(2,4)) {
    #   dev.off()
    # }
    
    
    ### then the same for log_model
    if (mod_idx %in% c(2,4)) {
      png(paste0("plots/final_results/scaled_sigmoid/rho_min_2/",subfolder,"/M",mod_idx,"/separate_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
      par(mfrow = c(3,3),
          mar = c(5,4,2,2))
    }
    
    for (i in 1:ncol(X)) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      if (mod_idx %in% c(1,3)) {
        if (mod_idx == 3) {
          res <- predict_beta_regression(mod_rho,predX,X,thinning,1,rho_modeled=TRUE,1000,0.5,FALSE,min_rho = 2) 
        } else {
          res <- predict_beta_regression(mod_rho,predX,X,thinning,1,rho_modeled=TRUE,1000,0.5,TRUE,min_rho = 2)
        }
      } else {
        if (mod_idx == 4) {
          res <- predict_ZI_beta_regression(mod_rho,predX,X,thinning,1,rho_modeled=TRUE,1000,0.5,FALSE,min_rho = 2) 
        } else {
          res <- predict_ZI_beta_regression(mod_rho,predX,X,thinning,1,rho_modeled=TRUE,1000,0.5,TRUE, min_rho = 2)
        }
      }
      res.EY <- res$EY_sam # n_rep x grid_length
      all_vars_preds_rho[[colnames(X)[i]]] <- res.EY
      
      if (use_median) {
        post.means <- apply(res.EY,2,median)
      } else {
        post.means <- colMeans(res.EY)
      }
      #plot the probability of presences and expectation given presence
      if (mod_idx %in% c(2,4)) {
        
        if (use_median) {
          pi <- apply(res$prob_suit_sam,2,median)
          EY_if_suitable <- apply(res$EY_if_suitable_sam,2,median)
        } else {
          pi <- colMeans(res$prob_suit_sam)
          EY_if_suitable <- colMeans(res$EY_if_suitable_sam)          
        }
        
        plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), main = colnames(X)[i],type="l")
        lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
        if (i == 1){
          legend("topright", legend = c("Pr(suitable)","E[Y|suitable]"), lty = 1, col = c("blue","red"), lwd = 1, cex = 1)
        }
        # add legend
        #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
        #legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
      }
      
      # save posterior means
      post.mean_rho_mat <- cbind(post.mean_rho_mat, post.means)
      prob.zero_rho_mat <- cbind(prob.zero_rho_mat, colMeans(res$probzero_sam))
    }
    
    if (mod_idx %in% c(2,4)) {
      dev.off()
    }
    
    
    # now we want to normalize such that largest posterior men gets value 1!
    # rho not modeled
    # max_post_mean <- max(post.mean_mat)
    # post.mean_mat <- post.mean_mat / max_post_mean #now max value is 1
    
    # rho modeled
    max_post_mean_rho <- max(post.mean_rho_mat)
    post.mean_rho_mat <- post.mean_rho_mat / max_post_mean_rho
    
    # save scaled posterior means
    # post.mean_list[[mod_idx]] <- post.mean_mat
    # prob.zero_list[[mod_idx]] <- prob.zero_mat
    
    post.mean_rho_list[[mod_idx]] <- post.mean_rho_mat
    prob.zero_rho_list[[mod_idx]] <- prob.zero_rho_mat
    
    # scale also all the other samples
    #all_vars_preds <- lapply(all_vars_preds, function(x) (x/max_post_mean))
    all_vars_preds_rho <- lapply(all_vars_preds_rho, function(x) (x/max_post_mean_rho))
    
    
    # draw plots
    # first for common rho model
    # png(paste0("plots/final_results/left_right_censored/",subfolder,"/M",mod_idx,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    # 
    # par(mfrow = c(3,3),
    #     mar = c(5,4,2,2))
    # 
    # for (i in 1:ncol(X)) { # index for variable
    #   plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
    #        xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "")
    #   
    #   n_rep <- nrow(res.EY)
    #   for (j in 1:n_rep) { # index for sample
    #     lines(grid_matrix[,i],all_vars_preds[[i]][j,], col = "lightgrey", lwd = 0.8)
    #   }
    #   lines(grid_matrix[,i], post.mean_mat[,i], col = "black", lwd = 1)
    # }
    # 
    # dev.off()
    
    # then for rho modeled version
    png(paste0("plots/final_results/scaled_sigmoid/rho_min_2/",subfolder,"/M",mod_idx,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    
    par(mfrow = c(3,3),
        mar = c(5,4,2,2))
    
    for (i in 1:ncol(X)) { # index for variable
      plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
           xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "")
      
      n_rep <- nrow(res.EY)
      for (j in 1:n_rep) { # index for sample
        lines(grid_matrix[,i],all_vars_preds_rho[[i]][j,], col = "lightgrey", lwd = 0.8)
      }
      lines(grid_matrix[,i], post.mean_rho_mat[,i], col = "black", lwd = 1)
    }
    
    dev.off()
    
    mod_idx <- mod_idx + 1
  }
  
  ### then also plot the colmeans on the same plot
  cols <- rainbow(length(post.mean_rho_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  
  col_names <- c("depth","nitrate","oxygen","phosphate","temperature","salinity","current","chlorophyll","light level")
  
  ### first the base models
  # png(paste0("plots/final_results/left_right_censored/",subfolder,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  # par(mfrow = c(3,3),
  #     mar = c(5,4,2,2))
  # for (i in 1:ncol(X)) {
  #   # plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
  #   #      xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "", lwd = 1)
  #   plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
  #        xlab = col_names[i], ylab = "relative exp. coverage", main = "", lwd = 1)
  #   for (j in 1:length(post.mean_list)) {
  #     lines(grid_matrix[,i],post.mean_list[[j]][,i], col = cols[j], lwd = 1)
  #     #lines(grid_matrix[,i],post.mean_rho_list[[j]][,i], col = cols[j], lty = 2)
  #   }
  #   if (i == 1){
  #     #legend("topright", legend = c("base","ZI","RE","ZI+RE","common rho","rho modeled"), lty = c(1,1,1,1,1,2), col = c(cols,"black","black"), lwd = 1, cex = 1)
  #     legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
  #   }
  # }
  # dev.off()
  
  ### then the rho_modeled models
  png(paste0("plots/final_results/scaled_sigmoid/rho_min_2/",subfolder,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  for (i in 1:ncol(X)) {
    # plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
    #      xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "", lwd = 1)
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = col_names[i], ylab = "relative exp. coverage", main = "", lwd = 1)
    for (j in 1:length(post.mean_rho_list)) {
      lines(grid_matrix[,i],post.mean_rho_list[[j]][,i], col = cols[j], lwd = 1)
    }
    if (i == 1){
      legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  dev.off()
  
  ### same for probability of zero?
  ### first base models
  # png(paste0("plots/final_results/left_right_censored/",subfolder,"/prob_zero_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  # par(mfrow = c(3,3),
  #     mar = c(5,4,2,2))
  # 
  # cols <- rainbow(length(prob.zero_list))
  # #cols <- c("red","royalblue","forestgreen","yellow2")
  # for (i in 1:ncol(X)) {
  #   plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
  #        xlab = colnames(X)[i], ylab = "prob. of zero", main = "", lwd = 1)
  #   for (j in 1:length(prob.zero_list)) {
  #     lines(grid_matrix[,i],prob.zero_list[[j]][,i], col = cols[j], lwd = 1)
  #     #lines(grid_matrix[,i],prob.zero_rho_list[[j]][,i], col = cols[j], lty = 2) # rho modeled version
  #   }
  #   if (i == 1){
  #     #legend("topright", legend = c("base","ZI","RE","ZI+RE","common rho","rho modeled"), lty = c(1,1,1,1,1,2), col = c(cols,"black","black"), lwd = 1, cex = 1)
  #     legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
  #   }
  # }
  # dev.off()
  
  ### then models with rho modeled
  png(paste0("plots/final_results/scaled_sigmoid/rho_min_2/",subfolder,"/prob_zero_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  
  cols <- rainbow(length(prob.zero_list))
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = colnames(X)[i], ylab = "prob. of zero", main = "", lwd = 1)
    for (j in 1:length(prob.zero_rho_list)) {
      lines(grid_matrix[,i],prob.zero_rho_list[[j]][,i], col = cols[j], lwd=1) # rho modeled version
    }
    if (i == 1){
      legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  dev.off()
}

