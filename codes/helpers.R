### utility functions to be downloaded in different scripts
library(rdist)

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
  means <- apply(X,2,mean)
  sds <- apply(X,2,sd)
  
  if (!is.null(X_new)) { # if X_new is given, return that instead of X
    X <- X_new
  }
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

fit_logistic_regression <- function(y,X,n_chains=4,n_iter=500,add_iid_noise=FALSE) {
  # y: 0/1 occurrences
  # X: data matrix nxp
  # n_chains: number of chains for MCMC to run
  # n_iter: number of iterations for MCMC
  # add_iid_noise: boolean telling whether to add random noise on top of f = alpha + XB 
  # RETURNS stan model fit
  data_list <- list(N=nrow(X),n_var=ncol(X),y=y,X=X)
  file <- "stan_files/logistic_regression.stan"
  if(add_iid_noise) {
    file <- "stan_files/logistic_regression_iid_noise.stan"
  }
  mod <- stan(file,data=data_list,chains=n_chains,iter=n_iter,seed=42)
}

fit_binomial_regression <- function(y,X,n_chains=4,n_iter=500,add_iid_noise=FALSE) {
  # y: 0/1 occurrences
  # X: data matrix nxp
  # n_chains: number of chains for MCMC to run
  # n_iter: number of iterations for MCMC
  # add_iid_noise: boolean telling whether to add random noise on top of f = alpha + XB 
  # RETURNS stan model fit
  data_list <- list(N=nrow(X),n_var=ncol(X),y=y,X=X)
  file <- "stan_files/binomial_regression.stan"
  if(add_iid_noise) {
    file <- "stan_files/binomial_regression_iid_noise.stan"
  }
  mod <- stan(file,data=data_list,chains=n_chains,iter=n_iter,seed=42)
}
check_convergence <- function(stan_fit,is_iid_noise=FALSE) {
  param_list <- c("alpha","beta")
  if(is_iid_noise) {
    param_list <- c(param_list,"s2")
  }
  print(summary(stan_fit, pars = param_list)$summary)
  stan_trace(stan_fit, pars = param_list)
}

calc_loo <- function(stan_fit) {
  fit.loo <- loo(stan_fit)
  return(fit.loo$estimates[1,1])
}

### plot the inference results
inv_logit <- function(lin.pred) {
  return(1/(1+exp(-lin.pred)))
}

plot_distributions <- function(stan_fit,X,is_iid_noise=FALSE,plot_nx=3,plot_ny=3) {
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta"))
  
  par(mfrow = c(plot_nx,plot_ny))
  hist(alpha_sam, main = "alpha", xlab = "value", col = "forestgreen", breaks = 20)
  abline(v=0, col = "red", lty = 2, lwd = 2)
  
  if (is_iid_noise) {
    s2_sam <- as.matrix(stan_fit, pars = c("s2"))
    hist(s2_sam, main = "s2", xlab = "value", col = "forestgreen", breaks = 20)
  }
  
  for (i in 1:ncol(X)) {
    hist(beta_sam[,i], main = colnames(X)[i], xlab = "value", col = "forestgreen", breaks = 20)
    abline(v=0,col="red",lty=2,lwd=2)
  }
}

plot_responses <- function(stan_fit, X.train, X.orig, is_standardized = TRUE, second_order = FALSE, xmin=-3, xmax=3, ymin=0, ymax=1, plot_nx=3, plot_ny=3) {
  ###
  # xmin: min depth
  # xmax: max depth
  
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta"))
  
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
    #x_grid <- cbind(x_grid, seq(0,100,length=500))
    x_grid <- cbind(x_grid, seq(min(X.orig[,i]),max(X.orig[,i]),length=500))
  }
  #x_grid[,depth_idx] <- seq(xmin,xmax,length=500)
  x_grid_orig <- x_grid
  
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
  
  for (i in 1:ncol(X.train)) {
    var_name <- colnames(X.train)[i]
    f <- mean(alpha_sam) + x_grid[,i] * mean(beta_sam[,i])
    if (second_order) {
      f <- f + x_grid[,i+ncol(X.train)] * mean(beta_sam[,i+ncol(X.train)])
    }
    
    plot(x_grid_orig[,i],inv_logit(f),type="l",
         ylab="prob. of occ.",xlab="standardized value",
         main=colnames(X.train)[i],ylim = c(ymin,ymax))
  }
}

plot_responses_hanko <- function(stan_fit, X.train, X.orig, is_standardized = TRUE, is_iid_noise = FALSE, ymin=0, ymax=1, plot_nx=3, plot_ny=3) {
  ###
  
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta"))
  
  ### response curves
  par(mfrow = c(plot_nx,plot_ny))
  
  ### prepare a grid matrix
  x_grid <- c()
  for (i in 1:ncol(X.train)) {
    x_grid <- cbind(x_grid, seq(min(X.orig[,i]),max(X.orig[,i]),length=500))
  }
  x_grid_orig <- x_grid
  
  ### scale if the X.train also standardized
  if (is_standardized) {
    x_grid <- scale_covariates(X.orig,x_grid)
  }
  
  for (i in 1:ncol(X.train)) {
    var_name <- colnames(X.train)[i]
    f <- mean(alpha_sam) + x_grid[,i] * mean(beta_sam[,i])
    
    plot(x_grid_orig[,i],inv_logit(f),type="l",
         ylab="prob. of occ.",xlab="standardized value",
         main=colnames(X.train)[i],ylim = c(ymin,ymax))
  }
}

plot_responses_categorical <- function(stan_fit, X.train, X.orig=NULL, is_standardized = FALSE, is_iid_noise = FALSE, xmin=0, xmax=60, ymin=0, ymax=1, plot_nx=3, plot_ny = 3) {
  ### plots the response curves
  # stan_fit: stanfit object
  # X.train: matrix X used to fit the model
  # X.orig: original data matrix (in case X.train is a scaled version)
  # is_standardized: boolean telling whether the X.train is standardized (1) or not (0)
  # is_iid_noise: boolean telling whether random noise was added on top of f = a + Xb
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta"))

  par(mfrow = c(plot_nx,plot_ny))
  for(i in which(colnames(X.train) != "depth")) {
    depth_idx <- which(colnames(X.train) == "depth")
    depth_grid <- seq(xmin,xmax,.1)
    depth_grid_orig <- depth_grid
    
    # if standardized regression, depth has to be turn into standardized version
    if (is_standardized) {
      depth_grid <- scale_covariates(as.matrix(X.orig[,depth_idx]),matrix(depth_grid,ncol=1))
    }
    
    var_name <- colnames(X.train)[i]
    f <- mean(alpha_sam) + mean(beta_sam[,i]) + depth_grid * mean(beta_sam[,depth_idx])
    
    plot(depth_grid_orig,inv_logit(f),type="l",
         ylab = "prob. of occ.", xlab = "depth",
         main = colnames(X.train)[i], ylim = c(ymin,ymax))
  }
  
}



### CHECK THIS!

plot_responses_highord <- function(stan_fit, X, is_iid_noise=FALSE, xmin = -3, xmax = 3, ymin = 0, ymax = 0.5, plot_nx=3, plot_ny=3) {
  ### plot the coefficient distributions
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta"))
  
  par(mfrow = c(plot_nx,plot_ny))
  hist(alpha_sam, main = "alpha", xlab = "value", col = "forestgreen", breaks = 20)
  abline(v=0, col = "red", lty = 2, lwd = 2)
  for (i in 1:ncol(X)) {
    hist(beta_sam[,i], main = colnames(X)[i], xlab = "value", col = "forestgreen", breaks = 20)
    abline(v=0,col="red",lty=2,lwd=2)
  }
  
  ### plot the response curves
  par(mfrow = c(plot_nx,plot_ny))
  n_first_ord_terms <- ncol(X) - sum(grepl("_\\^2", colnames(X))) #all columns except the ones that include "_^2" 
  for(i in 1:n_first_ord_terms) {
    x_grid <- seq(xmin,xmax,.1)
    var_name <- colnames(X)[i]
    idx_vec <- grep(var_name, colnames(X)) #take also higher order terms if included
    f <- mean(alpha_sam) + x_grid * mean(beta_sam[,idx_vec[1]]) # add first order term
    if (length(idx_vec) > 1) {
      f <- f + (x_grid^2) * mean(beta_sam[,idx_vec[2]]) # add second order term
    }
    plot(x_grid,inv_logit(f),type="l",
         ylab="prob. of occ.",xlab="standardized value",
         main=colnames(X)[i],ylim = c(ymin,ymax))
  }
}


plot_responses_binomial_regression <- function(stan_fit, X.train, X.orig, is_standardized = TRUE, second_order = FALSE, xmin=-3, xmax=3, ymin=0, ymax=1, plot_nx=3, plot_ny=3) {
  ###
  # xmin: min depth
  # xmax: max depth
  
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta"))
  
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
    x_grid <- cbind(x_grid, seq(0,100,length=500))
  }
  x_grid[,depth_idx] <- seq(xmin,xmax,length=500)
  x_grid_orig <- x_grid
  
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
  
  for (i in 1:ncol(X.train)) {
    var_name <- colnames(X.train)[i]
    f <- mean(alpha_sam) + x_grid[,i] * mean(beta_sam[,i])
    if (second_order) {
      f <- f + x_grid[,i+ncol(X.train)] * mean(beta_sam[,i+ncol(X.train)])
    }
    
    plot(x_grid_orig[,i],100*inv_logit(f),type="l",
         ylab="prob. of occ.",xlab="standardized value",
         main=colnames(X.train)[i],ylim = c(ymin,ymax))
  }
}

plot_responses_beta_regression <- function(stan_fit, X.train, X.orig, is_standardized = TRUE, second_order = TRUE, xmin=-3, xmax=3, ymin=0, ymax=1,
                                           plot_nx=3, plot_ny=3,a=1,type="expectation") {
  ###
  # xmin: min depth
  # xmax: max depth
  # type: either "expectation" or "probzero"
  
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta"))
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



