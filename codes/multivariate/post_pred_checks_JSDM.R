###########################################################################
### SCRIPT TO PERFORM POSTERIOR PREDICTIVE CHECKS FOR ALL FITTED MODELS ###
###########################################################################

# load in packages
library(terra)

# load in helper functions
source("codes/helpers.R")

# load in the training data
load("data/estonia_new/train/train_2020_2021_n100.Rdata")
df_sub <- train_n100
colnames(df_sub)
colSums(df_sub[,20:38] > 0)

###################### DATA PREPARATION ##############################

# remove species that are too rare (under 5 observations), and thus not used for modeling
train <- df_sub # more comfortable name
train <- train[,!(colnames(train) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# prepare the covariate matrix
X <- train[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd) #approximate amount of light reaching the bottom
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

# scale the covariates
X.scaled <- scale_covariates(X)
# add second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))

### load in spatial grid
spatial_grid <- vect("data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp")

# prepare a matrix of grid centers and their locations (TM35FIN coordinates)
grid_centers <- centroids(spatial_grid)
grid_centers.df <- as.data.frame(grid_centers, geom = "XY")
grid_centers.df <- grid_centers.df[,c("x","y")]
dim(grid_centers.df) # there are m = 191 grid cells

# turn the training data into terra-vector
train.vect <- vect(train, geom = c("x","y"), crs = "EPSG:3067")

# find the spatial grid cell for each sampling location
nearest_grid_center <- nearest(train.vect, grid_centers)
nearest_grid_center.df <- as.data.frame(nearest_grid_center)

# n-length vector indicating the ID of the nearest grid center
nearest_grid_center.vec <- nearest_grid_center.df$to_id

# take the indexes of grid cells that have observations in them
observed_grid_cells <- unique(nearest_grid_center.vec)
observed_grid_cells.df <- grid_centers.df[observed_grid_cells,c("x","y")]

# create matrix P (N x n_grid_cells), where i:th row indicates the spatial grid cell that i:th sampling point is located in
# so each row of matrix P sums up to 1
P <- matrix(0,ncol=length(observed_grid_cells),nrow=nrow(train))
colnames(P) <- rownames(observed_grid_cells.df)
for (i in 1:nrow(train)) {
  P[i,as.character(nearest_grid_center.vec[i])] <- 1
}

### turn the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000

######################### POSTERIOR PREDICTIVE CHECKS ############################

### BASE MODEL (LEFT-CENSORED BETA REGRESSION) ###

pp_check_beta_JSDM <- function(mod,y,X,n_rep=50,bin_width = 0.05,test_quantities=FALSE,rho_modeled=FALSE,C=1000,right_censoring=FALSE,a=1,b=0.5,min_rho=0) {
  ### FUNCTION TO PERFORM POSTERIOR PREDICTIVE CHECKS FOR MULTIVARIATE LEFT-CENSORED BETA REGRESSION MODEL ###
  ### INPUTS: ###
  # mod: stanfit object
  # y: observations (N-vector)
  # X: design matrix to fit the model (Nxp)
  # n_rep: how many datasets to create
  # bin_width: width for the bins for histograms
  # test_quantities: TRUE = plot the histograms of 4 test quantities (#zeros, mean, max, mean of positive obs.)
  #                  FALSE = plot the histograms of y
  # rho_modeled: TRUE = rho(x) modeled as function of environmental covariates
  #              FALSE = rho common over locations
  # C: maximum rho value => C used when rho modeled as rho(x) = C*inv_logit(a+xb)) 
  # right_censoring: TRUE if right-censoring was used to produce ones
  # a: left-censoring constant, latent beta variable is rescaled to (-a,1)
  # b: right-censoring constant, latent beta variable is rescaled to (-a,1+b)
  # min_rho: minimum value of rho, if > 0, then rho(x) = min_rho + (C-min_rho)*inv_logit(a+xb)
  
  # load the posterior samples of model parameters
  alpha.sam <- as.matrix(mod, pars = c("alpha"))
  beta.sam <- as.matrix(mod, pars = c("beta_1","beta_2"))
  
  if (rho_modeled) {
    # parameters related to rho(x)
    alpha_rho_sam <- as.matrix(mod, pars = c("alpha_rho"))
    beta_rho_sam <- as.matrix(mod, pars = c("beta_rho_1","beta_rho_2"))
  } else {
    # posterior sample of common rho
    rho.sam <- as.matrix(mod, pars = c("rho")) 
  }
  
  X <- as.matrix(X) # from df to matrix to work with matrix multiplication
  
  # initialize vectors to save test statistics from the replicated data sets
  T_probzero <- c() # proportion of zeros
  T_mean <- c() # mean(y)
  T_mean_positive <- c() #mean(y) for all y > 0
  T_max <- c() # max(y)
  
  # make n_rep replications of new data, compare to observed coverages
  rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE) #randomly take n_rep sets of posterior samples
  plot_idx <- sample(rep_idx,15,replace=FALSE) #take 15 to draw the histograms with observed dataset
  
  ### START PLOTTING
  # 4x4 grid for plotting
  par(mfrow = c(4,4),
      mar = c(2,4,2,0))
  
  # manual breaks
  breaks <- c(-bin_width, 0, seq(bin_width, 1, by=bin_width))
  
  # plot the observed data first
  if (!test_quantities) { 
    hist(y, breaks = breaks, xlim = c(0,1), main = "obs", ylim = c(0,length(y)), freq = TRUE)
  }
  
  # go through n_rep posterior samples to replicate n_rep datasets
  for (i in 1:n_rep) {
    
    # calculate the corresponding mu
    idx <- rep_idx[i]
    mu <- inv_logit(alpha.sam[idx,] + X %*% beta.sam[idx,])
    
    # calculate (or take) the corresponding rho
    if (rho_modeled) {
      rho <- min_rho + (C-min_rho)*inv_logit(alpha_rho_sam[idx,] + X %*% beta_rho_sam[idx,])
    } else {
      rho <- rep(rho.sam[idx,],nrow(X)) #common rho
    }
    
    # sample latent beta variables
    V <- rbeta(nrow(X),mu*rho,(1-mu)*rho)
    
    # scale & left-censor
    if (!right_censoring) {
      # without right-censoring, V is scaled to (-a,1)
      y_rep <- sapply(V,function(v) (max(0,(a+1)*v - a)))
    } else {
      # if also right-censoring is used, V is scaled to (-a,1+b)
      y_rep <- sapply(V,function(v) (max(0,min(1,(a+b+1)*v - a))))
    }
    
    if (!test_quantities) {
      # plot if the index was chosen to be plotted
      if (idx %in% plot_idx) {
        hist(y_rep, breaks = breaks, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,length(y)), freq = TRUE)
      }
    }
    
    # save the test statistics
    T_probzero <- c(T_probzero, mean(y_rep==0))
    T_mean <- c(T_mean, mean(y_rep))
    T_mean_positive <- c(T_mean_positive, mean(y_rep[y_rep > 0]))
    T_max <- c(T_max, max(y_rep))
  }
  
  ### draw histograms of test statistics
  if (test_quantities) {
    par(mfrow = c(2,2),
        mar = c(2,4,2,0))
    
    # 1) proportion of zeros
    obs_probzero <- mean(y == 0) # observed proportion of zeros
    probzero_pvalue <- mean(obs_probzero <= T_probzero) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_probzero, breaks = 40, main = "proportion of zeros",
         xlim = c(min(T_probzero,obs_probzero)-0.1*sd(T_probzero),max(T_probzero,obs_probzero)+0.1*sd(T_probzero)),
         xlab = "T1", probability = TRUE)
    # add legend and vertical line for the observed value
    abline(v = obs_probzero, col = "red", lty = 2, lwd = 2)
    legend("topright",legend = c("observed"),lty=2,lwd=2,col="red")
    legend("topleft",legend=paste0("p-value: ", round(probzero_pvalue,2)), bty = "n")
    
    # 2) maximum value
    obs_max <- max(y) # observed maximum value
    max_pvalue <- mean(obs_max <= T_max) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_max, breaks = 100, main = "max. value",
         xlim = c(min(T_max,obs_max)-0.1*sd(T_max),max(T_max,obs_max)+0.1*sd(T_max)),
         xlab = "T2", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_max, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(max_pvalue,2)), bty = "n")
    
    # 3) mean value
    obs_mean <- mean(y) # observed mean
    mean_pvalue <- mean(obs_mean <= T_mean) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_mean, breaks = 40, main = "sample mean",
         xlim = c(min(T_mean,obs_mean)-0.1*sd(T_mean),max(T_mean,obs_mean)+0.1*sd(T_mean)),
         xlab = "T3", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_mean, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(mean_pvalue,2)), bty = "n")
    
    # 4) mean value (for positive observations)
    obs_mean_positive <- mean(y[y>0]) # observed mean
    mean_positive_pvalue <- mean(obs_mean_positive <= T_mean_positive) # bayesian p-value
    
    # draw histograms of n_rep test quantities
    hist(T_mean_positive, breaks = 40, main = "sample mean (y>0)",
         xlim = c(min(T_mean_positive,obs_mean_positive)-0.1*sd(T_mean_positive),max(T_mean_positive,obs_mean_positive)+0.1*sd(T_mean_positive)),
         xlab = "T4", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_mean_positive, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(mean_positive_pvalue,2)), bty = "n")
    
    # return the observed bayesian p-values
    return(c(probzero_pval = probzero_pvalue,
             max_pval = max_pvalue,
             mean_pval = mean_pvalue,
             mean_positive_pval = mean_positive_pvalue))
  }
}

### test the function
set.seed(123)
# with common rho
mod.beta <- readRDS(file = "models/n_500/M1/Cladophora_glomerata.rds")
#mod.beta <- readRDS(file = "models/left_right_censored/n_500/M1/Cladophora_glomerata.rds")
pp_check_beta(mod.beta, train[,"Cladophora glomerata"]/100,X.sec_ord,200,0.05,test_quantities=TRUE,rho_modeled=FALSE,1000,FALSE,1,0.5)

# with rho(x) > 0
mod.beta_rho <- readRDS(file = "models/scaled_sigmoid/n_500/M1/Cladophora_glomerata.rds")
#mod.beta_rho <- readRDS(file = "models/left_right_censored/scaled_sigmoid/n_500/M1/Cladophora_glomerata.rds")
pp_check_beta(mod.beta_rho, train[,"Cladophora glomerata"]/100,X.sec_ord,200,0.05,test_quantities=TRUE,rho_modeled=TRUE,1000,FALSE,1,0.5,min_rho=0)

# with rho(x) > 2 (to avoid U-shaped response curves)
mod.beta_rho_min2 <- readRDS(file = "models/scaled_sigmoid/rho_min_2/n_500/M1/Cladophora_glomerata.rds")
pp_check_beta(mod.beta_rho, train[,"Cladophora glomerata"]/100,X.sec_ord,200,0.05,test_quantities=TRUE,rho_modeled=TRUE,1000,FALSE,1,0.5,min_rho=2)


pp_check_ZIbeta <- function(mod,y,X,n_rep,bin_width=0.05, test_quantities = FALSE, rho_modeled = FALSE, C = 1000, right_censoring=FALSE, a=1,b=0.5,min_rho=0) {
  ### FUNCTION TO PERFORM POSTERIOR PREDICTIVE CHECKS FOR ZERO-INFLATED LEFT-CENSORED BETA REGRESSION MODEL ###
  ### INPUTS: ###
  # mod: stanfit object
  # y: observations (N-vector)
  # X: design matrix to fit the model (Nxp)
  # n_rep: how many datasets to create
  # bin_width: width for the bins for histograms
  # test_quantities: TRUE = plot the histograms of 4 test quantities (#zeros, mean, max, mean of positive obs.)
  #                  FALSE = plot the histograms of y
  # rho_modeled: TRUE = rho(x) modeled as function of environmental covariates
  #              FALSE = rho common over locations
  # C: maximum rho value => C used when rho modeled as rho(x) = C*inv_logit(a+xb)) 
  # right_censoring: TRUE if right-censoring was used to produce ones
  # a: left-censoring constant, latent beta variable is rescaled to (-a,1)
  # b: right-censoring constant, latent beta variable is rescaled to (-a,1+b)
  # min_rho: minimum value of rho, if > 0, then rho(x) = min_rho + (C-min_rho)*inv_logit(a+xb)
  
  # load the posterior samples of model parameters
  alpha.sam <- as.matrix(mod, pars = c("alpha"))
  beta.sam <- as.matrix(mod, pars = c("beta_1","beta_2"))
  alpha_pi.sam <- as.matrix(mod, pars = c("alpha_pi"))
  beta_pi.sam <- as.matrix(mod, pars = c("beta_pi_1","beta_pi_2"))
  
  if (rho_modeled) {
    # parameters related to rho(x)
    alpha_rho_sam <- as.matrix(mod, pars = c("alpha_rho"))
    beta_rho_sam <- as.matrix(mod, pars = c("beta_rho_1","beta_rho_2"))
  } else {
    # posterior sample of common rho
    rho.sam <- as.matrix(mod, pars = c("rho")) 
  }
  
  X <- as.matrix(X) # # from df to matrix to work with matrix multiplication
  
  # initialize vectors to save test statistics from the replicated data sets
  T_probzero <- c() # proportion of zeros
  T_mean <- c() # mean(y)
  T_mean_positive <- c() #mean(y) for all y > 0
  T_max <- c() # max(y)
  
  # make n_rep replications of new data, compare to observed coverages
  rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE) #randomly take n_rep sets of posterior samples
  plot_idx <- sample(rep_idx,15,replace=FALSE) #take 15 to draw the histograms with observed dataset
  
  ### START PLOTTING
  # 4x4 grid for plotting
  par(mfrow = c(4,4),
      mar = c(2,4,2,0))
  
  # manual breaks
  breaks <- c(-bin_width, 0, seq(bin_width, 1, by=bin_width))
  
  # plot the observed data first
  if (!test_quantities) {
    hist(y, breaks = breaks, xlim = c(0,1), main = "obs", ylim = c(0,length(y)), freq = TRUE)
  }
  
  # go through n_rep posterior samples to replicate n_rep datasets
  for (i in 1:n_rep) {
    # take the index of posterior sample
    idx <- rep_idx[i]
    
    # calculate the corresponding mu and pi
    mu <- inv_logit(alpha.sam[idx,] + X %*% beta.sam[idx,])
    pi <- inv_logit(alpha_pi.sam[idx,] + X %*% beta_pi.sam[idx,])
    
    # calculate (or take) the corresponding rho
    if (rho_modeled) {
      rho <- min_rho + (C-min_rho)*inv_logit(alpha_rho_sam[idx,] + X %*% beta_rho_sam[idx,])
    } else {
      rho <- rep(rho.sam[idx,],nrow(X)) #common rho
    }
    
    # sample suitabilities from bernoulli(pi)
    Z <- rbinom(nrow(X),1,pi)
    
    # sample latent beta variables
    V <- rbeta(nrow(X),mu*rho,(1-mu)*rho)
    
    # scale & left-censor
    if (!right_censoring) {
      # without right-censoring, V is scaled to (-a,1)
      y_rep <- sapply(V,function(v) (max(0,(a+1)*v - a)))
    } else {
      # if also right-censoring is used, V is scaled to (-a,1+b)
      y_rep <- sapply(V,function(v) (max(0,min(1,(a+b+1)*v - a))))
    }
    
    # if location is not suitable (Z=0), automatically y=0
    y_rep <- Z*y_rep
    
    if (!test_quantities) {
      # plot if the index was chosen to be plotted
      if (idx %in% plot_idx) {
        hist(y_rep, breaks = breaks, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,length(y)), freq = TRUE)
      }
    }
    
    # save the test statistics
    T_probzero <- c(T_probzero, mean(y_rep==0))
    T_mean <- c(T_mean, mean(y_rep))
    T_mean_positive <- c(T_mean_positive, mean(y_rep[y_rep > 0]))
    T_max <- c(T_max, max(y_rep))
  }
  
  ### draw histograms of test statistics
  if (test_quantities) {
    par(mfrow = c(2,2),
        mar = c(2,4,2,0))
    
    # 1) proportion of zeros
    obs_probzero <- mean(y == 0) # observed proportion of zeros
    probzero_pvalue <- mean(obs_probzero <= T_probzero) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_probzero, breaks = 40, main = "proportion of zeros",
         xlim = c(min(T_probzero,obs_probzero)-0.1*sd(T_probzero),max(T_probzero,obs_probzero)+0.1*sd(T_probzero)),
         xlab = "T1", probability = TRUE)
    # add legend and vertical line for the observed value
    abline(v = obs_probzero, col = "red", lty = 2, lwd = 2)
    legend("topright",legend = c("observed"),lty=2,lwd=2,col="red")
    legend("topleft",legend=paste0("p-value: ", round(probzero_pvalue,2)), bty = "n")
    
    # 2) maximum value
    obs_max <- max(y) # observed maximum value
    max_pvalue <- mean(obs_max <= T_max) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_max, breaks = 100, main = "max. value",
         xlim = c(min(T_max,obs_max)-0.1*sd(T_max),max(T_max,obs_max)+0.1*sd(T_max)),
         xlab = "T2", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_max, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(max_pvalue,2)), bty = "n")
    
    # 3) mean value
    obs_mean <- mean(y) # observed mean
    mean_pvalue <- mean(obs_mean <= T_mean) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_mean, breaks = 40, main = "sample mean",
         xlim = c(min(T_mean,obs_mean)-0.1*sd(T_mean),max(T_mean,obs_mean)+0.1*sd(T_mean)),
         xlab = "T3", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_mean, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(mean_pvalue,2)), bty = "n")
    
    # 4) mean value (for positive observations)
    obs_mean_positive <- mean(y[y>0]) # observed mean
    mean_positive_pvalue <- mean(obs_mean_positive <= T_mean_positive) # bayesian p-value
    
    # draw histograms of n_rep test quantities
    hist(T_mean_positive, breaks = 40, main = "sample mean (y>0)",
         xlim = c(min(T_mean_positive,obs_mean_positive)-0.1*sd(T_mean_positive),max(T_mean_positive,obs_mean_positive)+0.1*sd(T_mean_positive)),
         xlab = "T4", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_mean_positive, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(mean_positive_pvalue,2)), bty = "n")
    
    # return the observed bayesian p-values
    return(c(probzero_pval = probzero_pvalue,
             max_pval = max_pvalue,
             mean_pval = mean_pvalue,
             mean_positive_pval = mean_positive_pvalue))
  }
}

# test that the function works
set.seed(123)
# common rho
mod.ZIbeta <- readRDS(file = "models/n_500/M2/Cladophora_glomerata.rds")
#mod.ZIbeta <- readRDS(file = "models/left_right_censored/n_500/M2/Cladophora_glomerata.rds")
pp_check_ZIbeta(mod.ZIbeta, train[,"Cladophora glomerata"]/100, X.sec_ord,100,0.05,test_quantities=TRUE,rho_modeled=FALSE,C=1000,right_censoring=FALSE,a=1,b=0.5)

# rho(x) > 0
mod.ZIbeta_rho <- readRDS(file="models/scaled_sigmoid/n_500/M2/Cladophora_glomerata.rds")
#mod.ZIbeta_rho <- readRDS(file="models/left_right_censored/scaled_sigmoid/n_500/M2/Cladophora_glomerata.rds")
pp_check_ZIbeta(mod.ZIbeta_rho, train[,"Cladophora glomerata"]/100, X.sec_ord,100,0.05,test_quantities=TRUE,rho_modeled=TRUE,C=1000,right_censoring=FALSE,a=1,b=0.5)

pp_check_beta_spat <- function(mod,y,X,P,n_rep,bin_width=0.05, test_quantities = FALSE, rho_modeled = FALSE, C=1000, a = 1) {
  ### FUNCTION TO PERFORM POSTERIOR PREDICTIVE CHECKS FOR LEFT-CENSORED BETA REGRESSION MODEL ###
  ### INPUTS: ###
  # mod: stanfit object
  # y: observations (N-vector)
  # X: design matrix to fit the model (Nxp)
  # P: matrix (NxN_grid_cells) that tells in which grid cell each sampling location belongs to (rows sum up to 1)
  # n_rep: how many datasets to create
  # bin_width: width for the bins for histograms
  # test_quantities: TRUE = plot the histograms of 4 test quantities (#zeros, mean, max, mean of positive obs.)
  #                  FALSE = plot the histograms of y
  # rho_modeled: TRUE = rho(x) modeled as function of environmental covariates
  #              FALSE = rho common over locations
  # C: maximum rho value => C used when rho modeled as rho(x) = C*inv_logit(a+xb)) 
  # a: left-censoring constant, latent beta variable is rescaled to (-a,1)
  
  
  # load the posterior samples of model parameters
  alpha.sam <- as.matrix(mod, pars = c("alpha"))
  beta.sam <- as.matrix(mod, pars = c("beta_1","beta_2"))
  
  if (rho_modeled) {
    # parameters related to rho(x)
    alpha_rho_sam <- as.matrix(mod, pars = c("alpha_rho"))
    beta_rho_sam <- as.matrix(mod, pars = c("beta_rho_1","beta_rho_2"))
  } else {
    # posterior sample of common rho
    rho.sam <- as.matrix(mod, pars = c("rho")) 
  }
  
  # posterior sample of spatial random effects
  phi.sam <- as.matrix(mod, pars = c("phi"))
  
  X <- as.matrix(X) # # from df to matrix to work with matrix multiplication
  
  # initialize vectors to save test statistics from the replicated data sets
  T_probzero <- c() # proportion of zeros
  T_mean <- c() # mean(y)
  T_mean_positive <- c() #mean(y) for all y > 0
  T_max <- c() # max(y)
  
  # make n_rep replications of new data, compare to observed coverages
  rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE) #randomly take n_rep sets of posterior samples
  plot_idx <- sample(rep_idx,15,replace=FALSE) #take 15 to draw the histograms with observed dataset
  
  ### START PLOTTING
  # 4x4 grid for plotting
  par(mfrow = c(4,4),
      mar = c(2,4,2,0))
  
  # manual breaks
  breaks <- c(-bin_width, 0, seq(bin_width, 1, by=bin_width))
  
  # plot the observed data first
  if (!test_quantities) {
    hist(y, breaks = breaks, xlim = c(0,1), main = "obs", ylim = c(0,length(y)), freq = TRUE)
  }
  
  # go through n_rep posterior samples to replicate n_rep datasets
  for (i in 1:n_rep) {
    # take the index of posterior sample
    idx <- rep_idx[i]
    
    # calculate the corresponding mu
    mu <- inv_logit(alpha.sam[idx,] + X %*% beta.sam[idx,] + P %*% phi.sam[idx,])
    
    # calculate (or take) the corresponding rho
    if (rho_modeled) {
      rho <- C*inv_logit(alpha_rho_sam[idx,] + X %*% beta_rho_sam[idx,])
    } else {
      rho <- rep(rho.sam[idx,],nrow(X)) #common rho
    }
    
    # sample latent beta variables
    V <- rbeta(nrow(X),mu*rho,(1-mu)*rho)
    # scale & left-censor
    y_rep <- sapply(V,function(v) (max(0,(a+1)*v - a)))
    
    if(!test_quantities) {
      # plot if the index was chosen to be plotted
      if (idx %in% plot_idx) {
        hist(y_rep, breaks = breaks, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,length(y)), freq = TRUE)
      }
    }
    
    # save the test statistics
    T_probzero <- c(T_probzero, mean(y_rep==0))
    T_mean <- c(T_mean, mean(y_rep))
    T_mean_positive <- c(T_mean_positive, mean(y_rep[y_rep > 0]))
    T_max <- c(T_max, max(y_rep))
  }
  
  ### draw histograms of test statistics
  if (test_quantities) {
    par(mfrow = c(2,2),
        mar = c(2,4,2,0))
    
    # 1) proportion of zeros
    obs_probzero <- mean(y == 0) # observed proportion of zeros
    probzero_pvalue <- mean(obs_probzero <= T_probzero) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_probzero, breaks = 40, main = "proportion of zeros",
         xlim = c(min(T_probzero,obs_probzero)-0.1*sd(T_probzero),max(T_probzero,obs_probzero)+0.1*sd(T_probzero)),
         xlab = "T1", probability = TRUE)
    # add legend and vertical line for the observed value
    abline(v = obs_probzero, col = "red", lty = 2, lwd = 2)
    legend("topright",legend = c("observed"),lty=2,lwd=2,col="red")
    legend("topleft",legend=paste0("p-value: ", round(probzero_pvalue,2)), bty = "n")
    
    # 2) maximum value
    obs_max <- max(y) # observed maximum value
    max_pvalue <- mean(obs_max <= T_max) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_max, breaks = 100, main = "max. value",
         xlim = c(min(T_max,obs_max)-0.1*sd(T_max),max(T_max,obs_max)+0.1*sd(T_max)),
         xlab = "T2", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_max, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(max_pvalue,2)), bty = "n")
    
    # 3) mean value
    obs_mean <- mean(y) # observed mean
    mean_pvalue <- mean(obs_mean <= T_mean) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_mean, breaks = 40, main = "sample mean",
         xlim = c(min(T_mean,obs_mean)-0.1*sd(T_mean),max(T_mean,obs_mean)+0.1*sd(T_mean)),
         xlab = "T3", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_mean, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(mean_pvalue,2)), bty = "n")
    
    # 4) mean value (for positive observations)
    obs_mean_positive <- mean(y[y>0]) # observed mean
    mean_positive_pvalue <- mean(obs_mean_positive <= T_mean_positive) # bayesian p-value
    
    # draw histograms of n_rep test quantities
    hist(T_mean_positive, breaks = 40, main = "sample mean (y>0)",
         xlim = c(min(T_mean_positive,obs_mean_positive)-0.1*sd(T_mean_positive),max(T_mean_positive,obs_mean_positive)+0.1*sd(T_mean_positive)),
         xlab = "T4", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_mean_positive, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(mean_positive_pvalue,2)), bty = "n")
    
    # return the observed bayesian p-values
    return(c(probzero_pval = probzero_pvalue,
             max_pval = max_pvalue,
             mean_pval = mean_pvalue,
             mean_positive_pval = mean_positive_pvalue))
  }
}

# test that the function works
set.seed(123)
# with common rho
mod.beta_spat <- readRDS(file = "models/n_500/M3/Cladophora_glomerata.rds")
pp_check_beta_spat(mod.beta_spat, train[,"Cladophora glomerata"]/100, X.sec_ord,P,100,0.05,test_quantities=TRUE,rho_modeled=FALSE,1000)

# rho(x) > 0
mod.beta_spat_rho <- readRDS(file = "models/scaled_sigmoid/n_500/M3/Cladophora_glomerata.rds")
pp_check_beta_spat(mod.beta_spat_rho, train[,"Cladophora glomerata"]/100, X.sec_ord,P,100,0.05,test_quantities=TRUE,rho_modeled=TRUE,1000)

pp_check_ZIbeta_spat <- function(mod,y,X,P,n_rep,bin_width,test_quantities = FALSE, rho_modeled = FALSE, C = 1000, a = 1) {
  ### FUNCTION TO PERFORM POSTERIOR PREDICTIVE CHECKS FOR LEFT-CENSORED BETA REGRESSION MODEL ###
  ### INPUTS: ###
  # mod: stanfit object
  # y: observations (N-vector)
  # X: design matrix to fit the model (Nxp)
  # P: matrix (NxN_grid_cells) that tells in which grid cell each sampling location belongs to (rows sum up to 1)
  # n_rep: how many datasets to create
  # bin_width: width for the bins for histograms
  # test_quantities: TRUE = plot the histograms of 4 test quantities (#zeros, mean, max, mean of positive obs.)
  #                  FALSE = plot the histograms of y
  # rho_modeled: TRUE = rho(x) modeled as function of environmental covariates
  #              FALSE = rho common over locations
  # C: maximum rho value => C used when rho modeled as rho(x) = C*inv_logit(a+xb)) 
  # a: left-censoring constant, latent beta variable is rescaled to (-a,1)
  
  # load the posterior samples of model parameters
  alpha.sam <- as.matrix(mod, pars = c("alpha"))
  beta.sam <- as.matrix(mod, pars = c("beta_1","beta_2"))
  alpha_pi.sam <- as.matrix(mod, pars = c("alpha_pi"))
  beta_pi.sam <- as.matrix(mod, pars = c("beta_pi_1","beta_pi_2"))
  
  if (rho_modeled) {
    # parameters related to rho(x)
    alpha_rho_sam <- as.matrix(mod, pars = c("alpha_rho"))
    beta_rho_sam <- as.matrix(mod, pars = c("beta_rho_1","beta_rho_2"))
  } else {
    # posterior sample of common rho
    rho.sam <- as.matrix(mod, pars = c("rho")) 
  }
  
  # posterior samples of spatial random effects
  phi_mu.sam <- as.matrix(mod, pars = c("phi_mu"))
  phi_pi.sam <- as.matrix(mod, pars = c("phi_pi"))
  
  X <- as.matrix(X) # # from df to matrix to work with matrix multiplication
  
  # initialize vectors to save test statistics from the replicated data sets
  T_probzero <- c() # proportion of zeros
  T_mean <- c() # mean(y)
  T_mean_positive <- c() #mean(y) for all y > 0
  T_max <- c() # max(y)
  
  # make n_rep replications of new data, compare to observed coverages
  rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE) #randomly take n_rep sets of posterior samples
  plot_idx <- sample(rep_idx,15,replace=FALSE) #take 15 to draw the histograms with observed dataset
  
  ### START PLOTTING
  # 4x4 grid for plotting
  par(mfrow = c(4,4),
      mar = c(2,4,2,0))
  
  # manual breaks
  breaks <- c(-bin_width, 0, seq(bin_width, 1, by=bin_width))
  
  # plot the observed data first
  if (!test_quantities) {
    hist(y, breaks = breaks, xlim = c(0,1), main = "obs", ylim = c(0,length(y)), freq = TRUE)
  }
  
  # go through n_rep posterior samples to replicate n_rep datasets
  for (i in 1:n_rep) {
    # take the index of posterior sample
    idx <- rep_idx[i]
    
    # calculate the corresponding pi,mu
    pi <- inv_logit(alpha_pi.sam[idx,] + X %*% beta_pi.sam[idx,] + P %*% phi_pi.sam[idx,])
    mu <- inv_logit(alpha.sam[idx,] + X %*% beta.sam[idx,] + P %*% phi_mu.sam[idx,])
    
    # calculate (or take) the corresponding rho
    if (rho_modeled) {
      rho <- C*inv_logit(alpha_rho_sam[idx,] + X %*% beta_rho_sam[idx,])
    } else {
      rho <- rep(rho.sam[idx,],nrow(X)) #common rho
    }
    
    # sample the suitability of location
    Z <- rbinom(nrow(X),1,pi)
    # sample latent beta variables
    V <- rbeta(nrow(X),mu*rho,(1-mu)*rho)
    # scale & left-censor
    y_rep <- sapply(V,function(v) (max(0,(a+1)*v - a)))
    
    # if location is not suitable (Z=0), automatically y=0
    y_rep <- Z*y_rep
    
    if (!test_quantities) {
      # plot if the index was chosen to be plotted
      if (idx %in% plot_idx) {
        hist(y_rep, breaks = breaks, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,length(y)), freq = TRUE)
      }
    }
    
    # save test statistics
    T_probzero <- c(T_probzero, mean(y_rep==0))
    T_mean <- c(T_mean, mean(y_rep))
    T_mean_positive <- c(T_mean_positive, mean(y_rep[y_rep > 0]))
    T_max <- c(T_max, max(y_rep))
  }
  
  ### draw histograms of test statistics
  if (test_quantities) {
    par(mfrow = c(2,2),
        mar = c(2,4,2,0))
    
    # 1) proportion of zeros
    obs_probzero <- mean(y == 0) # observed proportion of zeros
    probzero_pvalue <- mean(obs_probzero <= T_probzero) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_probzero, breaks = 40, main = "proportion of zeros",
         xlim = c(min(T_probzero,obs_probzero)-0.1*sd(T_probzero),max(T_probzero,obs_probzero)+0.1*sd(T_probzero)),
         xlab = "T1", probability = TRUE)
    # add legend and vertical line for the observed value
    abline(v = obs_probzero, col = "red", lty = 2, lwd = 2)
    legend("topright",legend = c("observed"),lty=2,lwd=2,col="red")
    legend("topleft",legend=paste0("p-value: ", round(probzero_pvalue,2)), bty = "n")
    
    # 2) maximum value
    obs_max <- max(y) # observed maximum value
    max_pvalue <- mean(obs_max <= T_max) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_max, breaks = 100, main = "max. value",
         xlim = c(min(T_max,obs_max)-0.1*sd(T_max),max(T_max,obs_max)+0.1*sd(T_max)),
         xlab = "T2", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_max, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(max_pvalue,2)), bty = "n")
    
    # 3) mean value
    obs_mean <- mean(y) # observed mean
    mean_pvalue <- mean(obs_mean <= T_mean) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_mean, breaks = 40, main = "sample mean",
         xlim = c(min(T_mean,obs_mean)-0.1*sd(T_mean),max(T_mean,obs_mean)+0.1*sd(T_mean)),
         xlab = "T3", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_mean, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(mean_pvalue,2)), bty = "n")
    
    # 4) mean value (for positive observations)
    obs_mean_positive <- mean(y[y>0]) # observed mean
    mean_positive_pvalue <- mean(obs_mean_positive <= T_mean_positive) # bayesian p-value
    
    # draw histograms of n_rep test quantities
    hist(T_mean_positive, breaks = 40, main = "sample mean (y>0)",
         xlim = c(min(T_mean_positive,obs_mean_positive)-0.1*sd(T_mean_positive),max(T_mean_positive,obs_mean_positive)+0.1*sd(T_mean_positive)),
         xlab = "T4", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_mean_positive, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(mean_positive_pvalue,2)), bty = "n")
    
    # return the observed bayesian p-values
    return(c(probzero_pval = probzero_pvalue,
             max_pval = max_pvalue,
             mean_pval = mean_pvalue,
             mean_positive_pval = mean_positive_pvalue))
  }
}

# test that the function works
set.seed(123)
mod.ZIBeta_spat <- readRDS(file = "models/n_500/M4/Cladophora_glomerata.rds")
pp_check_ZIbeta_spat(mod.ZIBeta_spat, train[,"Cladophora glomerata"]/100, X.sec_ord,P,100,0.05,test_quantities=TRUE,rho_modeled=FALSE,1000)

mod.ZIBeta_spat_rho <- readRDS(file = "models/scaled_sigmoid/n_500/M4/Cladophora_glomerata.rds")
pp_check_ZIbeta_spat(mod.ZIBeta_spat_rho, train[,"Cladophora glomerata"]/100, X.sec_ord,P,100,0.05,test_quantities=TRUE,rho_modeled=TRUE,1000)