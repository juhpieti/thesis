#################################################################################
### A script to fit logistic regression with given species y and covariates X ###
#################################################################################

### load in the packages
library(terra)
library(rstan)
library(loo)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

### load in the data
df <- read.csv("data/hanko1/velmu.to.Obama.2023-06-19.csv", header = TRUE)
df <- df[complete.cases(df), ] #drop rows with NAs


select_columns <- function(X,name_list) {
  return(X[,name_list])
}

var_names <- c("sal","surf_expo","turb","seashare","depth_classes","ruov","X_coord","Y_coord")
X_full <- select_columns(X,var_names)

# scale to 0 mean unit variances
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

X_full.scaled <- scale_covariates(X_full)

### add second order terms
add_second_order_terms <- function(X,name_list) {
  ### adds second order terms to data matrix X
  # X: a data matrix (nxp)
  # name_list: vector of column names to add second order terms
  # RETURNS a new data matrix (nx(p+m)) where m is the number of second order terms
  for (name in name_list) {
    col_name <- paste0(name,"_^2")
    X[,col_name] <- X[,name]^2
  }
  return(X)
}


### Start the analysis
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

plot_responses <- function(stan_fit, X, is_iid_noise=FALSE, xmin = -3, xmax = 3, ymin = 0, ymax = 0.5) {
  ### plot the coefficient distributions
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta"))
  
  par(mfrow = c(3,3))
  hist(alpha_sam, main = "alpha", xlab = "value", col = "forestgreen", breaks = 20)
  abline(v=0, col = "red", lty = 2, lwd = 2)
  for (i in 1:ncol(X)) {
    hist(beta_sam[,i], main = colnames(X)[i], xlab = "value", col = "forestgreen", breaks = 20)
    abline(v=0,col="red",lty=2,lwd=2)
  }
  
  ### plot the response curves
  par(mfrow = c(3,3))
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

### predict over the area of interest

build_pred_matrix <- function(X_new.rast, X.orig, res = 500) {
  # coarsen the raster to predict in
  X_new.rast_coarse <- rast(ext(X_new.rast),resolution=c(res,res),crs=crs(X_new.rast))
  X_new.rast_coarse <- resample(X_new.rast, X_new.rast_coarse) #values sampled from higher resolution grid

  X_new <- as.data.frame(X_new.rast_coarse, xy = TRUE)
  colnames(X_new)[1] <- "X_coord"
  colnames(X_new)[2] <- "Y_coord"
  X_new <- X_new[complete.cases(X_new), ]
  
  return(X_new[,colnames(X.orig)]) # return in same order
}

predict_logistic <- function(stan_fit, X_new.df,
                             first_order_names, second_order_names = NULL, is.iid.noise = FALSE,
                             res = 500, thinning = 10) {
  
  X_new.scaled <- X_new.df[,first_order_names] # take the same order as with original data matrix

  ### add the interactions
  if (!is.null(second_order_names)) {
    X_new.scaled <- add_second_order_terms(X_new.scaled, second_order_names)
  }
  
  X_new.scaled <- as.matrix(X_new.scaled)
  
  ### predict
  f_pred <- c()
  alpha.sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta.sam <- as.matrix(stan_fit, pars = c("beta"))
  if (is.iid.noise) {
    s2.sam <- as.matrix(stan_fit, pars = c("s2"))
  }
  
  for(i in seq(1,nrow(beta.sam),thinning)) {
    alpha <- alpha.sam[i]
    beta <- beta.sam[i,]
    if (is.iid.noise) {
      s2 <- s2.sam[i]
    }
    tmp <- as.vector(alpha + X_new.scaled %*% beta)
    if (is.iid.noise) {
      tmp + rnorm(nrow(X_new.scaled),s2)
    }
    f_pred <- rbind(f_pred, tmp)
  }
  return(f_pred)
}

### visualize
plot_prob_occ <- function(f_pred, x_coords, y_coords, coord_syst, res = 500, sp_name = "Zostera.marina") {
  df_pred <- data.frame(prob_occ_PPM = colMeans(inv_logit(f_pred)), x_coord = x_coords, y_coord = y_coords)
  pred.vect <- vect(df_pred, geom = c("x_coord","y_coord"), crs = coord_syst)
  pred.rast <- rast(ext(pred.vect), resolution = c(res,res), crs = coord_syst)
  pred.rast <- rasterize(pred.vect,pred.rast,field = "prob_occ_PPM")
  par(mfrow = c(1,1))
  plot(pred.rast, main = paste0("predicted abundance of ", sp_name))
}



### function that takes in list of variables and list of second order variables, runs the whole thing...?

# responses to 0/1
name_of_species <- "Zostera.marina"
y <- df[,name_of_species]
y.bin <- as.numeric(y > 0) # 1 = present, 0 = absent

### prediction X
# read in the raster layers of covariates
X_new.rast <- rast(c("data/hanko1/sal_rasters_Obama.2023-06-19.tif",
                     "data/hanko1/surf_expo_rasters_Obama.2023-06-19.tif",
                     "data/hanko1/turb_rasters_Obama.2023-06-19.tif",
                     "data/hanko1/seashare_rasters_Obama.2023-06-19.tif",
                     "data/hanko1/depth_classes_rasters_Obama.2023-06-19.tif",
                     "data/hanko1/ruov_rasters_Obama.2023-06-19.tif"))

### make a prediction matrix with certain resolution
X_pred <- build_pred_matrix(X_new.rast, X, 200)
X_pred.scaled <- scale_covariates(X,X_pred)


### FULL MODEL
m1_names <- c("sal","surf_expo","turb","seashare","depth_classes","ruov","X_coord","Y_coord")
X1 <- select_columns(X_full.scaled,m1_names)

### fit the model
m1.mod <- fit_logistic_regression(y.bin,X1,4,500,FALSE)
check_convergence(m1.mod,FALSE)

### calculate LOO-CV
m1.loo <- calc_loo(m1.mod)
m1.loo

plot_responses(m1.mod,X1,FALSE,-5,5,0,1)

m1.f_pred <- predict_logistic(m1.mod, X_pred.scaled, m1_names, NULL, FALSE, 500, 10)
plot_prob_occ(m1.f_pred, X_pred[,"X_coord"], X_pred[,"Y_coord"],crs(X_new.rast),500,"Zostera.marina")


### FULL MODEL WITH SECOND ORDER TERMS
m2_names <- c("sal","surf_expo","turb","seashare","depth_classes","ruov","X_coord","Y_coord")
m2_names_sec_ord <- c("sal", "surf_expo","turb","seashare","depth_classes","ruov")

X2 <- select_columns(X_full.scaled, m2_names)
X2 <- add_second_order_terms(X2, m2_names_sec_ord)

m2.mod <- fit_logistic_regression(y.bin,X2,4,500,FALSE)
check_convergence(m2.mod,FALSE)

m2.loo <- calc_loo(m2.mod)
m2.loo

plot_responses(m2.mod,X2,FALSE,-5,5,0,1)

m2.f_pred <- predict_logistic(m2.mod, X_pred.scaled, m2_names, m2_names_sec_ord, FALSE, 500, 10)
plot_prob_occ(m2.f_pred, X_pred[,"X_coord"], X_pred[,"Y_coord"],crs(X_new.rast),500,"Zostera.marina")

### REDUCED MODEL

