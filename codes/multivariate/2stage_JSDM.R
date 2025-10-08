#########################################################
### SCRIPT TO RUN 2-STAGE JSDM FOR PERCENT COVER DATA ###
#########################################################

# load in packages
library(terra)
library(loo)
library(ggplot2)
library(corrplot)

# load in stan
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in utility/helper functions
source("codes/helpers.R")

### load in the dataset(s)

load("data/estonia_new/train/train_2020_2021_all_species_n100.Rdata")
dim(train_n100_all_species)
train <- train_n100_all_species

load("data/estonia_new/train/train_2020_2021_all_species_n250.Rdata")
dim(train_n250_all_species)
train <- train_n250_all_species

load("data/estonia_new/train/train_2020_2021_all_species_n500.Rdata")
dim(train_n500_all_species)
train <- train_n500_all_species

### FULL DATA (n = 3221)
load("data/estonia_new/train/train_2020_2021_full_all_species.Rdata")
dim(df_all_species)
colnames(df_all_species)
sort(colSums(df_all_species[,20:71] > 0), decreasing = TRUE) # still quite many species are not observed at all
sum(colSums(df_all_species[,20:71] > 0) > 0) # so 27/51 species are observed in these 3221 sites
train <- df_all_species

### compare how total mass compares with 4 species to 21 species
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
Y_4 <- train[,sp_names]

colnames(train)
Y <- train[,20:71]
idx_positive <- colSums(Y > 0) > 0 # take only species that have some observations
sum(idx_positive) # it's 21 of those
Y_21 <- Y[,idx_positive]

### see the distributions both on normal scale and log-scale
total_mass_4sp <- rowSums(Y_4)
total_mass_21sp <- rowSums(Y_21)

par(mfrow = c(2,2))
hist(total_mass_4sp, breaks=25)
hist(total_mass_21sp, breaks = 25)
hist(log(total_mass_4sp), breaks = 25)
hist(log(total_mass_21sp), breaks = 25)

### examine how the mass behaves wrt. covariates
X <- train[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd) #approximate the light level at the bottom
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

par(mfrow = c(3,3),
    mar = c(4,4,2,0))
for (i in 1:ncol(X)) {
  plot(X[,i], total_mass_21sp, pch = 18, col = "red", xlab = colnames(X)[i], ylab = "total mass")
  points(X[,i], total_mass_4sp, pch = 18, col = "blue")
  if (i == 1) (legend("topright", legend = c("J=4","J=21"), col = c("blue","red"), pch = c(18,18)))
}

par(mfrow = c(3,3),
    mar = c(4,4,2,0))
for (i in 1:ncol(X)) {
  plot(X[,i], log(total_mass_21sp), pch = 18, col = "red", xlab = colnames(X)[i], ylab = "log(total mass)")
  points(X[,i], log(total_mass_4sp), pch = 18, col = "blue")
  if (i == 1) (legend("topright", legend = c("J=4","J=21"), col = c("blue","red"), pch = c(18,18)))
}

### fit the model
X.scaled <- scale_covariates(X) #scale the covariates
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled)) #add second order terms
Y_4sp.scaled <- Y_4/100

# prepare data for stan
data.list <- list(N = nrow(Y_4),
                  n_var = ncol(X.sec_ord),
                  J = ncol(Y_4),
                  Y = Y_4,
                  y_sum = rowSums(Y_21),
                  X = X.sec_ord)

# stan input parameters
n_chains <- 4
n_iter <- 500

# fit the model
fit.2stage <- stan("stan_files/multivariate/2stage_JSDM.stan",
                      data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                      pars = c("mu"), include = FALSE)

fit.2stage.gamma <- stan("stan_files/multivariate/2stage_JSDM_gamma.stan",
                   data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                   pars = c("mu"), include = FALSE)

fit.2stage.beta <- stan("stan_files/multivariate/2stage_JSDM_scaled_beta.stan",
                         data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                         pars = c("mu"), include = FALSE)

# save the models
subfolder <- paste0("n_",nrow(Y_21),"/")
saveRDS(fit.2stage, paste0("models/multivariate/",subfolder,"/two_stage_model/JSDM_2stage_21species.RDS"))
saveRDS(fit.2stage.gamma, paste0("models/multivariate/",subfolder,"/two_stage_model/JSDM_2stage_21species_gamma.RDS"))


### examine output
fit.2stage <- readRDS("models/multivariate/n_100/two_stage_model/JSDM_2stage_21species.RDS")

# response curves

# post pred checks

par(mfrow = c(4,4),
    mar = c(4,4,2,0))

# plot the observed cover
hist(total_mass_21sp, breaks = 20)
n <- length(total_mass_21sp)

# generate replicates
set.seed(42)
post.samples <- extract(fit.2stage)
n_samples <- nrow(post.samples$beta_1)
rand_idx <- sample(1:n_samples, 15, replace = FALSE)
for (idx in rand_idx) {
  alpha <- post.samples$alpha[idx]
  beta <- c(post.samples$beta_1[idx,], post.samples$beta_2[idx,])
  s <- post.samples$s[idx]
  
  f <- alpha + as.matrix(X.sec_ord) %*% beta
  log_y <- rnorm(n,f,s)
  y <- exp(log_y)
  
  hist(y, breaks = 20)
}

# plot map

### same for gamma model
fit.2stage.gamma <- readRDS("models/multivariate/n_100/two_stage_model/JSDM_2stage_21species_gamma.RDS")

# response curves

# post pred checks

par(mfrow = c(4,4),
    mar = c(4,4,2,0))

bin_width <- 5
breaks <- c(-bin_width,0,seq(bin_width,max(total_mass_21sp)+bin_width, by = bin_width))

# plot the observed cover
#hist(total_mass_21sp, breaks = 20)
hist(total_mass_21sp, breaks = breaks, ylim = c(0,30))
#hist(total_mass_4sp, breaks = 20)
n <- length(total_mass_21sp)

# generate replicates
set.seed(42)
post.samples <- extract(fit.2stage.gamma)
n_samples <- nrow(post.samples$beta_1)
rand_idx <- sample(1:n_samples, 15, replace = FALSE)
for (idx in rand_idx) {
  alpha <- post.samples$alpha[idx]
  beta <- c(post.samples$beta_1[idx,], post.samples$beta_2[idx,])
  s <- post.samples$s[idx]
  
  mu <- exp(alpha + as.matrix(X.sec_ord) %*% beta)
  y <- rgamma(n, mu/s, scale = s)
  #log_y <- rnorm(100,f,s)
  #y <- exp(log_y)
  
  breaks <- c(-bin_width,0,seq(bin_width,max(y)+bin_width, by = bin_width))
  #hist(y, breaks = 20)
  hist(y, breaks = breaks, ylim = c(0,30))
  
}


# plot maps (sum of expectations)

predict_gamma_regression <- function(stan_fit, X.pred, X.orig, thinning = 10) {
  ### function to make predictions with scaled beta regression
  # X.pred: prediction matrix (mxp)
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # thinning: use every thinning:th posterior sample in predictions
  # C: maximum percent cover (sum over species)
  
  # load the posterior samples
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta_1","beta_2"))
  
  s_sam <- as.matrix(stan_fit, pars = c("s"))
  
  # prepare prediction matrix (scale, add second order terms)
  X.pred.scaled <- scale_covariates(X.orig,X.pred)
  Xpred <- add_second_order_terms(X.pred.scaled, colnames(X.pred))
  Xpred <- as.matrix(Xpred) # df to matrix for matrix calculations
  
  # initialize matrices for output
  f_sam <- c()
  y_sam <- c()
  EY_sam <- c()
  s_sample <- c()
  
  # loop over posterior samples (take every thinning:th sample)
  for (i in seq(1,nrow(beta_sam),thinning)) {
    
    # corresponding coefficients
    alpha_i <- alpha_sam[i,]
    beta_i <- beta_sam[i,]
    s_i <- s_sam[i,]
    
    # save the rhos for output
    s_sample <- rbind(s_sample,s_i)
    
    # latent f (linear predictor a + xb)
    f_i <- as.vector(alpha_i + Xpred %*% beta_i)
    f_sam <- rbind(f_sam,f_i)
    
    # mean for beta distribution using logit-link
    mu_i <- exp(f_i)
    
    ### calculate expectations and probability of zeroes
    ### NOTE: mapply iterates over pairs of (mu,rho) with common a
    EYi <- mu_i
    
    # save the expectations and prob. of zeros
    EY_sam <- rbind(EY_sam,EYi)
    
    ### also predict Ys
    # sample latent Vs
    y_i <- rgamma(nrow(Xpred), mu_i/s_i, scale = s_i)
    #y_i <- C*rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
    
    # save the sample
    y_sam <- rbind(y_sam, y_i)
  }
  
  # return everything
  return(list(f_sam = f_sam,
              y_sam = y_sam,
              EY_sam = EY_sam,
              s_sample = s_sample))
}

par(mfrow = c(1,1))
pred_list.gamma <- predict_gamma_regression(fit.2stage.gamma, pred_grid_1km_2021_july_df[,colnames(X)], X, 20)
plot_map(pred_list.gamma$EY_sam,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","expected percent cover",0.7)



######################## SCALED BETA MODEL #########################


fit.2stage.beta <- readRDS("models/multivariate/n_100/two_stage_model/JSDM_2stage_21species_beta.RDS")

# response curves

# post pred checks

par(mfrow = c(4,4),
    mar = c(4,4,2,0))

bin_width <- 5
breaks <- c(-bin_width,0,seq(bin_width,max(total_mass_21sp)+bin_width, by = bin_width))

# bin_width <- 0.05
# breaks <- c(-bin_width,0,seq(bin_width,max(total_mass_21sp/100)+bin_width, by = bin_width))

# plot the observed cover
#hist(total_mass_21sp, breaks = 20)
hist(total_mass_21sp, breaks = breaks, ylim = c(0,0.3*length(total_mass_21sp)))
#hist(total_mass_4sp, breaks = breaks, ylim = c(0,0.3*length(total_mass_4sp)))
#hist(total_mass_21sp/100, breaks = breaks, ylim = c(0,0.3*length(total_mass_21sp)))
n <- length(total_mass_21sp)

C <- 200
#C <- 2

# generate replicates
set.seed(42)
post.samples <- extract(fit.2stage.beta)
n_samples <- nrow(post.samples$beta_1)
rand_idx <- sample(1:n_samples, 15, replace = FALSE)

for (idx in rand_idx) {
  alpha <- post.samples$alpha[idx]
  beta <- c(post.samples$beta_1[idx,], post.samples$beta_2[idx,])
  rho <- post.samples$rho[idx]
  
  mu <- inv_logit(alpha + as.matrix(X.sec_ord) %*% beta)
  y <- C*rbeta(n,mu*rho,(1-mu)*rho)
  #y <- rgamma(100, mu/s, scale = s)
  #log_y <- rnorm(100,f,s)
  #y <- exp(log_y)
  
  breaks <- c(-bin_width,0,seq(bin_width,max(y)+bin_width, by = bin_width))
  #hist(y, breaks = 20, xlim = c(0,200))
  hist(y, breaks = breaks, ylim = c(0,0.3*length(total_mass_21sp)))
}

# plot maps (sum of expectations?)

# load in the predictive grid
predictive_grid <- vect("data/estonia_new/predictive_grid_1km_all_variables_2021_july/predictive_grid_1km_all_variables_2021_july.shp")
load("data/estonia_new/predictive_grid_1km_all_variables_2021_july.Rdata")
dim(pred_grid_1km_2021_july_df)
colnames(pred_grid_1km_2021_july_df)

# add relative light level
pred_grid_1km_2021_july_df$light_bottom <- exp(-1.7*pred_grid_1km_2021_july_df$depth/pred_grid_1km_2021_july_df$zsd)

predict_scaled_beta <- function(stan_fit, X.pred, X.orig, thinning = 10, C = 200) {
  ### function to make predictions with scaled beta regression
  # X.pred: prediction matrix (mxp)
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # thinning: use every thinning:th posterior sample in predictions
  # C: maximum percent cover (sum over species)
  
  # load the posterior samples
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta_1","beta_2"))
  
  rho_sam <- as.matrix(stan_fit, pars = c("rho"))
  
  # prepare prediction matrix (scale, add second order terms)
  X.pred.scaled <- scale_covariates(X.orig,X.pred)
  Xpred <- add_second_order_terms(X.pred.scaled, colnames(X.pred))
  Xpred <- as.matrix(Xpred) # df to matrix for matrix calculations
  
  # initialize matrices for output
  f_sam <- c()
  y_sam <- c()
  EY_sam <- c()
  rho_sample <- c()
  
  # loop over posterior samples (take every thinning:th sample)
  for (i in seq(1,nrow(beta_sam),thinning)) {
    
    # corresponding coefficients
    alpha_i <- alpha_sam[i,]
    beta_i <- beta_sam[i,]
    rho_i <- rho_sam[i,]
    
    # save the rhos for output
    rho_sample <- rbind(rho_sample,rho_i)
    
    # latent f (linear predictor a + xb)
    f_i <- as.vector(alpha_i + Xpred %*% beta_i)
    f_sam <- rbind(f_sam,f_i)
    
    # mean for beta distribution using logit-link
    mu_i <- inv_logit(f_i)
    
    ### calculate expectations and probability of zeroes
    ### NOTE: mapply iterates over pairs of (mu,rho) with common a
    EYi <- mu_i*C

    # save the expectations and prob. of zeros
    EY_sam <- rbind(EY_sam,EYi)

    ### also predict Ys
    # sample latent Vs
    y_i <- C*rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
    
    # save the sample
    y_sam <- rbind(y_sam, y_i)
  }
  
  # return everything
  return(list(f_sam = f_sam,
              y_sam = y_sam,
              EY_sam = EY_sam,
              rho_sample = rho_sample))
}

plot_map <- function(pred_mat, locs, pred.grid.vect, type, title_chr = "", hotspot_proportion = 0.7) {
  ### type: "cover", "hotspot", 
  
  #vals <- colMeans(pred_list$EY_sam)
  vals <- colMeans(pred_mat)
  #vars <- apply(pred_list$EY_sam,2,var)
  vars <- apply(pred_mat,2,var)
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
    # plot(r, colNA = "lightgrey", main = title_chr, ,plg = list(title = "E[coverage]"),
    #      xlab = "Easting (m)", ylab = "Northing (m)")
    plot(r, colNA = "lightgrey", main = title_chr,
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

par(mfrow = c(1,1))
pred_list.scaled_beta <- predict_scaled_beta(fit.2stage.beta, pred_grid_1km_2021_july_df[,colnames(X)], X, 20, 200)
plot_map(pred_list.scaled_beta$EY_sam,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","expected percent cover",0.7)


