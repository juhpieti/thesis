### MODEL 2: zero-inflated left-censored beta-regression
### this script
# 1) fits
# 2) plots responses
# 3) plots predicted maps
# 4) ?

# load in packages
library(terra)
library(loo)
library(ggplot2)

# load in stan
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in helper functions
source("codes/helpers.R")

# load in the training data
load("data/estonia_new/train_2020_2021.Rdata")
colnames(df_sub)
colSums(df_sub[,20:38] > 0)

# remove species that are too rare (under 5 observations)
train <- df_sub # more comfortable name
train <- train[,!(colnames(train) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# prepare the covariate matrix
X <- train[,11:19]
X.scaled <- scale_covariates(X)
### add the second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))

# load in the predictive grid
load("data/estonia_new/predictive_grid_1km_all_variables_2021_july.Rdata")
dim(pred_grid_1km_2021_july_df)
colnames(pred_grid_1km_2021_july_df)

predictive_grid <- vect("data/estonia_new/predictive_grid_1km_all_variables_2021_july/predictive_grid_1km_all_variables_2021_july.shp")

# initialize loo table
loo_table <- c()


# loop over the species, save the models
sp_names <- colnames(train)[20:35]
n_chains <- 4
n_iter <- 50
for (sp_name in sp_names[1:1]) {
  y <- train[,sp_name]
  y.01 <- y/100
  
  dat.beta <- list(N = nrow(X.sec_ord),
                   n_var = ncol(X.sec_ord),
                   y = y.01,
                   X = X.sec_ord,
                   a = 1)
  
  mod.ZIbeta <- stan("stan_files/zero_inflated_left_censored_beta_regression.stan",data = dat.beta, chains = n_chains, iter = n_iter, seed = 42,
                           pars = c("mu","prob_suit","logneg_beta_2","logneg_beta_pi_2"), include = FALSE)
  
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  f_name <- paste0("models/M2/",sp_name_modified,".rds")
  
  saveRDS(mod.ZIbeta, f_name)
}


############################ TEST AREA ###################################
### this would make an own script: M3_analyze.R or just analyse_results.R
### which reads fitted models and does things

### map predictions
#pred_grid_2021_july_df$depth <- -1*pred_grid_2021_july_df$depth # note that depth should be positive! as in training data

### species 1: amphi
mod_amphi.ZI <- readRDS("models/M2/Amphibalanus_improvisus.rds")
mod_amphi.ZI <- readRDS("models/demo/M2/amphi.rds")
#mod_chara <- readRDS("models/demo/M1/chara.rds")

start_time <- Sys.time()
pred_list_m2 <- predict_ZI_beta_regression(mod_amphi.ZI,pred_grid_1km_2021_july_df[,2:10],X,10)
pred_list_m2 <- predict_ZI_beta_regression(mod_amphi.ZI,
                                           cbind(pred_grid_1km_2021_july_df[,2:10],pred_grid_1km_2021_july_df$depth / pred_grid_1km_2021_july_df$zsd),
                                           X,10)

#pred_list <- predict_ZI_beta_regression(mod_chara,pred_grid_1km_2021_july_df[,2:10],X,5)
end_time <- Sys.time()
print(end_time - start_time)

locs <- pred_grid_1km_2021_july_df[,c("x","y")]
vals <- colMeans(pred_list_m2$EY_sam)
vars <- apply(pred_list_m2$EY_sam,2,var)
#vals <- colMeans(pred_list$y_sam)
vals.scaled <- vals/max(vals)
test_df <- as.data.frame(cbind(locs, vals, vals.scaled,vars))

### try the cumulative thing
head(test_df)
test_df_ordered <- test_df[order(-test_df$vals), ]
test_df_ordered$hotspot <- 1

# find when 90% of the expected value is reached
tot_sum <- sum(test_df_ordered$vals)
idx <- which(cumsum(test_df_ordered$vals) > 0.7*tot_sum)[1]

# set all outside to be non-hotspot
test_df_ordered[idx:nrow(test_df_ordered),"hotspot"] <- 0


test_vect <- vect(test_df_ordered, geom=c("x","y"), crs = "EPSG:3067")
test_rast <- rast(ext = ext(predictive_grid), res = 1000, crs = "EPSG:3067")
r <- rasterize(test_vect, test_rast, field = "vals")
r.vars <- rasterize(test_vect, test_rast, field = "vars")
r.scaled <- rasterize(test_vect, test_rast, field = "vals.scaled")
r.hotspot <- rasterize(test_vect, test_rast, field = "hotspot")

par(mfrow = c(2,2))
plot(r, colNA = "lightgrey", main = "M2 (w/ ZI,  w/o RE)",plg = list(title = "E[coverage]"),
     xlab = "Easting (m)", ylab = "Northing (m)")
#plot(r.scaled, colNA = "lightgrey", main = "Amphi rel. exp. coverage: M1 (w/o ZI,  w/o RE)")
plot(r.hotspot, colNA = "lightgrey", main = "M2 (w/ ZI, w/o RE)", col = c("red","blue"),
     plg = list(title = "hotspot", legend = c("no","yes")),
     xlab = "Easting (m)", ylab = "Northing (m)")

# plot variance
plot(r.vars,colNA="lightgrey", main = "Amphi variance: M2 (w/ ZI, w/o RE)")



############################################################################################################################

### responses using prediction functions...!
grid_length <- 200
grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
for (i in 1:ncol(X)) {
  x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
  grid_matrix[,i] <- x_grid
}
colnames(grid_matrix) <- colnames(X)
#grid_matrix_scaled <- scale_covariates(X,grid_matrix)

colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
colnames(colmeans_matrix) <- colnames(X)

### go through variables and plot
par(mfrow = c(3,3),
    mar = c(4,4,2,0))


par(mfrow = c(3,4),
    mar = c(4,4,2,0))

mod_list <- list(mod_amphi.ZI, mod_amphi.spat.ZI)

#post.mean_list <- list() #list with n_model components, each are grid_length x p matrices

mod_idx <- 3
for(mod in mod_list) {
  # save the posterior means
  post.mean_mat <- c() #grid_length x p matrix 
  prob.zero_mat <- c()
  
  # save all the predictions for current model
  all_vars_preds <- list() #list with p components, each n_rep x grid_length
  
  par(mfrow = c(3,4),
      mar = c(4,4,2,0))
  
  # make predictions
  for (i in 1:ncol(X)) {
    predX <- colmeans_matrix # take all the variables to their mean in the training data
    predX[,i] <- grid_matrix[,i] # replace one covariate by grid
    predX <- as.data.frame(predX)
    print(i)
    res <- predict_ZI_beta_regression(mod,predX,X,10,1)
    res.EY <- res$EY_sam # n_rep x grid_length
    all_vars_preds[[colnames(X)[i]]] <- res.EY
    post.means <- colMeans(res.EY)
    
    #plot the probability of presences and expectation given presence
    pi <- colMeans(res$prob_suit_sam)
    EY_if_suitable <- colMeans(res$EY_if_suitable_sam)
    plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), main = colnames(X)[i],type="l")
    lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
    if(i == 1) {
      legend("topright",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1)
    }
    
    # save posterior means
    post.mean_mat <- cbind(post.mean_mat, post.means)
    prob.zero_mat <- cbind(prob.zero_mat, colMeans(res$probzero_sam))
  }
  
  # now we want to normalize such that largest posterior men gets value 1!
  max_post_mean <- max(post.mean_mat)
  post.mean_mat <- post.mean_mat / max_post_mean #now max value is 1
  
  # save scaled posterior means
  post.mean_list[[mod_idx]] <- post.mean_mat
  prob.zero_list[[mod_idx]] <- prob.zero_mat
  mod_idx <- mod_idx + 1
  
  # scale also all the other samples
  all_vars_preds <- lapply(all_vars_preds, function(x) (x/max_post_mean))
  
  par(mfrow = c(3,4),
      mar = c(4,4,2,0))
  
  # draw plots
  for (i in 1:ncol(X)) { # index for variable
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i])
    
    n_rep <- nrow(res.EY)
    for (j in 1:n_rep) { # index for sample
      lines(grid_matrix[,i],all_vars_preds[[i]][j,], col = "lightgrey", lwd = 0.8)
    }
    lines(grid_matrix[,i], post.mean_mat[,i], col = "black", lwd = 1)
  }
  
}

par(mfrow = c(3,4),
    mar = c(4,4,2,0))

### then also plot the colmeans on the same plot
cols <- rainbow(length(post.mean_list))
#cols <- c("red","royalblue","forestgreen","yellow2")
for (i in 1:ncol(X)) {
  plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
       xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i], lwd = 1)
  for (j in 1:length(post.mean_list)) {
    lines(grid_matrix[,i],post.mean_list[[j]][,i], col = cols[j], lwd = 1)
  }
}
plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
legend("center", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1.5, bty = "n")

par(mfrow = c(3,4),
    mar = c(4,4,2,0))

### same for probability of zero?
cols <- rainbow(length(prob.zero_list))
#cols <- c("red","royalblue","forestgreen","yellow2")
for (i in 1:ncol(X)) {
  plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
       xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i], lwd = 1)
  for (j in 1:length(prob.zero_list)) {
    lines(grid_matrix[,i],prob.zero_list[[j]][,i], col = cols[j], lwd = 1)
  }
}
plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
legend("center", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1.5, bty = "n")


#############################################################################################

### SEPARATELY FOR PR_SUIT AND EXPECT|SUIT



