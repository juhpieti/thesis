#####################################################################
### SCRIPT TO PERFORM POSTERIOR PREDICTIVE CHECKS FOR JSDM MODELS ###
#####################################################################

# load in packages
library(terra)

# load in helper functions
source("codes/helpers.R")

# load in the training data
load("data/estonia_new/train/train_2020_2021_n100.Rdata")
train <- train_n100
colnames(train)
colSums(train[,20:38] > 0)

load("data/estonia_new/train/train_2020_2021_all_species_n100.Rdata")
train <- train_n100_all_species
colSums(train[,20:71] > 0)

###################### DATA PREPARATION ##############################

# remove species that are too rare (under 5 observations), and thus not used for modeling
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

pp_check_beta_JSDM <- function(mod,Y_mat,X,sp_name_list,n_rep=50,bin_width = 0.05,test_quantities=FALSE,rho_modeled=FALSE,C=1000,right_censoring=FALSE,a=1,b=0.5,min_rho=0) {
  ### FUNCTION TO PERFORM POSTERIOR PREDICTIVE CHECKS FOR MULTIVARIATE LEFT-CENSORED BETA REGRESSION MODEL ###
  ### INPUTS: ###
  # mod: stanfit object
  # Y_mat: observations (NxJ matrix)
  # X: design matrix to fit the model (Nxp)
  # sp_name_list: list of the names of species modeled
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
  
  X <- as.matrix(X) # from df to matrix to work with matrix multiplication
  n_species <- length(sp_name_list)
  
  # extract posterior draws
  post.samples <- extract(mod)
  
  # make n_rep replications of new data, compare to observed coverages
  n_post_samples <- nrow(post.samples$alpha)
  rep_idx <- sample(1:n_post_samples,n_rep,replace = FALSE) #randomly take n_rep sets of posterior samples
  plot_idx <- sample(rep_idx,15,replace=FALSE) #take 15 to draw the histograms with observed dataset
  
  # initialize list to gather y_rep matrices (n_rep X n) for each species
  y_rep_list <- lapply(sp_name_list, function(x) matrix(0,nrow=n_rep,ncol=nrow(X)))
  names(y_rep_list) <- sp_name_list
  
  # loop over species
  for (j in 1:n_species) {
    # initialize vectors to save test statistics from the replicated data sets
    T_probzero <- c() # proportion of zeros
    T_mean <- c() # mean(y)
    T_mean_positive <- c() #mean(y) for all y > 0
    T_max <- c() # max(y)
    
    ### START PLOTTING
    # 4x4 grid for plotting
    par(mfrow = c(4,4),
        mar = c(2,4,2,0))
    
    # manual breaks
    breaks <- c(-bin_width, 0, seq(bin_width, 1, by=bin_width))
    
    # observations from j:th species
    y <- Y_mat[,j]
    
    # plot the observed data first
    if (!test_quantities) { 
      hist(y, breaks = breaks, xlim = c(0,1), main = "obs", ylim = c(0,length(y)), freq = TRUE)
    }
    
    # go through n_rep posterior samples to replicate n_rep datasets
    for (i in 1:n_rep) {
      idx <- rep_idx[i] # index for posterior sample
      
      # load the model parameters for j:th species and i:th posterior sample
      alpha <- post.samples$alpha[idx,j]
      beta <- c(post.samples$beta_1[idx,,j], post.samples$beta_2[idx,,j])
      lambda <- post.samples$Lambda[idx,,j]
      # latent factors
      Z <- post.samples$Z[idx,,]
      
      # calculate latent f 
      f <- alpha + X %*% beta + Z %*% lambda
      
      # calculate the corresponding mu
      mu <- inv_logit(f)
      
      # prepare rho
      if (rho_modeled) {
        alpha_rho <- post.samples$alpha_rho[idx,j]
        beta_rho <- c(post.samples$beta_rho_1[idx,,j], post.samples$beta_rho_2[idx,,j])
        rho <- min_rho + (C-min_rho)*inv_logit(alpha_rho + X %*% beta_rho)
      } else {
        rho <- rep(post.samples$rho[idx,j],nrow(X))
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
      
      # save also the predicted dataset
      y_rep_list[[j]][i,] <- y_rep
    }
    
    ### draw histograms of test statistics (if test_quantities == TRUE)
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
    }
  }
  return(y_rep_list)
}

plot_community_checks <- function(y_rep_list, Y_mat, plot_richness = FALSE, plot_euclidian_distances = FALSE, test_quantities = FALSE) {
  ### DESCRIPTION OF THE FUNCTION
  # y_rep_list: list (with component for each species) of matrices (n_rep x n) that represent n_rep replicated datasets
  # Y_mat: observed percent covers
  # type: either "richness" or "sum" to tell wether species richness or sum of coverages is plotted
  # test_quantities: logical, if TRUE, summaries using test quantities are plotted, if FALSE, 15 datasets are plotted next to observed dataset
  
  # initialize vectors to save test statistics from the replicated data sets
  T_probzero <- c() # proportion of zeros (no species present)
  T_mean_richness <- c() # mean species richness
  T_mean_sum <- c() # mean sum of covers
  T_mean_shannon <- c() # mean shannon diversity
  
  mean_eucl_dists <- c() # mean euclidian distances ||y_i - ỹ_i||
  
  ### calculate the summaries for observed data sets
  sum_obs <- rowSums(Y_mat)
  sp_richness_obs <- rowSums(Y_mat > 0)
  
  ### calculate diversity metrics for observed data sets
  P_mat_obs <- Y_mat / (sum_obs + 1e-12) # relative percent covers
  
  ### Shannon diversity
  shannon_obs <- rep(0,nrow(Y_mat))
  # loop over species
  n_species <- length(y_rep_list)
  for(j in 1:n_species) {
    shannon_obs <- shannon_obs + P_mat_obs[,j]*log(P_mat_obs[,j] + 1e-12)
  }
  shannon_obs <- -1*shannon_obs
  
  n_rep <- nrow(y_rep_list[[1]]) #number of replications
  plot_idx <- sample(1:n_rep,15,replace=FALSE) #sample ones to show in plots
  
  # set the bin width for histograms
  if (plot_richness) {
    bin_width <- 1
  } else {
    bin_width <- 0.05
  }
  
  if (!test_quantities) {
    par(mfrow = c(4,4),
        mar = c(4,4,2,0))
    # which type of metric used for plotting?
    if (plot_richness) {
      metrics_obs <- sp_richness_obs
    } else {
      metrics_obs <- sum_obs
    }
    
    # histogram with manually placed bins
    if (plot_richness) {
      breaks <- c(-bin_width, 0, seq(bin_width, max(metrics_obs), by=bin_width))
    } else {
      breaks <- c(-bin_width, 0, seq(bin_width, max(metrics_obs)+bin_width, by=bin_width))
    }
    x_label <- ifelse(plot_richness,"sp. richness","total cover")
    #hist(metrics_obs, ylim = c(0,0.5*nrow(Y_mat)), breaks = breaks, main = "obs", xlab = x_label)
    hist(metrics_obs, breaks = breaks, main = "obs", xlab = x_label)
    
  }
  
  # loop over replicated datasets
  for (i in 1:n_rep) {
    ### gather replicated data Y by looping over species
    Y_rep <- c()
    for (j in 1:length(y_rep_list)) {
      Y_rep <- cbind(Y_rep, y_rep_list[[j]][i,])
    }
    ### calculate the summaries for replicated datasets
    sum_rep <- rowSums(Y_rep)
    sp_richness_rep <- rowSums(Y_rep > 0)
    
    ### calculate diversity metrics for observed data sets
    P_mat_rep <- Y_rep / (sum_rep + 1e-12) # relative percent covers
    
    ### Shannon diversity
    shannon_rep <- rep(0,nrow(Y_mat))
    # loop over species
    n_species <- length(y_rep_list)
    for(j in 1:n_species) {
      shannon_rep <- shannon_rep + P_mat_rep[,j]*log(P_mat_rep[,j] + 1e-12)
    }
    shannon_rep <- -1*shannon_rep
    
    ### Calculate also distances ||y_i - ỹ_i||
    eucl_dists <- sqrt(rowSums((Y_mat - Y_rep)^2))
    mean_eucl_dists <- c(mean_eucl_dists, mean(eucl_dists))
    
    if (i %in% plot_idx) {
      if (plot_richness) {
        metrics_rep <- sp_richness_rep
      } else {
        metrics_rep <- sum_rep
      }
      
      if (!test_quantities) {
        # histogram with manually placed bins
        if (plot_richness) {
          breaks <- c(-bin_width, 0, seq(bin_width, max(metrics_rep), by=bin_width)) 
        } else {
          breaks <- c(-bin_width, 0, seq(bin_width, max(metrics_rep)+bin_width, by=bin_width))
        }
        
        x_label <- ifelse(plot_richness,"sp. richness","total cover")
        #hist(metrics_rep, ylim = c(0,0.5*nrow(Y_mat)), breaks = breaks, main = paste0("rep",i), xlab = x_label)
        hist(metrics_rep, breaks = breaks, main = paste0("rep",i), xlab = x_label)
        
      }
    }
    
    # save test quantities
    T_probzero <- c(T_probzero, mean(sp_richness_rep == 0)) # proportion of zeros (no species present)
    T_mean_richness <- c(T_mean_richness, mean(sp_richness_rep)) # mean species richness
    T_mean_sum <- c(T_mean_sum, mean(sum_rep)) # mean sum of covers
    T_mean_shannon <- c(T_mean_shannon, mean(shannon_rep)) # mean shannon diversity
  }
  
  # plot histograms of test quantities (if wanted)
  if (test_quantities) {
    par(mfrow = c(2,2),
        mar = c(2,4,2,0))
    
    # 1) proportion of zeros
    obs_probzero <- mean(sp_richness_obs == 0) # observed proportion of zeros
    probzero_pvalue <- mean(obs_probzero <= T_probzero) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_probzero, breaks = 40, main = "prop. of zero species present",
         xlim = c(min(T_probzero,obs_probzero)-0.1*sd(T_probzero),max(T_probzero,obs_probzero)+0.1*sd(T_probzero)),
         xlab = "T1", probability = TRUE)
    # add legend and vertical line for the observed value
    abline(v = obs_probzero, col = "red", lty = 2, lwd = 2)
    legend("topright",legend = c("observed"),lty=2,lwd=2,col="red")
    legend("topleft",legend=paste0("p-value: ", round(probzero_pvalue,2)), bty = "n")
    
    # 2) mean richness
    obs_richness <- mean(sp_richness_obs) # observed mean richness
    richness_pvalue <- mean(obs_richness <= T_mean_richness) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_mean_richness, breaks = 100, main = "mean sp. richness",
         xlim = c(min(T_mean_richness,obs_richness)-0.1*sd(T_mean_richness),max(T_mean_richness,obs_richness)+0.1*sd(T_mean_richness)),
         xlab = "T2", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_richness, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(richness_pvalue,2)), bty = "n")
    
    # 3) mean sum of covers
    obs_sum <- mean(sum_obs) # observed mean sum of covers
    sum_pvalue <- mean(obs_sum <= T_mean_sum) # bayesian p-value
    
    # draw histogram of n_rep test quantities
    hist(T_mean_sum, breaks = 40, main = "mean total cover",
         xlim = c(min(T_mean_sum,obs_sum)-0.1*sd(T_mean_sum),max(T_mean_sum,obs_sum)+0.1*sd(T_mean_sum)),
         xlab = "T3", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_sum, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(sum_pvalue,2)), bty = "n")
    
    # 4) mean value (for positive observations)
    obs_shannon <- mean(shannon_obs) # observed mean
    shannon_pvalue <- mean(obs_shannon <= T_mean_shannon) # bayesian p-value
    
    # draw histograms of n_rep test quantities
    hist(T_mean_shannon, breaks = 40, main = "mean Shannon diversity",
         xlim = c(min(T_mean_shannon,obs_shannon)-0.1*sd(T_mean_shannon),max(T_mean_shannon,obs_shannon)+0.1*sd(T_mean_shannon)),
         xlab = "T4", probability = TRUE)
    # add legend and vertical line for observed value
    abline(v = obs_shannon, col = "red", lty = 2, lwd = 2)
    legend("topleft",legend=paste0("p-value: ", round(shannon_pvalue,2)), bty = "n")
  }
  
  if(plot_euclidian_distances) {
    par(mfrow = c(1,1))
    hist(mean_eucl_dists, breaks = 40, main = "mean euclidian distances ||y_i - ỹ_i||", xlab = "||y_i - ỹ_i||")
    
    ### draw also a red vertical line for the mean
    avg_mean_eucl_dist <- mean(mean_eucl_dists)
    abline(v = avg_mean_eucl_dist, lty = 2, col = "red", lwd = 2)
    
    # add also the value as text
    text(avg_mean_eucl_dist, 0.9*par("usr")[4], labels = signif(avg_mean_eucl_dist,3), col = "red", cex = 1, pos = 4)
    
  }
}

### test the function
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
Y <- train[,sp_names]

### first with JSDM
set.seed(123)
mod.JSDM_4species <- readRDS("models/multivariate/n_100/9_covariates/JSDM_test.RDS")
Y_rep_list_JSDM_base <- pp_check_beta_JSDM(mod.JSDM_4species,Y/100,X.sec_ord,sp_names,200,0.05,test_quantities = TRUE,rho_modeled = FALSE,1000,FALSE,1,0,0)

set.seed(42)
plot_community_checks(Y_rep_list_JSDM_base,Y/100,plot_richness = FALSE, plot_euclidian_distances = TRUE, test_quantities = TRUE)


### then with stacked SDMs
Y_rep_list_stacked_SDM_base <- list()
set.seed(123)
for(sp_name in sp_names) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/",subfolder,"/M1/",sp_name_modified,".rds"))
  Y_mat <- pp_check_beta(mod, Y[,sp_name]/100,X.sec_ord,200,0.05,test_quantities = TRUE, rho_modeled = FALSE, 1000,FALSE,1,0,0)
  Y_rep_list_stacked_SDM_base[[sp_name]] <- Y_mat
}

set.seed(42)
plot_community_checks(Y_rep_list_stacked_SDM_base,Y/100,plot_richness = TRUE,plot_euclidian_distances = TRUE, test_quantities = TRUE)


### do the same with 21 species?
Y <- train[,20:71]
Y_21 <- Y[,colSums(Y > 0) > 0]
sp_names <- colnames(Y_21)

### pp-checks with JSDM
set.seed(123)
mod.JSDM_21_species <- readRDS("models/multivariate/n_100/M1/JSDM_21species.RDS")
Y_rep_list_JSDM_base <- pp_check_beta_JSDM(mod.JSDM_21_species,Y_21/100,X.sec_ord,sp_names,100,0.05,test_quantities = FALSE,rho_modeled = FALSE,1000,right_censoring = FALSE,1,0,0)

set.seed(42)
plot_community_checks(Y_rep_list_JSDM_base,Y/100,plot_richness = FALSE,plot_euclidian_distances = FALSE,test_quantities = FALSE)


pp_check_ZIbeta_JSDM <- function(mod,y,X,n_rep,bin_width=0.05, test_quantities = FALSE, rho_modeled = FALSE, C = 1000, right_censoring=FALSE, a=1,b=0.5,min_rho=0) {
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

pp_check_beta_spat_JSDM <- function(mod,Y_mat,X,sp_name_list,P,n_rep,bin_width=0.05, test_quantities = FALSE, rho_modeled = FALSE, C=1000, right_censoring = FALSE, a = 1,b = 0.5, min_rho=0) {
  ### FUNCTION TO PERFORM POSTERIOR PREDICTIVE CHECKS FOR LEFT-CENSORED BETA REGRESSION MODEL ###
  ### INPUTS: ###
  # mod: stanfit object
  # Y_mat: observations (NxJ matrix)
  # X: design matrix to fit the model (Nxp)
  # sp_name_list: list of the names of species modeled
  # P: matrix (NxN_grid_cells) that tells in which grid cell each sampling location belongs to (rows sum up to 1)
  # n_rep: how many datasets to create
  # bin_width: width for the bins for histograms
  # test_quantities: TRUE = plot the histograms of 4 test quantities (#zeros, mean, max, mean of positive obs.)
  #                  FALSE = plot the histograms of y
  # rho_modeled: TRUE = rho(x) modeled as function of environmental covariates
  #              FALSE = rho common over locations
  # C: maximum rho value => C used when rho modeled as rho(x) = C*inv_logit(a+xb)) 
  # a: left-censoring constant, latent beta variable is rescaled to (-a,1)
  
  X <- as.matrix(X) # # from df to matrix to work with matrix multiplication
  n_species <- length(sp_name_list)
  
  # extract posterior draws
  post.samples <- extract(mod)
  
  # make n_rep replications of new data, compare to observed coverages
  rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE) #randomly take n_rep sets of posterior samples
  plot_idx <- sample(rep_idx,15,replace=FALSE) #take 15 to draw the histograms with observed dataset
  
  # initialize list to gather y_rep matrices (n_rep X n) for each species
  y_rep_list <- lapply(sp_name_list, function(x) matrix(0,nrow=n_rep,ncol=nrow(X)))
  names(y_rep_list) <- sp_name_list
  
  # loop over species
  for (j in 1:n_species) {
    # initialize vectors to save test statistics from the replicated data sets
    T_probzero <- c() # proportion of zeros
    T_mean <- c() # mean(y)
    T_mean_positive <- c() #mean(y) for all y > 0
    T_max <- c() # max(y)
    
    ### START PLOTTING
    # 4x4 grid for observed dataset + 15 replicated datasets
    par (mfrow = c(4,4),
         mar = c(2,4,2,0))
    
    # manual breaks
    breaks <- c(-bin_width, 0, seq(bin_width, 1, by=bin_width))
    
    # observations for the j:th species
    y <- Y_mat[,j]
    
    # plot the observed data first
    if (!test_quantities) {
      hist(y, breaks = breaks, xlim = c(0,1), main = "obs", ylim = c(0,length(y)), freq = TRUE)
    }
    
    # go through n_rep posterior samples to replicate n_rep datasets
    for (i in 1:n_rep) {
      # take the index of posterior sample
      idx <- rep_idx[i]
      
      # load the model parameters for j:th species and i:th posterior sample
      alpha <- post.samples$alpha[idx,j]
      beta <- c(post.samples$beta_1[idx,,j], post.samples$beta_2[idx,,j])
      lambda <- post.samples$Lambda[idx,,j]
      # latent factors (spatially correlated)
      Z <- post.samples$phi[idx,,]
      
      # calculate latent f 
      f <- alpha + X %*% beta + P %*% Z %*% lambda
      
      # calculate the corresponding mu
      mu <- inv_logit(f)
      
      # prepare rho
      if (rho_modeled) {
        alpha_rho <- post.samples$alpha_rho[idx,j]
        beta_rho <- c(post.samples$beta_rho_1[idx,,j], post.samples$beta_rho_2[idx,,j])
        rho <- min_rho + (C-min_rho)*inv_logit(alpha_rho + X %*% beta_rho)
      } else {
        rho <- rep(post.samples$rho[idx,j],nrow(X))
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
      
      # save also the predicted dataset
      y_rep_list[[j]][i,] <- y_rep
    }
    
    ### draw histograms of test statistics (if test_quantities == TRUE)
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
    }
  }
  return(y_rep_list)
}

### test the function
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
Y <- train[,sp_names]set.seed(123)

### first with JSDM
set.seed(123)
mod.JSDM_spat <- readRDS("models/multivariate/n_100/9_covariates/JSDM_spatial_test.RDS")
Y_rep_list_JSDM_spat <- pp_check_beta_spat_JSDM(mod.JSDM_spat,Y/100,X.sec_ord,sp_names,P,200,0.05,test_quantities = TRUE,rho_modeled = FALSE,1000,FALSE,1,0,0)

set.seed(42)
plot_community_checks(Y_rep_list_JSDM_spat,Y/100,plot_richness = TRUE, plot_euclidian_distances = TRUE, test_quantities = TRUE)


### then with stacked SDMs
Y_rep_list_stacked_SDM_spat <- list()
set.seed(123)
for(sp_name in sp_names) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/",subfolder,"/M3/",sp_name_modified,".rds"))
  Y_mat <- pp_check_beta_spat(mod, Y[,sp_name]/100,X.sec_ord,P,200,0.05,test_quantities = TRUE, rho_modeled = FALSE, 1000,1)
  Y_rep_list_stacked_SDM_spat[[sp_name]] <- Y_mat
}

set.seed(42)
plot_community_checks(Y_rep_list_stacked_SDM_base,Y/100,plot_richness = TRUE,plot_euclidian_distances = TRUE, test_quantities = TRUE)




pp_check_ZIbeta_spat_JSDM <- function(mod,y,X,P,n_rep,bin_width,test_quantities = FALSE, rho_modeled = FALSE, C = 1000, a = 1) {
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