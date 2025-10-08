###########################################################
### SCRIPT TO FIT SPATIAL LEFT-CENSORED BETA REGRESSION ###
###########################################################

# load in packages
library(terra)
library(loo)
library(ggplot2)

# load in stan
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in helper/utility functions
source("codes/helpers.R")

# load in the training data (there are different sizes of training data)
load("data/estonia_new/train/train_2020_2021_n500.Rdata")
df_sub <- train_n500

# load("data/estonia_new/train/train_2020_2021_n1000.Rdata")
# df_sub <- train_n1000
# 
# load("data/estonia_new/train/train_2020_2021_n2000.Rdata")
# df_sub <- train_n2000
# 
# load("data/estonia_new/train/train_2020_2021_n100.Rdata")
# df_sub <- train_n100

# this includes all species
load("data/estonia_new/train/train_2020_2021_all_species_n500.Rdata")
train <- train_n500_all_species

colnames(train)
sort(colSums(train[,20:71] > 0), decreasing = TRUE)
sum(colSums(train[,20:71] > 0) > 0) # 26 observed species
sum(colSums(train[,20:71] > 0) > 2) # 21 observed species

# take only species with at least 3 presences
Y <- train[,20:71]
Y <- Y[,colSums(Y > 0) > 2]

# remove one species to get nicer J = 20
colnames(Y)
Y <- Y[,!(colnames(Y) == "Ranunculus peltatus subsp_ Baudotii")]

# remove species that are too rare (under 5 observations)
# train <- df_sub # more comfortable name
# train <- train[,!(colnames(train) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

### prepare the covariate matrix
X <- train[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd) #approximate light level at the bottom
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

X.scaled <- scale_covariates(X) #scale the covariates
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled)) #add second order terms

### PREPARE THE DATA FOR SPATIAL MODELING
### load in the coarse spatial random effect grid
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
for (i in 1:nrow(estonia_sub)) {
  P[i,as.character(nearest_grid_center.vec[i])] <- 1
}

# turn the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000

# save P and locations for remote computer (which does not have terra)
save(P, observed_grid_cells.df, file = "data/estonia_new/for_remote_computer/P_and_coords_n500.Rdata")

### Prepare for stan
n_chains <- 4
n_iter <- 500

# loop over the species and save the models
#sp_names <- colnames(train)[20:35]
sp_names <- colnames(Y)
subfolder <- paste0("n_",nrow(X))
#model_subfolder <- "" # model where rho is fixed
#model_subfolder <- "scaled_sigmoid/" # model where log(rho) = a + XB w/ second order terms

#for (model_subfolder in c("","scaled_sigmoid/")) {
for (model_subfolder in c("")) {
  
  # location of the stan file depends on whether rho is common (stan_files/MODEL.stan) or modeled with covariates (stan_files/scaled_sigmoid/MODEL.stan)
  stan_file_loc <- paste0("stan_files/",model_subfolder,"left_censored_beta_regression_spatial.stan")
  
  for (sp_name in sp_names) { #possibly select a subset of species
    # select and scale percent cover data
    y <- train[,sp_name]
    y.01 <- y/100
    
    # prepare data for stan
    dat.beta_spat <- list(N = nrow(X.sec_ord),
                          n_var = ncol(X.sec_ord),
                          n_obs_grid = ncol(P),
                          y = y.01,
                          X = X.sec_ord,
                          s = observed_grid_cells.df,
                          a = 1,
                          P = P)
    
    # fit the model
    mod.beta_spat <- stan(stan_file_loc,
                          data = dat.beta_spat, chains = n_chains, iter = n_iter, seed = 42,
                          pars = c("mu","z"), include = FALSE)
    
    sp_name_modified <- gsub(" ","_",sp_name)
    sp_name_modified <- gsub("/","_",sp_name_modified)
    
    # set the location to save the model output
    f_name <- paste0("models/",model_subfolder,subfolder,"/M3/",sp_name_modified,".rds")
    #f_name <- paste0("models/",model_subfolder,subfolder,"/9_covariates/M3/",sp_name_modified,".rds")
    
    saveRDS(mod.beta_spat, f_name)
  }
}


