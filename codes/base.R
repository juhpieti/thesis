###################################################
### SCRIPT TO FIT LEFT-CENSORED BETA REGRESSION ###
###################################################

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

# # remove species that are too rare (under 5 observations)
# train <- df_sub # more comfortable name
# train <- train[,!(colnames(train) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# prepare the covariate matrix
X <- train[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd) # approximate the light level at bottom
X <- X[,-which(colnames(X) == "zsd")] # remove secchi depth since it is not interesting for modeling in itself

#take a subset of covariates
#cov_names <- c("depth","o2_bottom","bottomT","so_bottom","light_bottom")
#X <- X[,cov_names]

X.scaled <- scale_covariates(X) #scale the covariates
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled)) #add second order terms

# parameters for stan
n_chains <- 4
n_iter <- 500

# loop over the species and save the models
#sp_names <- colnames(train)[20:35]
sp_names <- colnames(Y)
subfolder <- paste0("n_",nrow(X))
#model_subfolder <- "" # model where rho is fixed
#model_subfolder <- "scaled_sigmoid/" # model where log(rho) = a + XB w/ second order terms

#for (model_subfolder in c("","scaled_sigmoid/")) { # for running both models with and without varying rho
for (model_subfolder in c("")) {  
  
  # location of the stan file depends on whether rho is common (stan_files/MODEL.stan) or modeled with covariates (stan_files/scaled_sigmoid/MODEL.stan)
  stan_file_loc <- paste0("stan_files/",model_subfolder,"left_censored_beta_regression.stan")
  
  for (sp_name in sp_names) { # possibly select a subset of species
    # select and scale the percent cover data
    y <- train[,sp_name]
    y.01 <- y/100
    
    # prepare data for stan
    dat.beta <- list(N = nrow(X.sec_ord),
                           n_var = ncol(X.sec_ord),
                           y = y.01,
                           X = X.sec_ord,
                           a = 1)
    
    # run the model
    mod.beta <- stan(stan_file_loc,
                     data = dat.beta, chains = n_chains, iter = n_iter, seed = 42,
                           pars = c("mu"), include = FALSE)
    
    sp_name_modified <- gsub(" ","_",sp_name)
    sp_name_modified <- gsub("/","_",sp_name_modified)
    
    # set the location for saving the model output
    f_name <- paste0("models/",model_subfolder,subfolder,"/M1/",sp_name_modified,".rds")
    #f_name <- paste0("models/",model_subfolder,subfolder,"/9_covariates/M1/",sp_name_modified,".rds")
    
    saveRDS(mod.beta, f_name)
  }
}

# ### LEFT & RIGHT CENSORED
# for (model_subfolder in c("","scaled_sigmoid/")) {
#   
#   stan_file_loc <- paste0("stan_files/",model_subfolder,"left_right_censored_beta_regression.stan")
#   
#   for (sp_name in sp_names[4]) {
#     y <- train[,sp_name]
#     y.01 <- y/100
#     
#     dat.beta <- list(N = nrow(X.sec_ord),
#                      n_var = ncol(X.sec_ord),
#                      y = y.01,
#                      X = X.sec_ord,
#                      a = 1,
#                      b = 0.5)
#     
#     mod.beta <- stan(stan_file_loc,
#                      data = dat.beta, chains = n_chains, iter = n_iter, seed = 42,
#                      pars = c("mu"), include = FALSE)
#     
#     sp_name_modified <- gsub(" ","_",sp_name)
#     sp_name_modified <- gsub("/","_",sp_name_modified)
#     
#     f_name <- paste0("models/left_right_censored/",model_subfolder,subfolder,"/M1/",sp_name_modified,".rds")
#     
#     saveRDS(mod.beta, f_name)
#   }
# }
