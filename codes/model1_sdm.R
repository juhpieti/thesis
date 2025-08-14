### LOGISTIC REGRESSION WITH
### 1) only first order terms
### 2) including second order terms

# download packages
library(terra)
library(rstan)
library(ggplot2)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in the subset data
load("data/estonia_new/train_2020_2021.Rdata")
estonia_sub <- df_sub

# run the script including functions
source("codes/helpers.R")

# in the data there are three species with less than 10 observations, delete those
colSums(estonia_sub[,20:38]>0)
estonia_sub <- estonia_sub[,!(colnames(estonia_sub) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# 1) regression using only first order terms
res_loo <- c()

# scale the covariates
X <- estonia_sub[,11:19]
X.scaled <- scale_covariates(X)

### add depth to secchi ratio
X_new <- X
X_new$depth_to_secchi <- X$depth / X$zsd
X.scaled.new <- scale_covariates(X_new)

### add exp(-depth/secchi)
X_new2 <- X
X_new2$light_level <- exp(-X_new2$depth/X_new2$zsd)
X.scaled.new2 <- scale_covariates(X_new2)

for (sp_name in colnames(estonia_sub)[20:35]) {
  y <- estonia_sub[,sp_name]
  y.bin <- as.numeric(y > 0)
  
  ### fit the model
  mod <- fit_logistic_regression(y.bin,X.scaled,4,2000,FALSE)

  ### check convergence
  check_convergence(mod,FALSE)

  ### calculate loo
  mod.loo <- calc_loo(mod)
  # save to a table
  res_loo <- c(res_loo, mod.loo)
  
  ### examine the coefficient distributions
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  png(paste0("plots/SDM/M1/log_reg/coefficient_distributions/",sp_name_modified,".png"))
  plot_distributions(mod,X.scaled,FALSE,3,4)
  dev.off()

  ### examine the responses
  png(paste0("plots/SDM/M1/log_reg/response_curves/",sp_name_modified,".png"))
  plot_responses(mod,X.scaled,X,TRUE,FALSE,0,50,0,1,3,4)
  dev.off()
}

# 2) regression using second order terms

### add second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))
X.sec_ord.new <- add_second_order_terms(X.scaled.new,colnames(X.scaled.new))
X.sec_ord.new2 <- add_second_order_terms(X.scaled.new2,colnames(X.scaled.new2))

res_loo <- c()
for (sp_name in colnames(estonia_sub)[20:24]) {
  y <- estonia_sub[,sp_name]
  y.bin <- as.numeric(y > 0)
  
  ### fit the model (try different set of covariates)
  mod1 <- fit_logistic_regression(y.bin,X.sec_ord,3,300,FALSE)
  mod2 <- fit_logistic_regression(y.bin,X.sec_ord.new,3,300,FALSE)
  mod3 <- fit_logistic_regression(y.bin,X.sec_ord.new2,3,300,FALSE)
  
  ### calculate loo
  mod1.loo <- calc_loo(mod1)
  mod2.loo <- calc_loo(mod2)
  mod3.loo <- calc_loo(mod3)
  
  # save to a table
  res_loo<- rbind(res_loo, c(mod1.loo,mod2.loo,mod3.loo))
  
  ### examine the responses
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  png(paste0("plots/SDM/M1/log_reg_sec_ord/without_extra_covariate_for_depth_secchi/",sp_name_modified,".png"))
  plot_responses(mod1,X.sec_ord,X,TRUE,TRUE,0,50,0,1,3,4)
  dev.off()
  
  png(paste0("plots/SDM/M1/log_reg_sec_ord/with_depth_to_secchi/",sp_name_modified,".png"))
  plot_responses(mod2,X.sec_ord.new,X_new,TRUE,TRUE,0,50,0,1,3,4)
  dev.off()
  
  png(paste0("plots/SDM/M1/log_reg_sec_ord/with_light_level/",sp_name_modified,".png"))
  plot_responses(mod3,X.sec_ord.new2,X_new2,TRUE,TRUE,0,50,0,1,3,4)
  dev.off()
}

### Examine loo-values
res_loo
