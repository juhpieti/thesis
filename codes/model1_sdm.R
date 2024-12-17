### MODEL 1: logistic regression with two settings
### 1) only first order terms
### 2) including second order terms

library(terra)
library(rstan)
library(ggplot2)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in the subset data
#load("data/estonia_sub/estonia_sub_df.RData")
load("data/estonia_new/train_2020_2021.Rdata")
estonia_sub <- df_sub
source("codes/helpers.R")

# in the data there are three species with less than 10 observations, I will delete them
estonia_sub <- estonia_sub[,!(colnames(estonia_sub) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# 1) regression using only first order terms
res_loo <- c()
#coeff_mat <- c()

X <- estonia_sub[,11:19]
### scale the covariates
X.scaled <- scale_covariates(X)


for (sp_name in colnames(estonia_sub)[20:35]) {
  y <- estonia_sub[,sp_name]
  y.bin <- as.numeric(y > 0)
  
  ### fit the model
  mod <- fit_logistic_regression(y.bin,X.scaled,4,2000,FALSE)
  mod.iid <- fit_logistic_regression(y.bin,X.scaled,4,2000,TRUE)
  
  ### check convergence
  check_convergence(mod,FALSE)
  check_convergence(mod.iid,TRUE)
  
  ### calculate loo
  mod.loo <- calc_loo(mod)
  mod.iid.loo <- calc_loo(mod.iid)
  # save to a table
  res_loo <- rbind(res_loo, c(mod.loo,mod.iid.loo))
  
  ### examine the coefficient distributions
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  png(paste0("plots/SDM/M1/log_reg/coefficient_distributions/",sp_name_modified,".png"))
  plot_distributions(mod,X.scaled,FALSE,3,4)
  dev.off()

  png(paste0("plots/SDM/M1/log_reg_noise/coefficient_distributions/",sp_name_modified,".png"))
  plot_distributions(mod.iid,X.scaled,TRUE,4,4)
  dev.off()
  
  ### save the coefficients
  
  
  ### examine the responses
  png(paste0("plots/SDM/M1/log_reg/response_curves/",sp_name_modified,".png"))
  plot_responses(mod,X.scaled,X,TRUE,FALSE,0,50,0,1,3,4)
  dev.off()
  
  png(paste0("plots/SDM/M1/log_reg_noise/response_curves/",sp_name_modified,".png"))
  plot_responses(mod.iid,X.scaled,X,TRUE,FALSE,0,50,0,1,3,4)
  dev.off()
}



# 2) regression using second order terms

### add second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))
#X.sec_ord <- add_interactions(X.sec_ord,"depth","zsd")
res_loo.secord <- c()
coeff_mat_list <- list()

for (sp_name in colnames(estonia_sub)[20:35]) {
  y <- estonia_sub[,sp_name]
  y.bin <- as.numeric(y > 0)
  
  ### fit the model
  mod <- fit_logistic_regression(y.bin,X.sec_ord,4,1000,FALSE)
  mod.iid <- fit_logistic_regression(y.bin,X.sec_ord,4,1000,TRUE)
  
  ### check convergence
  check_convergence(mod,FALSE)
  check_convergence(mod.iid,TRUE)
  
  ### calculate loo
  mod.loo <- calc_loo(mod)
  mod.iid.loo <- calc_loo(mod.iid)
  # save to a table
  res_loo.secord <- rbind(res_loo.secord, c(mod.loo,mod.iid.loo))
  
  ### examine the coefficient distributions
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  png(paste0("plots/SDM/M1/log_reg_sec_ord/new_covariate_distributions/",sp_name_modified,".png"))
  plot_distributions(mod,X.sec_ord,FALSE,5,5)
  dev.off()
  
  png(paste0("plots/SDM/M1/log_reg_sec_ord_noise/new_covariate_distributions/",sp_name_modified,".png"))
  plot_distributions(mod.iid,X.sec_ord,TRUE,5,5)
  dev.off()
  
  ### save the coefficients
  coeff_mat <- c()
  coeff_mat <- cbind(coeff_mat, c(get_posterior_mean(mod, pars = c("alpha","beta"))[,5], s2=0))
  coeff_mat <- cbind(coeff_mat, c(get_posterior_mean(mod.iid, pars = c("alpha","beta","s2"))))
  coeff_mat_list[[sp_name]] <- coeff_mat
  
  ### examine the responses
  png(paste0("plots/SDM/M1/log_reg_sec_ord/new_covariate_responses/",sp_name_modified,".png"))
  plot_responses(mod,X.sec_ord,X,TRUE,TRUE,0,50,0,1,3,4)
  dev.off()
  
  png(paste0("plots/SDM/M1/log_reg_sec_ord_noise/new_covariate_responses/",sp_name_modified,".png"))
  plot_responses(mod.iid,X.sec_ord,X,TRUE,TRUE,0,50,0,1,3,4)
  dev.off()
}



### examine the loo's and save !!!
res_loo
res_loo.secord


res_loo_combined <- cbind(res_loo,res_loo.secord)
rownames(res_loo_combined) <- colnames(estonia_sub)[17:32]
colnames(res_loo_combined) <- c("1st ord","1st ord + noise","2nd ord","2nd ord + noise")

res_loo_combined
apply(res_loo_combined,2,median)

save(res_loo_combined, file = "results/M1/loo_table.Rdata")
