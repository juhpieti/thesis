### MODEL 3: left-censored beta-regression
### 2nd order polynomials

library(terra)
library(rstan)
library(ggplot2)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in the subset data
load("data/estonia_sub/estonia_sub_df.RData")
source("codes/helpers.R")

# in the data there are three species with less than 10 observations, I will delete them
estonia_sub <- estonia_sub[,!(colnames(estonia_sub) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# 1) regression using only first order terms
res_loo <- c()

X <- estonia_sub[,6:16]
### scale the covariates
X.scaled <- scale_covariates(X)


for (sp_name in colnames(estonia_sub)[17:32]) {
  y <- estonia_sub[,sp_name]
  
  ### fit the model
  mod <- fit_binomial_regression(y,X.scaled,4,2000,FALSE)
  mod.iid <- fit_binomial_regression(y,X.scaled,4,2000,TRUE)
  
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
  
  png(paste0("plots/SDM/M2/bin_reg/coefficient_distributions/",sp_name_modified,".png"))
  plot_distributions(mod,X.scaled,FALSE,3,4)
  dev.off()
  
  png(paste0("plots/SDM/M2/bin_reg_noise/coefficient_distributions/",sp_name_modified,".png"))
  plot_distributions(mod.iid,X.scaled,TRUE,4,4)
  dev.off()
  
  ### examine the responses
  png(paste0("plots/SDM/M2/bin_reg/response_curves/",sp_name_modified,".png"))
  plot_responses(mod,X.scaled,X,TRUE,FALSE,0,50,0,1,3,4)
  dev.off()
  
  png(paste0("plots/SDM/M2/bin_reg_noise/response_curves/",sp_name_modified,".png"))
  plot_responses(mod.iid,X.scaled,X,TRUE,FALSE,0,50,0,1,3,4)
  dev.off()
}



# 2) regression using second order terms

### add second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))
res_loo.secord <- c()

### load the coefficients from early models

for (sp_name in colnames(estonia_sub)[17:32]) {
  y <- estonia_sub[,sp_name]
  
  ### fit the model
  mod <- fit_binomial_regression(y,X.sec_ord,4,1000,FALSE)
  mod.iid <- fit_binomial_regression(y,X.sec_ord,4,1000,TRUE)
  
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
  
  png(paste0("plots/SDM/M2/bin_reg_sec_ord/coefficient_distributions/",sp_name_modified,".png"))
  plot_distributions(mod,X.sec_ord,FALSE,5,5)
  dev.off()
  
  png(paste0("plots/SDM/M2/bin_reg_sec_ord_noise/coefficient_distributions/",sp_name_modified,".png"))
  plot_distributions(mod.iid,X.sec_ord,TRUE,5,5)
  dev.off()
  
  ### save the coefficients
  coeff_mat <- coeff_mat_list[[sp_name]]
  coeff_mat <- cbind(coeff_mat, c(get_posterior_mean(mod, pars = c("alpha","beta"))[,5], s2=0))
  coeff_mat <- cbind(coeff_mat, c(get_posterior_mean(mod.iid, pars = c("alpha","beta","s2"))))
  coeff_mat_list[[sp_name]] <- coeff_mat
  
  ### examine the responses
  png(paste0("plots/SDM/M2/bin_reg_sec_ord/response_curves/",sp_name_modified,".png"))
  plot_responses(mod,X.sec_ord,X,TRUE,TRUE,0,50,0,1,3,4)
  dev.off()
  
  png(paste0("plots/SDM/M2/bin_reg_sec_ord_noise/response_curves/",sp_name_modified,".png"))
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

###
save(coeff_mat_list, file = "results/coeff_matrix.Rdata")
load("results/coeff_matrix.Rdata")

for (sp_name in names(coeff_mat_list)) {
  colnames(coeff_mat_list[[sp_name]]) <- c("M1","M1+RE","M2","M2+RE")
  rownames(coeff_mat_list[[sp_name]]) <- c("alpha",colnames(X),paste0(colnames(X),"^2"),"s2")
}


#### test area

y <- estonia_sub[,17]
y.01 <- y/100

test_dat <- list(N = nrow(X.sec_ord),
                 n_var = ncol(X.sec_ord),
                 y = y.01,
                 X = X.sec_ord,
                 a = 1)

test_mod <- stan("stan_files/left_censored_beta_regression.stan",data = test_dat, chains = 4, iter = 1000, seed = 42)

