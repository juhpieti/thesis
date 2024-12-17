### Script to test two different models with Estonia vegetation data
### 1) model using covariates as continuous variables
### 2) model using categorized bottom vegetation types

library(terra)
library(rstan)
library(ggplot2)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in the subset data
load("data/estonia_sub/estonia_sub_df.RData")
source("codes/helpers.R")

# in the data there are three species with less than 10 observations, I will delete them
estonia_sub <- estonia_sub[,!(colnames(estonia_sub) %in% c("Furcellaria lumbricalis loose form","Tolypella midifica","Chara tomentosa"))]

# 1) regression using variables as continuous
res_loo <- c()

X <- estonia_sub[,6:16]
### scale the covariates
X.scaled <- scale_covariates(X)


for (sp_name in colnames(estonia_sub)[17:33]) {
  y <- estonia_sub[,sp_name]
  y.bin <- as.numeric(y > 0)
  
  ### fit the model
  mod <- fit_logistic_regression(y.bin,X.scaled,4,1000,FALSE)
  #mod.iid <- fit_logistic_regression(y.bin,X.scaled,4,1000,TRUE)
  
  ### check convergence
  check_convergence(mod,FALSE)
  
  ### calculate loo
  mod.loo <- calc_loo(mod)
  # save to a table
  res_loo <- rbind(res_loo, c(mod.loo))
  
  ### examine the coefficient distributions
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  png(paste0("plots/estonia/simple_regression/coefficient_distributions/",sp_name_modified,".png"))
  plot_distributions(mod,X.scaled,FALSE,3,4)
  dev.off()
  #ggsave(filename = paste0("plots/estonia/simple_regression/coefficient_distributions/",gsub(" ","_",sp_name),".png"), plot = dist_plot)
  
  ### examine the responses
  png(paste0("plots/estonia/simple_regression/response_curves/",sp_name_modified,".png"))
  plot_responses(mod,X.scaled,X,TRUE,FALSE,0,50,0,1,3,4)
  dev.off()
  #ggsave(filename = paste0("plots/estonia/simple_regression/response_curves/",gsub(" ","_",sp_name),".png"))
}



# 2) regression using categorized bottom type classes
par(mfrow = c(1,1))

k_vec <- 1:15
tot_wss <- c()

set.seed(13)

for(k in k_vec) {
  km <- kmeans(X.scaled[,-1], k, nstart = 10)
  tot_wss <- c(tot_wss, km$tot.withinss)
}

plot(k_vec, tot_wss, type = "b",ylim=c(0,max(tot_wss)), main = "Selecting k in k-means",
     xlab = "k", ylab = "total within SS")

### k = 4
k_4 <- kmeans(X[,-1],4,nstart=10)
k_4_clust <- k_4$cluster
cols_4 <- rainbow(length(unique(k_4_clust)))
for (k in unique(k_4_clust)) {
  par(mfrow = c(3,4))
  for (i in 2:11) {
    hist(X[k_4_clust == k,i], breaks = 10, col = cols_4[k], main = colnames(X)[i],
         xlab = "value", xlim = c(0,100), probability = TRUE)
  }
}

k_4$centers

### k = 6
k_6 <- kmeans(X[,-1],6,nstart=10)
k_6_clust <- k_6$cluster
cols_6 <- rainbow(length(unique(k_6_clust)))
for (k in unique(k_6_clust)) {
  par(mfrow = c(3,4))
  for (i in 2:11) {
    hist(X[k_6_clust == k,i], breaks = 10, col = cols_6[k], main = colnames(X)[i],
         xlab = "value", xlim = c(0,100), probability = TRUE)
  }
}

k_6$centers

### k = 8
k_8 <- kmeans(X[,-1],8,nstart=10)
k_8_clust <- k_8$cluster
cols_8 <- rainbow(length(unique(k_8_clust)))
for (k in unique(k_8_clust)) {
  par(mfrow = c(3,4))
  for (i in 2:11) {
    hist(X[k_8_clust == k,i], breaks = 10, col = cols_8[k], main = colnames(X)[i],
         xlab = "value", xlim = c(0,100), probability = TRUE)
  }
}

k_8$centers


### I'll stick with k = 6

set.seed(13)
km <- kmeans(X[,-1], 6, nstart = 10)
km$centers
### 1: rock
### 2: gravel/small boulder
### 3: silt
### 4: clay/gravel
### 5: clayplate/big boulder
### 6: finesand

bottom_types <- c("rock","gravel/small boulder","silt","clay/gravel","clayplate/big boulder","fine sand")

### add a column to the data set 

estonia_sub <- cbind(estonia_sub[,1:16],
                     bottom_type = factor(km$cluster, levels = 1:6, labels = bottom_types),
                     estonia_sub[17:ncol(estonia_sub)])

X.categorical <- estonia_sub[,c("depth","bottom_type")]

### scale the depth
mu.depth <- mean(X.categorical[,1])
sd.depth <- sd(X.categorical[,1])
X.categorical.scaled <- X.categorical
X.categorical.scaled[,1] <- (X.categorical[,1] - mu.depth ) / sd.depth

### turn into a dummy matrix
for (b_type in bottom_types) {
  X.categorical.scaled[,b_type] <- as.numeric(X.categorical[,"bottom_type"] == b_type)
}

### delete the bottom type column
X.categorical.scaled <- X.categorical.scaled[,-2]


res_loo.cat <- c()
for (sp_name in colnames(estonia_sub)[18:34]) {
  y <- estonia_sub[,sp_name]
  y.bin <- as.numeric(y > 0)
  
  mod.cat <- fit_logistic_regression(y.bin,X.categorical.scaled,4,1000)
  check_convergence(mod.cat,FALSE)
  mod.cat.loo <- calc_loo(mod.cat)
  res_loo.cat <- rbind(res_loo.cat, c(mod.cat.loo))
  
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  png(paste0("plots/estonia/categorized_regression/coefficient_distributions/",sp_name_modified,".png"))
  plot_distributions(mod.cat,X.categorical.scaled,FALSE,3,3)
  dev.off()
  
  png(paste0("plots/estonia/categorized_regression/response_curves/",sp_name_modified,".png"))
  plot_responses_categorical(mod.cat,X.categorical.scaled,X.categorical,TRUE,FALSE,0,50,0,1,2,3)
  dev.off()
  
}

### examine the loo's
res_loo_combined <- cbind(res_loo,res_loo.cat)
rownames(res_loo_combined) <- colnames(estonia_sub)[17:33]
colnames(res_loo_combined) <- c("regression","categorized")

res_loo_combined
apply(res_loo_combined,2,median)

save(res_loo_combined, file = "results/estonia/loo_table.Rdata")


### visualize
par(mfrow = c(1,1))
hist(rowSums(estonia_sub[,18:34]), breaks=40)
hist(rowSums(estonia_sub[,18:34] > 0))

