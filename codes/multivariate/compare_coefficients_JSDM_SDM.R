####################################################################################################
### THIS SCRIPT CREATES PLOTS THAT COMPARE THE COEFFICIENTS FROM JDSM MODEL WITH THOSE FROM SDMs ###
####################################################################################################

### MODEL 1 (base)
# load in JSDM
mod.JSDM <- readRDS("models/multivariate/n_100/9_covariates/JSDM_test.RDS")
mod.JSDM.hier.priors <- readRDS("models/multivariate/n_100/9_covariates/JSDM_test_hierarchical_priors.RDS")

### plot the coefficients

# extract posterior sample for JSDM
post.samples <- extract(mod.JSDM)
#post.samples <- extract(mod.JSDM.hier.priors)
beta.samples <- post.samples$beta_1

sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
n_species <- length(sp_names)

cols <- rainbow(n_species)

par(mfrow = c(2,2),
    mar = c(4,4,2,0))
subfolder <- paste0("n_",nrow(X))

#plot(NULL, xlim = c(0,10), ylim = c(-2,2), xlab = "covariate index", ylab = "value", main = "coefficients from base model (SDM vs JSDM)")
for (j in 1:n_species) {
  # calculate means and sds
  means.jsdm <- colMeans(beta.samples[,,j])
  sds.jsdm <- apply(beta.samples[,,j],2,sd)
  # first plot JSDM coefficients for species j
  #points(1:9-0.05, means.jsdm, col = cols[j], pch = 16)
  plot(1:9-0.05, means.jsdm, col = "red", pch = 16, xlab = "covariate index", ylab = "value", main = sp_names[j], xlim = c(0,10), ylim = c(-3,3))
  segments(x0 = 1:9-0.05, y0 = means.jsdm - 2*sds.jsdm, y1 = means.jsdm + 2*sds.jsdm, col = "red") 
  
  # then plot SDM coefficients for the corresponding species
  sp_name <- sp_names[j]
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/",subfolder,"/M1/",sp_name_modified,".rds"))
  samples.sdm <- extract(mod)
  beta.samples.sdm <- samples.sdm$beta_1
  means.sdm <- colMeans(beta.samples.sdm)
  sds.sdm <- apply(beta.samples.sdm, 2, sd)
  points(1:9+0.05, means.sdm, col = "blue", pch = 16)
  segments(x0 = 1:9+0.05, y0 = means.sdm - 2*sds.sdm, y1 = means.sdm + 2*sds.sdm, col = "blue") 
  
  abline(h=0, lty = 2, col = "red")

  
  if (j == 1) {
    legend("bottomright", legend = c("SDM","JSDM"), col = c("blue","red"), pch = c(16,16))
  }
}

# load in SDMs

### MODEL 3 (spatial REs)

mod.JSDM_spat <- readRDS("models/multivariate/n_100/9_covariates/JSDM_spatial_test.RDS")

### plot the coefficients

# extract posterior sample for JSDM
post.samples <- extract(mod.JSDM)
beta.samples <- post.samples$beta_1

sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
n_species <- length(sp_names)

par(mfrow = c(2,2),
    mar = c(4,4,2,0))
#plot(NULL, xlim = c(0,10), ylim = c(-2,2), xlab = "covariate index", ylab = "value", main = "coefficients from base model (SDM vs JSDM)")

subfolder <- paste0("n_",nrow(Y.scaled))

for (j in 1:n_species) {
  # calculate means and sds
  means.jsdm <- colMeans(beta.samples[,,j])
  sds.jsdm <- apply(beta.samples[,,j],2,sd)
  # first plot JSDM coefficients for species j
  #points(1:9-0.05, means.jsdm, col = cols[j], pch = 16)
  plot(1:9-0.05, means.jsdm, col = "red", pch = 16, xlab = "covariate index", ylab = "value", main = sp_names[j], xlim = c(0,10), ylim = c(-3,3))
  segments(x0 = 1:9-0.05, y0 = means.jsdm - 2*sds.jsdm, y1 = means.jsdm + 2*sds.jsdm, col = "red") 
  
  # then plot SDM coefficients for the corresponding species
  sp_name <- sp_names[j]
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/",subfolder,"/M3/",sp_name_modified,".rds"))
  samples.sdm <- extract(mod)
  beta.samples.sdm <- samples.sdm$beta_1
  means.sdm <- colMeans(beta.samples.sdm)
  sds.sdm <- apply(beta.samples.sdm, 2, sd)
  points(1:9+0.05, colMeans(beta.samples.sdm), col = "blue", pch = 16)
  segments(x0 = 1:9+0.05, y0 = means.sdm - 2*sds.sdm, y1 = means.sdm + 2*sds.sdm, col = "blue") 
  
  abline(h=0, lty = 2, col = "red")
  
  
  if (j == 1) {
    legend("bottomright", legend = c("SDM","JSDM"), col = c("blue","red"), pch = c(16,16))
  }
}
