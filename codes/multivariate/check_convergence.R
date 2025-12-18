##################################################################
### SCRIPT TO CHECK THE CONVERGENCE (R-HATS) FOR FITTED MODELS ###
##################################################################

# load in packages
load(rstan)

# load in list of species names (used for modeling)
load("data/estonia_new/train/sp_names.Rdata")

### SDMs (non-spatial)
for(sp_name in sp_names) {
  print(sp_name)
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/n_500/M1/",sp_name_modified,".rds"))
  Rhats <- summary(mod)$summary[,"Rhat"]
  print(sort(Rhats,TRUE)[1:10])
}

### SDMs (spatial)
for(sp_name in sp_names) {
  print(sp_name)
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/n_500/M3/",sp_name_modified,".rds"))
  Rhats <- summary(mod)$summary[,"Rhat"]
  print(sort(Rhats,TRUE)[1:10])
}

### JSDM
### examine the convergence
mod.JSDM <- readRDS("models/multivariate/n_500/M1/JSDM_hier_priors.RDS")
Rhats <- summary(mod.JSDM)$summary[,"Rhat"]
sort(Rhats,decreasing = TRUE)[1:10] # there are many non-converging chains

# examine the trace plots
stan_trace(mod.JSDM, pars = paste0("Lambda[2,",5:20,"]")) # some chains seem to converge towards both a and -a for a>0

### check the rest of the parameters
Rhats <- Rhats[!grepl("^(Z|Lambda)\\[[0-9]+,[0-9]+\\]$", names(Rhats))]
sort(Rhats, TRUE)[1:10]

# instead of Zs and Lambdas, check the convergence of correlations
check_correlation_convergence <- function(stan_fit, sp_names, n_factors) {
  ### checks the correlation between Lambda*Z (random effects)
  
  post.samples.array <- as.array(stan_fit)
  n_species <- length(sp_names)
  n_samples <- dim(post.samples.array)[1]
  n_chains <- dim(post.samples.array)[2]
  R_hats <- c()
  
  par(mfrow = c(5,5),
      mar = c(2,2,2,0))
  for(i in 1:(n_species-1)) {
    for (j in (i+1):n_species) {
      Sigma_ij <- matrix(0,nrow=n_samples,ncol=n_chains)
      Sigma_ii <- matrix(0,nrow=n_samples,ncol=n_chains)
      Sigma_jj <- matrix(0,nrow=n_samples,ncol=n_chains)
      for (k in 1:n_factors) {
        L_ki <- post.samples.array[,,paste0("Lambda[",k,",",i,"]")]
        L_kj <- post.samples.array[,,paste0("Lambda[",k,",",j,"]")]
        Sigma_ij <- Sigma_ij + L_ki*L_kj
        Sigma_ii <- Sigma_ii + L_ki^2
        Sigma_jj <- Sigma_jj + L_kj^2
      }
      # from covariance to correlation
      R_ij <- Sigma_ij / (sqrt(Sigma_ii)*sqrt(Sigma_jj))
      
      # calculate R-hat
      Rh <- Rhat(R_ij)
      R_hats[paste0("R_",i,"_",j)] <- Rh
      
      # plot chains
      matplot(Sigma_ij, type = "l", lty = 1, lwd = 2, col = 1:n_chains, main = paste0("Cor(LR_",i,",LR_",j,")"))
      legend("bottomleft", legend = paste0("Rhat: ", round(Rh,2)), bty = "n", col = "red")
    }
  }
  return(R_hats)
}

### check the convergence in terms of covariance matrix S = L^T %*% L
Rhats <- check_correlation_convergence(mod.JSDM, sp_names, 2)
sort(Rhats,decreasing = TRUE)[1:10]

### JSDM (spatial)
mod.JSDM.spat <- readRDS("models/multivariate/n_500/M3/JSDM_spatial_hier_priors.RDS")

# check the convergence of correlations
Rhats <- check_correlation_convergence(mod.JSDM.spat, sp_names, 2)
sort(Rhats, decreasing = TRUE)[1:10]

# check the other parameters
Rhats <- summary(mod.JSDM.spat)$summary[,"Rhat"]
Rhats <- Rhats[!grepl("^(phi|Lambda|Z)\\[[0-9]+,[0-9]+\\]$", names(Rhats))]
sort(Rhats, TRUE)[1:10]

### 2-stage (non-spatial)
mod.2stage <- readRDS("models/two_stage/n_500/M1/JSDM_NegBin_DirMult_hier_priors.RDS")

# check the convergence of correlations
Rhats <- check_correlation_convergence(mod.2stage, sp_names, 2)
sort(Rhats, decreasing = TRUE)[1:10]

# check the other parameters
Rhats <- summary(mod.2stage)$summary[,"Rhat"]
Rhats <- Rhats[!grepl("^(Lambda|Z)\\[[0-9]+,[0-9]+\\]$", names(Rhats))]
sort(Rhats, TRUE)[1:10]

### 2-stage (spatial)
mod.2stage.spat <- readRDS("models/two_stage/n_500/M3/JSDM_NegBin_DirMult_spatial_hier_priors.RDS")

# check the convergence of correlations
Rhats <- check_correlation_convergence(mod.2stage.spat, sp_names, 2)
sort(Rhats, decreasing = TRUE)[1:10]

# check the other parameters
Rhats <- summary(mod.2stage.spat)$summary[,"Rhat"]
Rhats <- Rhats[!grepl("^(Lambda|Z_spat|phi|Z)\\[[0-9]+,[0-9]+\\]$", names(Rhats))]
sort(Rhats, TRUE)[1:10]

### CHECK THE CHAIN PLOTS FOR BETAs WITH LARGER R-HAT
stan_trace(mod.2stage.spat, pars = paste0("beta_2[6,",1:20,"]")) #does not look that bad actually, probably would need a bit more iterations since the late chains look fine


mod.2stage.spat.nonhier <- readRDS("models/two_stage/n_500/M3/JSDM_NegBin_DirMult_spatial.RDS")
Rhats <- summary(mod.2stage.spat.nonhier)$summary[,"Rhat"]
Rhats <- Rhats[!grepl("^(Lambda|Z_spat|phi|Z)\\[[0-9]+,[0-9]+\\]$", names(Rhats))]
sort(Rhats, TRUE)[1:10]
