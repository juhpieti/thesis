##################################################
### SCRIPT TO DRAW CORRELATION PLOTS FOR JSDMs ###
##################################################
library(rstan)
library(corrplot)

draw_corr_plot <- function(stan_fit, sp_names, n_factors, binary = FALSE) {
  ### binary: TRUE ==> instead of posterior mean correelations, plots positive with blue and negative with red (>95 post. probability)
  post.samples.array <- as.array(stan_fit)
  n_species <- length(sp_names)
  n_samples <- dim(post.samples.array)[1]
  n_chains <- dim(post.samples.array)[2]
  R_hats <- c()
  res_list <- list("mean(R)" = matrix(0,nrow=n_species,ncol=n_species,dimnames = list(sp_names,sp_names)),
                   "Pr(R>0)" = matrix(0,nrow=n_species,ncol=n_species,dimnames = list(sp_names,sp_names)),
                   "Pr(R<0)" = matrix(0,nrow=n_species,ncol=n_species,dimnames = list(sp_names,sp_names)))
  
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
      R_hats[paste0("S_",i,j)] <- Rh
      
      # save posterior mean, Pr(R_ij > 0) and Pr(R_ij < 0)
      res_list[[1]][i,j] <- mean(R_ij)
      res_list[[2]][i,j] <- mean(R_ij > 0)
      res_list[[3]][i,j] <- mean(R_ij < 0)
    }
  }
  
  # 1) just posterior means
  par(mfrow = c(1,1))
  if (binary) {
    # 2) different colors for 1) positive correlation 2) negative correlation 3) no correlation (overlaps with 0)
    
    R_mat_categorical <- matrix(0,nrow=n_species,ncol=n_species,dimnames=list(sp_names,sp_names))
    R_mat_categorical[which(res_list$`Pr(R>0` > 0.95)] <- 1
    R_mat_categorical[which(res_list$`Pr(R<0` > 0.95)] <- -1
    corrplot(R_mat_categorical, type = "upper", order = "original")
  } else {
    corrplot(res_list$`mean(R`, type = "upper", order = "original")
  }
}

### load in models
mod.JSDM <- readRDS("models/multivariate/n_500/M1/JSDM_hier_priors.RDS")
mod.JSDM.spat <- readRDS("models/multivariate/n_500/M3/JSDM_spatial_hier_priors.RDS")
mod.2stage <- readRDS("models/two_stage/n_500/M1/JSDM_NegBin_DirMult_hier_priors.RDS")
mod.2stage.spat <- readRDS("models/two_stage/n_500/M3/JSDM_NegBin_DirMult_spatial_hier_priors.RDS")

### start plotting
im_width <- 800
im_height <- 600

### continuous correlation plots
png(paste0("plots/final_results/JSDM/",subfolder,"/M1/correlation_plot_continuous.png"), width = im_width, height = im_height)
draw_corr_plot(mod.JSDM,sp_names,2,binary=FALSE)
dev.off()

png(paste0("plots/final_results/JSDM/",subfolder,"/M3/correlation_plot_continuous.png"), width = im_width, height = im_height)
draw_corr_plot(mod.JSDM.spat,sp_names,2,binary=FALSE)
dev.off()

png(paste0("plots/final_results/two_stage/",subfolder,"/M1/correlation_plot_continuous.png"), width = im_width, height = im_height)
draw_corr_plot(mod.2stage,sp_names,2,binary=FALSE)
dev.off()

png(paste0("plots/final_results/two_stage/",subfolder,"/M3/correlation_plot_continuous.png"), width = im_width, height = im_height)
draw_corr_plot(mod.2stage.spat,sp_names,2,binary=FALSE)
dev.off()

### binary correlation plots
png(paste0("plots/final_results/JSDM/",subfolder,"/M1/correlation_plot_binary.png"), width = im_width, height = im_height)
draw_corr_plot(mod.JSDM,sp_names,2,binary=TRUE)
dev.off()

png(paste0("plots/final_results/JSDM/",subfolder,"/M3/correlation_plot_binary.png"), width = im_width, height = im_height)
draw_corr_plot(mod.JSDM.spat,sp_names,2,binary=TRUE)
dev.off()

png(paste0("plots/final_results/two_stage/",subfolder,"/M1/correlation_plot_binary.png"), width = im_width, height = im_height)
draw_corr_plot(mod.2stage,sp_names,2,binary=TRUE)
dev.off()

png(paste0("plots/final_results/two_stage/",subfolder,"/M3/correlation_plot_binary.png"), width = im_width, height = im_height)
draw_corr_plot(mod.2stage.spat,sp_names,2,binary=TRUE)
dev.off()
