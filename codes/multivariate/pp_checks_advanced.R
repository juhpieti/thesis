###################################################################
### SCRIPT TO PERFORM MORE ADVANCED POSTERIOR PREDICTIVE CHECKS ###
### THE GOAL IS TO LOCATE SOURCES FOR WEIRD BEHAVIOR IN         ###
### SPECIES-LEVEL PREDICTIONS FOR 2-STAGE MODEL                 ###
###################################################################

# load in libraries
library(rstan)
library(terra)

# load in helper functions
load("codes/helpers.R")

# prepare data
load("data/estonia_new/train/train_2020_2021_all_species_n500.Rdata")
train <- train_n500_all_species

# prepare the covariate matrix
X <- train[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd) #approximate amount of light reaching the bottom
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

# scale the covariates
X.scaled <- scale_covariates(X)
# add second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))

### take the correct species
Y <- train[,20:71]
idx_positive <- colSums(Y > 0) > 2 # take only species that have some observations
sum(idx_positive) # it's 21 of those
Y_21 <- Y[,idx_positive]
Y <- Y_21[,!colnames(Y_21) == "Ranunculus peltatus subsp_ Baudotii"] # drop 1 species for convenient J = 20

pp_check_2stage_advanced <- function(mod,Y_mat,X,sp_name_list,n_rep=50,bin_width = 0.05,test_quantities=FALSE,save_plot=FALSE,save_loc="",im_width=800,im_height=600) {
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
  # save_p: TRUE = plot gets saved
  # save_loc: path to the saving location
  # im_height: height of the saved plot
  # im_width: width of the saved plot
  
  X <- as.matrix(X) # from df to matrix to work with matrix multiplication
  n_species <- length(sp_name_list)
  
  # extract posterior draws
  post.samples <- extract(mod)
  
  # make n_rep replications of new data, compare to observed coverages
  n_post_samples <- nrow(post.samples$alpha)
  rep_idx <- sample(1:n_post_samples,n_rep,replace = FALSE) #randomly take n_rep sets of posterior samples
  plot_idx <- sample(rep_idx,15,replace=FALSE) #take 15 to draw the histograms with observed dataset
  
  # initialize list to gather y_rep matrices (n_rep X n) for each species
  y_rep_list <- lapply(sp_name_list, function(x) list(Y = matrix(0,nrow=n_rep,ncol=nrow(X)),
                                                      Pi = matrix(0,nrow=n_rep,ncol=nrow(X)),
                                                      L_fe = matrix(0,nrow=n_rep,ncol=nrow(X)),
                                                      L_re = matrix(0,nrow=n_rep,ncol=nrow(X))))
  names(y_rep_list) <- sp_name_list
  
  # save also y_tot and alpha_0
  y_rep_list[["y_tot"]] <- matrix(0,nrow=n_rep,ncol=nrow(X))
  y_rep_list[["alpha_0"]] <- matrix(0,nrow=n_rep,ncol=nrow(X))
  y_rep_list[["Mu_tot"]] <- matrix(0,nrow=n_rep,ncol=nrow(X))
  
  
  ### MAKE PREDICTIONS
  # go through n_rep posterior samples to replicate n_rep datasets
  for (i in 1:n_rep) {
    idx <- rep_idx[i] # index for posterior sample
    
    # load the model parameters for j:th species and i:th posterior sample
    ### Negative-Binomial part of the model
    alpha_M <- post.samples$alpha_M[idx]
    beta_M <- c(post.samples$beta_M_1[idx,],post.samples$beta_M_2[idx,])
    scale_M <- post.samples$scale_NegBin[idx]
    
    f_M <- alpha_M + X %*% beta_M
    mu_M <- exp(f_M)
    
    ### Dirichlet-Multinomial part of the model
    alpha <- post.samples$alpha[idx, ]
    beta <- rbind(post.samples$beta_1[idx,,], post.samples$beta_2[idx,,])
    lambda <- post.samples$Lambda[idx,,]
    scale_Dir <- post.samples$scale_Dir[idx]
    
    # latent factors
    Z <- post.samples$Z[idx,,]
    
    # calculate latent f 
    f_fe <- alpha + X %*% beta
    f_re <- Z %*% lambda
    f <- f_fe + f_re # n x J matrix
    
    # calculate the corresponding mu
    mu <- t(apply(f,1,softmax)) ### n_pred x J matrix
    
    # sample Y
    M <- rnbinom(nrow(X),size=scale_M,mu=mu_M) # total coverage from Negative-Binomial
    Y_rep_mat <- t(sapply(1:length(M), function(i) rdirmn(1,M[i],scale_Dir*mu[i,]))) #Y_i from Dirichlet-Multinomial
    
    # save the results
    for (j in 1:n_species) {
      y_rep_list[[j]]$Y[i,] <- Y_rep_mat[,j]
      y_rep_list[[j]]$Pi[i,] <- mu[,j]
      y_rep_list[[j]]$L_fe[i,] <- f_fe[,j]
      y_rep_list[[j]]$L_re[i,] <- f_re[,j]
    }
    y_rep_list[["y_tot"]][i,] <- M
    y_rep_list[["alpha_0"]][i,] <- scale_Dir
    y_rep_list[["Mu_tot"]][i,] <- mu_M
  }
  
  ### go through species to produce plots
  for (j in 1:n_species) {
    # initialize vectors to save test statistics from the replicated data sets
    T_probzero <- c() # proportion of zeros
    T_mean <- c() # mean(y)
    T_mean_positive <- c() #mean(y) for all y > 0
    T_max <- c() # max(y)
    
    ### START PLOTTING
    # save if asked
    if(save_plot) {
      sp_name <- sp_names[j]
      sp_name_modified <- gsub(" ","_",sp_name)
      sp_name_modified <- gsub("/","_",sp_name_modified)
      png(paste0(save_loc,"pp_check_reps_",sp_name_modified,".png"),width = im_width, height = im_height)
    }
    
    # 4x4 grid for plotting
    par(mfrow = c(4,4),
        mar = c(2,4,2,0))
    
    # observations from j:th species
    y <- Y_mat[,j]
    
    # plot the observed data first
    if (!test_quantities) { 
      hist(y, breaks = c(-5,0,seq(5,100,by=5)), xlim = c(0,100), main = "obs", ylim = c(0,length(y)), freq = TRUE)
    }
    
    # go through n_rep posterior samples again
    for (i in 1:n_rep) {
      idx <- rep_idx[i] # index for posterior sample
      # take the i:th replication of j:th species
      y_rep <- y_rep_list[[j]]$Y[i,]
      
      # manual breaks
      breaks <- c(-bin_width, 0, seq(bin_width, 5*ceiling(max(100,y_rep)/5), by=bin_width))
      
      if (!test_quantities) {
        # plot if the index was chosen to be plotted
        if (idx %in% plot_idx) {
          hist(y_rep, breaks = breaks, xlim = c(0,5*ceiling(max(c(100,y_rep))/5)), main = paste0("rep",i), ylim = c(0,length(y)), freq = TRUE)
        }
      }
      
      # save the test statistics
      T_probzero <- c(T_probzero, mean(y_rep==0))
      T_mean <- c(T_mean, mean(y_rep))
      T_mean_positive <- c(T_mean_positive, ifelse(any(y_rep > 0), mean(y_rep[y_rep > 0]), 0))
      T_max <- c(T_max, max(y_rep))
    }
    
    if(save_plot) {
      dev.off()
    }
    
    ### draw histograms of test statistics (if test_quantities == TRUE)
    if (test_quantities) {
      # save if asked
      if(save_plot) {
        sp_name <- sp_names[j]
        sp_name_modified <- gsub(" ","_",sp_name)
        sp_name_modified <- gsub("/","_",sp_name_modified)
        png(paste0(save_loc,"pp_check_",sp_name_modified,".png"),width = im_width, height = im_height)
      }
      
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
      
      if(save_plot) {
        dev.off()
      }
    }
  }
  return(y_rep_list)
}


n_rep <- 200
set.seed(123)
mod.2stage <- readRDS(paste0("models/two_stage/",subfolder,"/M1/JSDM_NegBin_DirMult_hier_priors.RDS"))
tmp <- pp_check_2stage_advanced(mod.2stage,Y,X.sec_ord,sp_names,n_rep,5,test_quantities = TRUE,save_plot = FALSE)

# examine the results

# check that actually sum_j(y_ij) = y_tot?
sp_names <- colnames(Y)
for (i in 1:2) {
  y_tot <- tmp$y_tot[i,]
  sum_yij <- rep(0,nrow(X))
  for(sp_name in sp_names) {
    sum_yij <- sum_yij + tmp[[sp_name]]$Y[i,]
  }
  print(cbind(y_tot,sum_yij)) ### works
}


# take the species that show awkward behavior
sp_of_int <- c("Chara aspera","Chara canescens","Myriophyllum spicatum","Potamogeton perfoliatus")
#sp_of_int <- c("Chara aspera")

# if y_ij > y.treshold, plot some of the interesting quantities
y.treshold <- 450

res.tbl <- c()

for (s in 1:length(sp_of_int)) {
  sp_name <- sp_of_int[s]
  species_idx <- which(sp_names == sp_name)
  for (r in 1:n_rep) {
    for (i in 1:nrow(X)) {
      if (tmp[[sp_name]]$Y[r,i] >= y.treshold) {
        res.tbl <- rbind(res.tbl, data.frame(sp_idx=species_idx,
                                    mu_tot=tmp$Mu_tot[r,i],
                                    y_tot=tmp$y_tot[r,i],
                                    y_j=tmp[[sp_name]]$Y[r,i],
                                    pi=tmp[[sp_name]]$Pi[r,i],
                                    a0=tmp$alpha_0[r,i],
                                    X[i,,drop=FALSE]))
      }
    }
  }
}

colnames(res.tbl)[ncol(res.tbl)] <- "light level"
res.tbl

# plot these conditions wrt. distribution of the data (look at the response curves)
par(mfrow = c(3,3))
for(i in 1:ncol(X)) {
  hist(X[,i], breaks = 40, main = colnames(X)[i], xlab = "value")
  points(res.tbl[,colnames(X)[i]],rep(1,nrow(res.tbl)),pch = 19, cex = 0.5, col = "red")
}


### EXAMINE ROWS 15, 26, 36
### ONE FOR EACH SPECIES, TAKE POSTERIOR SAMPLES AND CALCULATE THE SOFT-MAX VALUES
### JUST USING THE ENVIRONMENTAL RESPONSES
row_idx <- c(15,26,36,40)
#X.examples <- res.tbl[row_idx,colnames(X)]
#X.ex.scaled <- scale_covariates(X,X.examples)
post.sam <- extract(mod.2stage, pars = c("alpha","beta_1","beta_2"))

for (i in row_idx) {
  f.vals <- c()
  for (j in 1:ncol(Y)) {
    # extract posterior means
    alpha <- mean(post.sam$alpha[,j])
    beta1 <- colMeans(post.sam$beta_1[,,j])
    beta2 <- colMeans(post.sam$beta_1[,,j])
    X.sc <- scale_covariates(X,res.tbl[i,colnames(X)])
    f.vals <- c(f.vals, alpha + sum(X.sc*beta1) + sum(beta2*X.sc^2))
  }
  print("f")
  print(f.vals)
  print("softmax(f)")
  print(softmax(f.vals))
  print("pi_j in sample")
  print(res.tbl[i,"pi"])
  print("pi_j using post means")
  print(softmax(f.vals)[res.tbl[i,"sp_idx"]])
}

### look more closely at these locations: what are the other posterior samples saying about E[pi_ij]?
loc_covariates <- res.tbl[15,colnames(X)]
loc_idx <- which(apply(X,1,function(x) sum(x == loc_covariates) == ncol(X)))
loc_idx

# print all the Pi samples from this locations for Chara aspera, which showed weird behavior
tmp$`Chara aspera`$Pi[,loc_idx] #almost all of it is \approx 1
tmp$`Myriophyllum spicatum`$Pi[,loc_idx] #why is this not the dominant species, using posterior means it is?

# examine the fixed effect part for f
tmp$`Chara aspera`$L_fe[,loc_idx]
tmp$`Myriophyllum spicatum`$L_fe[,loc_idx]

tmp$`Chara aspera`$L_re[,loc_idx]
tmp$`Myriophyllum spicatum`$L_re[,loc_idx]

# examine the expected cover in that location
tmp$y_tot[,loc_idx]


### PLOT FIXED EFFECT PART FOR EACH SPECIES?
fe.mat <- c()
for(sp_name in sp_names) {
  fe.mat <- cbind(fe.mat, tmp[[sp_name]]$L_fe[,loc_idx])
}
colnames(fe.mat) <- sp_names
colMeans(fe.mat)
boxplot(fe.mat)
abline(h=0,col = "red", lty = 2)

### PLOT RANDOM EFFECT PART FOR EACH SPECIES?
re.mat <- c()
for(sp_name in sp_names) {
  re.mat <- cbind(re.mat, tmp[[sp_name]]$L_re[,loc_idx])
}
colnames(re.mat) <- sp_names
colMeans(re.mat)
boxplot(re.mat)
abline(h=0,col = "red", lty = 2)

### PLOT ALSO VARIATION IN EXPECTED COVER IN THIS LOCATION AND LETS SAY 9 OTHERS? IS THERE CONSIDERABLE VARIATION?
EYtot.mat <- c()
set.seed(123)
rand_loc_idx <- sample(1:nrow(X),8,TRUE)
for (i in c(loc_idx,rand_loc_idx)) {
  EYtot.mat <- cbind(EYtot.mat, tmp$Mu_tot[,i])
}
colnames(EYtot.mat) <- c(loc_idx,rand_loc_idx)
par(mfrow = c(1,1))
boxplot(EYtot.mat, xlab = "location idx")

### MAKE THIS SAME EXPERIMENT FOR SOME AMOUNT OF LOCATIONS, HOW IT SEEMS CONTRADICTIONALLY THAT THE RARE SPECIES WINS IN LOCATIONS WITH HIGH EXPECTED TOTAL COVER?
pdf("plots/final_results/two_stage/n_500/M1/pp_checks/adv_checks_FE_RE.pdf")
par(mfrow = c(3,1))
set.seed(123)
for(i in sample(1:500,40,FALSE)) {
  ###
  y_tot_mean <- mean(tmp$Mu_tot[,i])
  y_tot_sd <- sd(tmp$Mu_tot[,i])
  
  ### PLOT THE PROPORTION
  pi.mat <- c()
  for(sp_name in sp_names) {
    pi.mat <- cbind(pi.mat, tmp[[sp_name]]$Pi[,i])
  }
  colnames(pi.mat) <- sp_names
  colMeans(pi.mat)
  boxplot(pi.mat, main = paste0("location ",i,", total cover: ", round(y_tot_mean,2)," +- 2*", round(y_tot_sd)), ylab = "E[Pi]", las = 2, cex.axis = 0.5)
  
  ### PLOT FIXED EFFECT PART FOR EACH SPECIES?
  fe.mat <- c()
  for(sp_name in sp_names) {
    fe.mat <- cbind(fe.mat, tmp[[sp_name]]$L_fe[,i])
  }
  colnames(fe.mat) <- sp_names
  colMeans(fe.mat)
  #boxplot(fe.mat, main = paste0("location ",i,", total cover: ", round(y_tot_mean,2)," +- 2*", round(y_tot_sd)), ylab = "FE", las = 2, cex.axis = 0.5)
  boxplot(fe.mat, ylab = "FE", las = 2, cex.axis = 0.5)
  abline(h=0,col = "red", lty = 2)
  
  ### PLOT RANDOM EFFECT PART FOR EACH SPECIES?
  re.mat <- c()
  for(sp_name in sp_names) {
    re.mat <- cbind(re.mat, tmp[[sp_name]]$L_re[,i])
  }
  colnames(re.mat) <- sp_names
  colMeans(re.mat)
  boxplot(re.mat, ylab = "RE", las = 2, cex.axis = 0.5)
  abline(h=0,col = "red", lty = 2)
}
dev.off()

### EXAMINE THE VARIANCES IN COEFFICIENTS WITH E.G CHARA ASPERA AND MYTILUS TROSSULUS
colnames(Y)
sp_idx_chara <- 3
sp_idx_myrio <- 10
sp_idx_mytilus <- 11

summary(mod.2stage, pars = c(paste0("beta_1[",1:9,",",sp_idx_chara,"]"),paste0("beta_2[",1:9,",",sp_idx_chara,"]")))$summary
summary(mod.2stage, pars = c(paste0("beta_1[",1:9,",",sp_idx_myrio,"]"),paste0("beta_2[",1:9,",",sp_idx_myrio,"]")))$summary
summary(mod.2stage, pars = c(paste0("beta_1[",1:9,",",sp_idx_mytilus,"]"),paste0("beta_2[",1:9,",",sp_idx_mytilus,"]")))$summary

scale_covariates(X,loc_covariates)

# remove to save space
rm(tmp)