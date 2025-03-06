### models fitted with M1-M4.R scripts
### this script loads fitted models and draws different maps and curves, as well as produces tables

# load in the training data
load("data/estonia_new/train_2020_2021.Rdata")
colnames(df_sub)
colSums(df_sub[,20:38] > 0)

# remove species that are too rare (under 5 observations)
train <- df_sub # more comfortable name
train <- train[,!(colnames(train) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# prepare the covariate matrix
X <- train[,11:19]
#X$depth_to_secchi <- X$depth / X$zsd # add secchi/depth for a variable representing seafloor light level
X$light_bottom <- exp(-1.7*X$depth / X$zsd)
### leave-one out CV

# load in the predictive grid
predictive_grid <- vect("data/estonia_new/predictive_grid_1km_all_variables_2021_july/predictive_grid_1km_all_variables_2021_july.shp")
load("data/estonia_new/predictive_grid_1km_all_variables_2021_july.Rdata")
dim(pred_grid_1km_2021_july_df)
colnames(pred_grid_1km_2021_july_df)

# add relative light level
pred_grid_1km_2021_july_df$light_bottom <- exp(-1.7*pred_grid_1km_2021_july_df$depth/pred_grid_1km_2021_july_df$zsd)

### load in spatial grid
spatial_grid <- vect("data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp")

grid_centers <- centroids(spatial_grid)
grid_centers.df <- as.data.frame(grid_centers, geom = "XY")
grid_centers.df <- grid_centers.df[,c("x","y")]

dim(grid_centers.df) # there are m = 191 grid cells

### find the observed grid cells
estonia_sub.vect <- vect(train, geom = c("x","y"), crs = "EPSG:3067")

nearest_grid_center <- nearest(estonia_sub.vect, grid_centers)
nearest_grid_center.df <- as.data.frame(nearest_grid_center)

# n-length vector indicating the ID of the nearest grid center
nearest_grid_center.vec <- nearest_grid_center.df$to_id

# take the indexes of grid cells that has observations in them
observed_grid_cells <- unique(nearest_grid_center.vec)
observed_grid_cells.df <- grid_centers.df[observed_grid_cells,c("x","y")]

# create the P matrix
P <- matrix(0,ncol=length(observed_grid_cells),nrow=nrow(train))
colnames(P) <- rownames(observed_grid_cells.df)
for (i in 1:nrow(estonia_sub)) {
  P[i,as.character(nearest_grid_center.vec[i])] <- 1
}

### put the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000


############### START ANALYZING THE RESULTS ####################

### calculate loo-CV

### MODEL 1 (BASE)
loo.beta <- c()

sp_names <- colnames(train)[20:35]

for(sp_name in sp_names[1:16]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  mod.beta <- readRDS(paste0("models/M1/",sp_name_modified,".rds"))
  
  # calculate leave-one-out LPD
  loo.beta <- c(loo.beta,calc_loo(mod.beta))
  
}

### MODEL 2
loo.ZIbeta <- c()

sp_names <- colnames(train)[20:35]

for(sp_name in sp_names[1:16]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  mod.ZIbeta <- readRDS(paste0("models/M2/",sp_name_modified,".rds"))
  
  # calculate leave-one-out LPD
  loo.ZIbeta <- c(loo.ZIbeta, calc_loo(mod.ZIbeta))

}

### MODEL 3
loo.beta_spat <- c()

sp_names <- colnames(train)[20:35]

for(sp_name in sp_names[1:16]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  mod.beta_spat <- readRDS(paste0("models/M3/",sp_name_modified,".rds"))
  
  # calculate leave-one-out LPD
  #loo.beta_spat <- c(loo.beta_spat,calc_loo(mod.beta_spat))
  
  print(summary(mod.beta_spat, pars = c("rho","l","s2_cf"))$summary)
}

### MODEL 4
loo.ZIBeta_spat <- c()

sp_names <- colnames(train)[20:35]

for(sp_name in sp_names[1:16]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  mod.ZIBeta_spat <- readRDS(paste0("models/M4/",sp_name_modified,".rds"))
  
  # calculate leave-one-out LPD
  #loo.ZIBeta_spat <- c(loo.ZIBeta_spat, calc_loo(mod.ZIBeta_spat))
  print(summary(mod.ZIBeta_spat, pars = c("rho","l_mu","s2_mu","l_pi","s2_pi"))$summary)
}

### Combine results
loo_table <- cbind(loo.beta,loo.ZIbeta,loo.beta_spat,loo.ZIBeta_spat)

colnames(loo_table) <- c("base","ZI","RE","ZI+RE")
rownames(loo_table) <- sp_names

save(loo_table,file = "results/final_results/loo_table.Rdata")

### to get the latex output
library(xtable)
print(xtable(loo_table, type = "latex"))


### PRODUCE MAPS (expected coverages, hotspots)

sp_names <- colnames(train)[20:35]

thinning <- 10
hotspot_limit <- 0.7

im_width <- 800
im_height <- 600

for(sp_name in sp_names[1:2]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  ### load in models
  mod.beta <- readRDS(paste0("models/M1/",sp_name_modified,".rds"))
  mod.ZIbeta <- readRDS(paste0("models/M2/",sp_name_modified,".rds"))
  mod.beta_spat <- readRDS(paste0("models/M3/",sp_name_modified,".rds"))
  mod.ZIBeta_spat <- readRDS(paste0("models/M4/",sp_name_modified,".rds"))
  
  
  ### print summaries of key parameters?
  #print(summary(mod.beta_spat, pars = c("rho","l","s2_cf"))$summary)
  
  ### make predictions over predictive grid
  pred_list.beta <- predict_beta_regression(mod.beta,
                                            pred_grid_1km_2021_july_df[,c(2:10,13)],
                                            X,thinning)
  pred_list.ZIbeta <- predict_ZI_beta_regression(mod.ZIbeta,
                                             pred_grid_1km_2021_july_df[,c(2:10,13)],
                                             X,thinning)
  pred_list.beta_spat <- predict_spatial_beta_regression(mod.beta_spat,
                                                         pred_grid_1km_2021_july_df[,c(2:10,13)],
                                                         X,
                                                         pred_grid_1km_2021_july_df[,c("x","y")],
                                                         grid_centers.df/1000, observed_grid_cells.df,thinning,1)
  pred_list.ZIBeta_spat <- predict_spatial_ZI_beta_regression(mod.ZIBeta_spat,
                                                              pred_grid_1km_2021_july_df[,c(2:10,13)],
                                                              X,
                                                              pred_grid_1km_2021_july_df[,c("x","y")],
                                                              grid_centers.df/1000, observed_grid_cells.df,thinning,1)
  
  ### Plot coverages
  png(paste0("plots/final_results/coverage_maps/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_map(pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","BASE")
  plot_map(pred_list.ZIbeta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","ZI")
  plot_map(pred_list.beta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","RE")
  plot_map(pred_list.ZIBeta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","ZI+RE")
  dev.off()
  
  
  ### Plot hotspots
  png(paste0("plots/final_results/hotspot_maps/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_map(pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "BASE",hotspot_limit)
  plot_map(pred_list.ZIbeta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "ZI",hotspot_limit)
  plot_map(pred_list.beta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "RE",hotspot_limit)
  plot_map(pred_list.ZIBeta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "ZI+RE",hotspot_limit)
  dev.off()
}


plot_map <- function(pred_list, locs, pred.grid.vect, type, title_chr = "", hotspot_proportion = 0.7) {
  ### type: "cover", "hotspot", 
  
  vals <- colMeans(pred_list$EY_sam)
  vars <- apply(pred_list$EY_sam,2,var)
  test_df <- as.data.frame(cbind(locs, vals, vars))
  
  ### accumulative
  test_df_ordered <- test_df[order(-test_df$vals), ]
  test_df_ordered$hotspot <- 1
  
  # find when 70% of the expected value is reached
  tot_sum <- sum(test_df_ordered$vals)
  idx <- which(cumsum(test_df_ordered$vals) > hotspot_proportion*tot_sum)[1]
  
  # set all outside to be non-hotspot
  test_df_ordered[idx:nrow(test_df_ordered),"hotspot"] <- 0
  
  test_vect <- vect(test_df_ordered, geom=c("x","y"), crs = "EPSG:3067")
  test_rast <- rast(ext = ext(pred.grid.vect), res = 1000, crs = "EPSG:3067")
  r <- rasterize(test_vect, test_rast, field = "vals")
  r.vars <- rasterize(test_vect, test_rast, field = "vars")
  r.hotspot <- rasterize(test_vect, test_rast, field = "hotspot")

  if (type == "cover") {
    plot(r, colNA = "lightgrey", main = title_chr, ,plg = list(title = "E[coverage]"),
         xlab = "Easting (m)", ylab = "Northing (m)")
  } else if (type == "hotspot") {
    plot(r.hotspot, colNA = "lightgrey", main = title_chr, col = c("red","blue"),
         plg = list(title = "hotspot", legend = c("no","yes")),
         xlab = "Easting (m)", ylab = "Northing (m)")
  }
}

plot_and_save_responses <- function(mod_list,X,grid_length = 200,im_width,im_height) {
  
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)

  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  post.mean_list <- list() #list with n_model components, each are grid_length x p matrices
  prob.zero_list <- list()
  
  mod_idx <- 1
  for(mod in mod_list) {
    # save the posterior means
    post.mean_mat <- c() #grid_length x p matrix 
    prob.zero_mat <- c()
    
    # save all the predictions for current model
    all_vars_preds <- list() #list with p components, each n_rep x grid_length
    
    # make predictions
    par(mfrow = c(3,4),
        mar = c(4,4,2,0))
    ### plots if ZI model
    if (mod_idx %in% c(2:4)) {
      png(paste0("plots/final_results/M",mod_idx,"/separate_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    }

    for (i in 1:ncol(X)) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      if (mod_idx %in% c(1,3)) {
        res <- predict_beta_regression(mod,predX,X,thinning,1)
      } else {
        res <- predict_ZI_beta_regression(mod,predX,X,thinning,1)
      }
      res.EY <- res$EY_sam # n_rep x grid_length
      all_vars_preds[[colnames(X)[i]]] <- res.EY
      post.means <- colMeans(res.EY)
      
      #plot the probability of presences and expectation given presence
      if (mod_idx %in% c(2,4)) {
        
        pi <- colMeans(res$prob_suit_sam)
        EY_if_suitable <- colMeans(res$EY_if_suitable_sam)
        
        plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), main = colnames(X)[i],type="l")
        lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
        # add legend
        plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
        legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
      }
      if (mod_idx %in% c(2,4)) {
        dev.off()
      }
      
      # save posterior means
      post.mean_mat <- cbind(post.mean_mat, post.means)
      prob.zero_mat <- cbind(prob.zero_mat, colMeans(res$probzero_sam))
    }
    

    # now we want to normalize such that largest posterior men gets value 1!
    max_post_mean <- max(post.mean_mat)
    post.mean_mat <- post.mean_mat / max_post_mean #now max value is 1
    
    # save scaled posterior means
    post.mean_list[[mod_idx]] <- post.mean_mat
    prob.zero_list[[mod_idx]] <- prob.zero_mat

    # scale also all the other samples
    all_vars_preds <- lapply(all_vars_preds, function(x) (x/max_post_mean))
  
    par(mfrow = c(3,4),
        mar = c(4,4,2,0))
    
    
    # draw plots
    png(paste0("plots/final_results/M",mod_idx,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    for (i in 1:ncol(X)) { # index for variable
      plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
           xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i])
      
      n_rep <- nrow(res.EY)
      for (j in 1:n_rep) { # index for sample
        lines(grid_matrix[,i],all_vars_preds[[i]][j,], col = "lightgrey", lwd = 0.8)
      }
      lines(grid_matrix[,i], post.mean_mat[,i], col = "black", lwd = 1)
    }
    dev.off()
    
    mod_idx <- mod_idx + 1
  }
  
  ### then also plot the colmeans on the same plot
  cols <- rainbow(length(post.mean_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  
  png(paste0("plots/final_results/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i], lwd = 1)
    for (j in 1:length(post.mean_list)) {
      lines(grid_matrix[,i],post.mean_list[[j]][,i], col = cols[j], lwd = 1)
    }
  }
  plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
  legend("center", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1.5, bty = "n")
  dev.off()

  par(mfrow = c(3,4),
      mar = c(4,4,2,0))
  
  png(paste0("plots/final_results/prob_zero_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  
  ### same for probability of zero?
  cols <- rainbow(length(prob.zero_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = "value", ylab = "prob. of zero", main = colnames(X)[i], lwd = 1)
    for (j in 1:length(prob.zero_list)) {
      lines(grid_matrix[,i],prob.zero_list[[j]][,i], col = cols[j], lwd = 1)
    }
  }
  plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
  legend("center", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1.5, bty = "n")
  dev.off()
}


### PRODUCE RESPONSE CURVES

sp_names <- colnames(train)[20:35]

thinning <- 25
hotspot_limit <- 0.7

im_width <- 800
im_height <- 600

for(sp_name in sp_names[1:1]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  ### load in models
  mod.beta <- readRDS(paste0("models/M1/",sp_name_modified,".rds"))
  mod.ZIbeta <- readRDS(paste0("models/M2/",sp_name_modified,".rds"))
  mod.beta_spat <- readRDS(paste0("models/M3/",sp_name_modified,".rds"))
  mod.ZIBeta_spat <- readRDS(paste0("models/M4/",sp_name_modified,".rds"))
  
  plot_and_save_responses(list(mod.beta,mod.ZIbeta,mod.beta_spat,mod.ZIBeta_spat),X,200,im_width,im_height)
}



### difference maps

par(mfrow = c(2,2))

locs <- pred_grid_1km_2021_july_df[,c("x","y")]
vals <- colMeans(pred_list_m2$EY_sam) - colMeans(pred_list_m1$EY_sam)

test_df <- as.data.frame(cbind(locs, vals))
test_vect <- vect(test_df, geom = c("x","y"), crs = "EPSG:3067")
test_rast <- rast(ext = ext(predictive_grid), res = 1000, crs = "EPSG:3067")
r.diff <- rasterize(test_vect, test_rast, field = "vals")

plot(r.diff, colNA = "lightgrey", main = "M2-M1 (ZI-base)",
     plg = list(title = "diff."),
     xlab = "Easting (m)", ylab = "Northing (m)")

locs <- pred_grid_1km_2021_july_df[,c("x","y")]
vals <- colMeans(pred_list_m3$EY_sam) - colMeans(pred_list_m1$EY_sam)
test_df <- as.data.frame(cbind(locs, vals))
test_vect <- vect(test_df, geom = c("x","y"), crs = "EPSG:3067")
test_rast <- rast(ext = ext(predictive_grid), res = 1000, crs = "EPSG:3067")
r.diff <- rasterize(test_vect, test_rast, field = "vals")

plot(r.diff, colNA = "lightgrey", main = "M3-M1 (spat-base)",
     plg = list(title = "diff."),
     xlab = "Easting (m)", ylab = "Northing (m)")

locs <- pred_grid_1km_2021_july_df[,c("x","y")]
vals <- colMeans(pred_list_m4$EY_sam) - colMeans(pred_list_m2$EY_sam)
test_df <- as.data.frame(cbind(locs, vals))
test_vect <- vect(test_df, geom = c("x","y"), crs = "EPSG:3067")
test_rast <- rast(ext = ext(predictive_grid), res = 1000, crs = "EPSG:3067")
r.diff <- rasterize(test_vect, test_rast, field = "vals")

plot(r.diff, colNA = "lightgrey", main = "M4-M2 ((ZI+spat)-ZI)",
     plg = list(title = "diff."),
     xlab = "Easting (m)", ylab = "Northing (m)")

locs <- pred_grid_1km_2021_july_df[,c("x","y")]
vals <- colMeans(pred_list_m4$EY_sam) - colMeans(pred_list_m3$EY_sam)
test_df <- as.data.frame(cbind(locs, vals))
test_vect <- vect(test_df, geom = c("x","y"), crs = "EPSG:3067")
test_rast <- rast(ext = ext(predictive_grid), res = 1000, crs = "EPSG:3067")
r.diff <- rasterize(test_vect, test_rast, field = "vals")

plot(r.diff, colNA = "lightgrey", main = "M4-M3 ((ZI+spat)-spat)",
     plg = list(title = "diff."),
     xlab = "Easting (m)", ylab = "Northing (m)")


par(mfrow = c(1,1))

locs <- pred_grid_1km_2021_july_df[,c("x","y")]
vals <- colMeans(pred_list_m4$EY_sam) - colMeans(pred_list_m1$EY_sam)
test_df <- as.data.frame(cbind(locs, vals))
test_vect <- vect(test_df, geom = c("x","y"), crs = "EPSG:3067")
test_rast <- rast(ext = ext(predictive_grid), res = 1000, crs = "EPSG:3067")
r.diff <- rasterize(test_vect, test_rast, field = "vals")

plot(r.diff, colNA = "lightgrey", main = "Amphi M4-M1 expected coverage")
