### script to visualize the data used for modeling
### this is a subset (n = 500) using year 2022 in the estonia data

#load in packages
library(terra)

#load in the training data
load("data/estonia_new/train_2020_2021.Rdata") #df_sub

# check how many presences in the data
colSums(df_sub[,20:38] > 0)
estonia_sub <- df_sub

# remove species with less than 5 (0, 1, 2 in this case) observations
estonia_sub <- estonia_sub[,!(colnames(estonia_sub) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

par(mfrow = c(4,4),
    mar = c(4,4,2,0))

sp_list <- colnames(estonia_sub)[20:35]

# plot the prevalences and the distribution of coverages for each species
prevalences <- round(100*colMeans(estonia_sub[,sp_list] > 0),2)

for (sp_name in sp_list) {
  hist(estonia_sub[,sp_name], main = sp_name, xlab = "cover (%)", ylab = "count", xlim = c(0,100), breaks = 15)
  legend("topright", bty = "n", legend = paste0(prevalences[sp_name],"%"), col = "red")
  #text(50,0.05,paste0(prevalences[sp_name],"%"), col = "red")
}

### the same without zeros
par(mfrow = c(4,4),
    mar = c(4,4,2,0))

for (sp_name in sp_list) {
  hist(estonia_sub[estonia_sub[,sp_name] > 0, sp_name], main = sp_name, xlab = "%", xlim = c(0,100), probability = TRUE, breaks = 20)
}

### pairplots for each species
X <- estonia_sub[,c(11:19)]
X$depth_to_secchi <- X$depth/X$zsd

for (i in 20:35) {
  y <- estonia_sub[,i]
  sp_name <- colnames(estonia_sub)[i]
  sp_name <- gsub(" ","_",sp_name)
  sp_name <- gsub("/","_",sp_name)
  
  png(paste0("plots/estonia_new/training_data/species_against_covariates_2/",sp_name))
  par(mfrow = c(3,4),
      mar = c(4,4,2,0))
  for(j in 1:ncol(X)) {
    plot(X[,j],y, xlab = colnames(X)[j], ylab = colnames(estonia_sub)[i], pch = 18, cex = 1,ylim=c(0,100))
  }
  dev.off()
}

### maps of coverages!!!

# turn the training data into shapefile
load("data/estonia_new/train_2020_2021.Rdata")
df_sub <- df_sub[,!(colnames(df_sub) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]
est.vect <- vect(df_sub, geom = c("x","y"), crs = "EPSG:3067")

# load in the coastline shapefile
eur.coast <- vect("europe_coastline_shapefile/Europe_coastline_poly.shp")
eur.coast <- project(eur.coast, "EPSG:3067")

# cut for visualization
ext_vec <- c(170000,550000,6410000,6640000)
eur.coast <- crop(eur.coast,ext_vec)

own_palette <- colorRampPalette(c("red","yellow","green"))

for (sp_name in colnames(df_sub)[20:35]) {
  n_cols <- 1+max(est.vect[,sp_name])
  cols <- own_palette(n_cols)
  
  par(mfrow = c(1,1))
  plot(eur.coast, main = sp_name)
  plot(est.vect, col = cols[1+df_sub[,sp_name]], add = TRUE, cex = 1)
}



##### visualize predictive grid values
dim(pred_grid_2021_july_df)
colnames(pred_grid_2021_july_df)
par(mfrow = c(3,3))

vars <- colnames(pred_grid_2021_july_df)[2:10]
for (var_name in vars) {
  hist(pred_grid_2021_july_df[,var_name], probability = TRUE, breaks = 50,
       xlab = "value", ylab = var_name, main = var_name, col = "forestgreen")
  hist(estonia_sub[,var_name], probability = TRUE, breaks = 50,
       xlab = "value", ylab = var_name, main = var_name, col = "royalblue", add = TRUE)
}

for (var_name in vars) {
  hist(estonia_sub[,var_name], probability = TRUE, breaks = 50,
       xlab = "value", ylab = var_name, main = var_name)
}



##### visualize study area and observations on map
### load in the training data
load("data/estonia_new/train_2020_2021.Rdata")

spatial_grid <- vect("data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp")

estonia_sub.vect <- vect(df_sub, geom = c("x","y"), crs = "EPSG:3067")
#ext(estonia_sub.vect) <- ext(spatial_grid)

### visualize the data spatially
europe.vect <- vect("europe_coastline_shapefile/Europe_coastline_poly.shp")
europe.vect <- project(europe.vect, "EPSG:3067") #reproject to TM35FIN

par(mfrow = c(1,1))
plot(spatial_grid, xlab = "Easting (m)", ylab = "Northing (m)")
plot(europe.vect.cut, col = "lightgrey", add = TRUE)
plot(estonia_sub.vect, add = TRUE, col = "red", cex = 0.8)

### add scale bars?
sbar(xy = "bottomright", type = "bar", divs = 10, scaleby = 1000, below = "km")
north(xy = c(600000,6370000),type = 2)


# crop for some region of interest
#ext_vec <- c(0,566613.8,6000000,6650000)
#ext_vec <- c(100000,600000,6300000,6650000)
ext_vec_large <- c(-1000000,900000,6100000,8000000)

europe.vect.cut <- crop(europe.vect,ext_vec_large)

par(mfrow = c(1,1))
plot(europe.vect.cut, col = "lightgrey", xlab = "easting (m)", ylab = "northing (m)")
plot(estonia_sub.vect, cex = 1, col = "red", add = TRUE)


ext_vec <- c(100000,630000,6300000,6650000)
europe.vect.cut <- crop(europe.vect,ext_vec)

par(mfrow = c(1,1))
plot(spatial_grid, alpha = 100)
plot(estonia_sub.vect, cex = 1, col = "red", xlab  = "easting (m)", ylab = "northing (m)")
plot(europe.vect.cut, add = TRUE, col = "lightgrey")
plot(estonia_sub.vect, cex = 1, col = "red", main = "observations", add = TRUE)
