##########################################################
### SCRIPT TO MAKE VISUALIZATIONS OF THE TRAINING DATA ###
### n=500 from summer months in 2020-2021              ###
##########################################################

#load in packages
library(terra)

#load in the training data
load("data/estonia_new/train/train_2020_2021_n500.Rdata")
df_sub <- train_n500

# check how many presences in the data
colSums(df_sub[,20:38] > 0)
estonia_sub <- df_sub

# remove species with less than 5 (0, 1, 2 in this case) observations
estonia_sub <- estonia_sub[,!(colnames(estonia_sub) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]


sp_list <- colnames(estonia_sub)[20:35]
sp_list <- sp_list[c(2:6,10,13:14,7,9,11:12,15:16,1,8)] # order by the type of species (Vascular plants,  Algae etc.)

# plot the prevalences and the distribution of coverages for each species
prevalences <- round(100*colMeans(estonia_sub[,sp_list] > 0),2)

#manual breaks
bin_width <- 5
breaks <- c(-bin_width, 0, seq(bin_width, 100, by=bin_width)) #create own bar for zero observations by starting at -bin_width

im_width <- 800
im_height <- 600
subfolder <- paste0("n_",nrow(df_sub))

png(paste0("plots/estonia_new/training_data/",subfolder,"/counts_with_prevalences.png"), width = im_width, height = im_height)
par(mfrow = c(4,4),
    mar = c(4,4,2,1))
for (sp_name in sp_list) {
  hist(estonia_sub[,sp_name], main = sp_name, xlab = "coverage (%)", ylab = "count", ylim = c(0,nrow(df_sub)), xlim = c(0,100), breaks = breaks)
  legend("topright", bty = "n", legend = paste0(prevalences[sp_name],"%"), col = "red")
}
dev.off()

### the same without zeros
#manual breaks
bin_width <- 5
breaks <- c(0, seq(bin_width, 100, by=bin_width))

png(paste0("plots/estonia_new/training_data/",subfolder,"/positive_coverages.png"), width = im_width, height = im_height)
par(mfrow = c(4,4),
    mar = c(4,4,2,1))
for (sp_name in sp_list) {
  hist(estonia_sub[estonia_sub[,sp_name] > 0, sp_name], main = sp_name, ylab = "count", xlab = "coverage (%)", xlim = c(0,100), breaks = breaks)
}
dev.off()

### plot each species coverages against each covariate
X <- estonia_sub[,c(11:19)]
# X$depth_to_secchi <- X$depth/X$zsd
X$light_bottom <- exp(-1.7*X$depth/X$zsd)

for (i in 20:35) {
  y <- estonia_sub[,i]
  sp_name <- colnames(estonia_sub)[i]
  sp_name <- gsub(" ","_",sp_name)
  sp_name <- gsub("/","_",sp_name)
  
  png(paste0("plots/estonia_new/training_data/",subfolder,"/species_against_covariates/",sp_name,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  for(j in c(1:4,6:10)) {
    plot(X[,j],y, xlab = colnames(X)[j], ylab = colnames(estonia_sub)[i], pch = 18, cex = 1,ylim=c(0,100))
  }
  dev.off()
}


##### visualize predictive grid values by histograms
load("data/estonia_new/predictive_grid_1km_all_variables_2021_july.Rdata")
dim(pred_grid_1km_2021_july_df)
colnames(pred_grid_1km_2021_july_df)

# predictive grid vs training data
par(mfrow = c(3,3))
vars <- colnames(pred_grid_1km_2021_july_df)[2:10]
for (var_name in vars) {
  hist(pred_grid_1km_2021_july_df[,var_name], probability = TRUE, breaks = 50,
       xlab = "value", ylab = var_name, main = var_name, col = "forestgreen")
  hist(estonia_sub[,var_name], probability = TRUE, breaks = 50,
       xlab = "value", ylab = var_name, main = var_name, col = "royalblue", add = TRUE)
}

# just the training data
for (var_name in vars) {
  hist(estonia_sub[,var_name], probability = TRUE, breaks = 50,
       xlab = "value", ylab = var_name, main = var_name)
}


### visualize study area and observations on map
# load in the training data
# either n = 2000
load("data/estonia_new/train/train_2020_2021_n2000.Rdata")
df_sub <- train_n2000

# or n = 500 (used in modeling)
load("data/estonia_new/train/train_2020_2021_n500.Rdata")
df_sub <- train_n500

# load in the spatial grid
spatial_grid <- vect("data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp")

estonia_sub.vect <- vect(df_sub, geom = c("x","y"), crs = "EPSG:3067")
train.vect.lonlat <- project(estonia_sub.vect, "EPSG:4326")
#ext(estonia_sub.vect) <- ext(spatial_grid)

### load in the coastline as polygon
europe.vect <- vect("europe_coastline_shapefile/Europe_coastline_poly.shp")
europe.vect.lonlat <- project(europe.vect, "EPSG:4326")
europe.vect <- project(europe.vect, "EPSG:3067") #reproject to TM35FIN

### visualize spatial random effect grid
subfolder <- paste0("n_",nrow(df_sub))
png(paste0("plots/estonia_new/training_data/",subfolder,"/spatial_effect_grid.png"), width = im_width, height = im_height)
par(mfrow = c(1,1))
plot(spatial_grid, xlab = "Easting (m)", ylab = "Northing (m)")
plot(europe.vect, col = "lightgrey", add = TRUE)
plot(estonia_sub.vect, add = TRUE, col = "red", cex = 0.7, alpha = 0.8)

### add scale bar and north arrow
sbar(xy = "bottomright", type = "bar", divs = 10, scaleby = 1000, below = "km")
north(xy = c(600000,6370000),type = 2)

dev.off()

### same but without the grid on the background (just observations)
subfolder <- paste0("n_",nrow(df_sub))
png(paste0("plots/estonia_new/training_data/",subfolder,"/study_area_close.png"), width = im_width, height = im_height)
par(mfrow = c(1,1))
#plot(spatial_grid, xlab = "Easting (m)", ylab = "Northing (m)", border = NA, background = "lightblue1")
plot(spatial_grid, xlab = "Easting (m)", ylab = "Northing (m)", border = NA)
#plot(europe.vect.cut, col = "darkseagreen3", add = TRUE)
plot(europe.vect, col = "lightgrey", add = TRUE)
plot(estonia_sub.vect, add = TRUE, col = "red", cex = 0.7, alpha = 0.8)

### add scale bars
sbar(d = 100000, xy = "bottomright", type = "bar", divs = 4, scaleby = 1000, below = "km")
north(xy = c(600000,6370000),type = 2)
dev.off()

### create a large scale map to see better where the study are is located
# crop for some region of interest in longitudes and latitudes
#ext_vec_lonlat <- c(-25,60,45,80)
ext_vec_lonlat <- c(-25,60,47,73)
#ext_vec_lonlat <- c(-25,55,50,75)

# cut the coastline shapefile
europe.vect.lonlat.cut <- crop(europe.vect.lonlat, ext_vec_lonlat)

# plot
subfolder <- paste0("n_",nrow(df_sub))
png(paste0("plots/estonia_new/training_data/",subfolder,"/study_area.png"), width = im_width, height = im_height)
par(mfrow = c(1,1))
plot(europe.vect.lonlat.cut, col = "lightgrey", xlab = "longitude (degrees)", ylab = "latitude (degrees)")
plot(train.vect.lonlat, cex = 0.4, alpha = 0.8, col = "red", add = TRUE)

# add the scale bar and north arrow
sbar(d = 1000, xy = "bottomright", type = "bar", divs = 4, below = "km", lonlat = TRUE, lwd = 1, adj = c(1,-1.5))
north(xy = c(-23,71),type = 2)

### draw extent of spatial grid (extent of the study area)
spatial_grid.lonlat <- project(spatial_grid, "EPSG:4326")
plot(as.polygons(ext(spatial_grid.lonlat)), add = TRUE, border = "blue", lwd = 2)
dev.off()
