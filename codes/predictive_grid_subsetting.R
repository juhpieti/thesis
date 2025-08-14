####################################################################################
### 1) CREATES A COARSE (20km x 20km) GRID FOR SPATIAL RANDOM EFFECTS            ###
### 2) CREATES THE PREDICTIVE GRID OVER REGION OF INTERST (1km x 1km resolution) ###
####################################################################################

# load in package to manage spatial data
library(terra)

### load in the training data
load("data/estonia_new/train_2020_2021.Rdata")

### visualize the data spatially
europe.vect <- vect("europe_coastline_shapefile/Europe_coastline_poly.shp") #coastline of Europe as polygon
europe.vect <- project(europe.vect, "EPSG:3067") #reproject to TM35FIN

# crop for some region of interest
#ext_vec <- c(0,566613.8,6000000,6650000)
#ext_vec <- c(100000,600000,6300000,6650000)
ext_vec <- c(100000,630000,6300000,6650000)

europe.vect.cut <- crop(europe.vect,ext_vec)

# turn the dataframe into GIS-type vector
estonia_sub.vect <- vect(df_sub, geom = c("x","y"), crs = "EPSG:3067")

# plot the data
par(mfrow = c(1,1))
plot(estonia_sub.vect, cex = 1, col = "red", main = "observations")
plot(europe.vect.cut, add = TRUE, col = "lightgrey")
plot(estonia_sub.vect, cex = 1, col = "red", main = "observations", add = TRUE)

### 1) Create coarse grid for spatial random effects

### even grid over the study area
ext_vec_grid <- ext(c(134000,614000,6318000,6638000))
#ext_vec_grid <- ext(c(140000,560000,6320000,6640000)) #original!!!
length_grid_cell <- 20000 #20km x 20km grid
ncols <- ( ext_vec_grid[2] - ext_vec_grid[1] ) / length_grid_cell
nrows <- ( ext_vec_grid[4] - ext_vec_grid[3] ) / length_grid_cell

grid_raster <- rast(ext_vec_grid,ncols=ncols,nrows=nrows, crs = "EPSG:3067")

# erase the land areas
grid_vector <- as.polygons(grid_raster)
spatial_grid <- erase(grid_vector, europe.vect.cut)

# plot the grid
par(mfrow = c(1,1))
plot(spatial_grid, main = "prediction grid")
plot(europe.vect.cut, add = TRUE, col = "lightgrey")
plot(estonia_sub.vect, add = TRUE, col = "red")

# save the grid
writeVector(spatial_grid, filename = "data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp", overwrite = TRUE)


### 2) PREPARE DEAPTH LAYER FOR THE EXTENT OF THIS SPATIAL GRID
depth <- readRDS("data/estonia_new/depth_SWM_SWM-Bekkby.rds")

# vector version
depth.vect <- vect(depth, geom = c("x","y"), crs = "EPSG:3035")
depth.vect <- project(depth.vect, "EPSG:3067") #reproject to TM35FIN

## crop by the extent of spatial grid
depth.vect.3067.cropped <- crop(depth.vect, ext(spatial_grid))

# # raster version
# res <- 1000 # 1km x 1km grid?
# ext_depth <- c(min(depth$x),max(depth$x),min(depth$y),max(depth$y))
# r <- rast(extent = ext_depth, crs = "EPSG:3035", res = res) ### which coordinate system is this?
# depth.rast <- rasterize(as.matrix(depth[,c("x","y")]),r,values=depth$depth.mean)
# depth.rast_3067 <- project(depth.rast, "EPSG:3067")

# data frame version
depth.df <- as.data.frame(depth.vect, geom = "XY")
depth_grid_3067 <- depth.df[,c("id_bs1km","x","y","depth.mean")]
save(depth_grid_3067, file = "data/estonia_new/depth_grid_3067_coords.Rdata")

### test that the locations match...
# depth.rast_3067_cropped <- crop(depth.rast_3067, ext(spatial_grid))
# plot(depth.rast_3067_cropped)
# plot(estonia_sub.vect, col = "red", add = TRUE)


### ONE IDEA TO EASE THE COMPUTATIONAL BURDEN IS TO DELETE THE DEEP SEA, WHERE WATER PLANTS SHOULD NOT APPEAR ###
### PREPARE A DEEP SEA LAYER
deep.sea <- depth.rast_3067_cropped < -70 # denote everything deeper than 70m as "deep"
deep.sea.vect <- as.polygons(deep.sea)
deep.sea.vect <- deep.sea.vect[deep.sea.vect$last == 1] #take only deep sea
test.spat.grid <- erase(spatial_grid, deep.sea.vect) #remove spatial grid cells that are denoted "deep sea"

### DOES NOT REALLY EFFECT TO THE NUMBER OF SPATIAL GRID CELLS:
spatial_grid
test.spat.grid


### PREPARE COPERNICUS LAYER (environmental covariates) FOR THE EXTENT OF THIS SPATIAL GRID

# load in the data
load("data/estonia_new/copernicus_grid_2021_July.Rdata")

# read in the copernicus station locations and add TM35FIN coordinates
copernicus_grid_points <- readRDS("data/estonia_new/BGC_003_012 and PHY_003_011 Copernicus grid points.rds")
copernicus_grid_points.vect <- vect(copernicus_grid_points[,1:3], geom = c("longitude_wgs84","latitude_wgs84"), crs = "EPSG:4326") # original grid points are in WGS84
copernicus_grid_points.vect <- project(copernicus_grid_points.vect, "EPSG:3067") # reproject to TM35FIN

# take a subset of only those stations appearing in the predictive grid
copernicus_grid_points.vect.sub <- terra::subset(copernicus_grid_points.vect, id_copernicus_phy_bgc %in% unique(grid_2021_july$id_copernicus_phy_bgc), NSE = T)

# turn into a data frame
CGP <- as.data.frame(copernicus_grid_points.vect, geom = "XY")

### now go through the copernicus grid and add the corresponding TM35FIN (EPSG:3067) coordinates
grid_2021_july[,c("x_3067","y_3067")] <- 0
for (i in 1:nrow(grid_2021_july)) {
  cop_id <- grid_2021_july[i,"id_copernicus_phy_bgc"] #ID of copernicus station
  coords <- CGP[cop_id,c("x","y")] # coordinates of that station
  grid_2021_july[i,c("x_3067","y_3067")] <- coords
}

head(grid_2021_july) #now the new coordinates should appear

# take only variables of interest
vars <- c("no3_bottom","o2_bottom","po4_bottom","zsd","bottomT","so_bottom","current_bottom","chl_bottom")
grid_2021_july <- grid_2021_july[,c("id_copernicus_phy_bgc","year","month","x_3067","y_3067",vars)]

# make it terra vector
grid_2021_july_vect <- vect(grid_2021_july, geom = c("x_3067","y_3067"), crs = "EPSG:3067")

# crop using the extent of predictive grid
grid_2021_july_vect_cropped <- crop(grid_2021_july_vect, ext(spatial_grid))


### NOW WE HAVE DEPTH AND COPERNICUS COVARIATES FOR THE AREA OF INTEREST (EXT(SPATIAL_GRID))
par(mfrow = c(1,1))
depth.vect.3067.cropped
plot(depth.vect.3067.cropped, cex = 0.1, col = "red")
plot(europe.vect.cut, add = TRUE)

grid_2021_july_vect_cropped
plot(grid_2021_july_vect_cropped, cex = 0.1, col = "red")
plot(europe.vect.cut, add = TRUE)


### what are their resolutions? check by looking at a small area
### depth should be 1km x 1km!
ext_vec_small_area <- ext(c(200000,340000,6425000,6510000))
ext_vec_small_area <- ext(c(330000,570000,6570000,6637900))
ext_vec_small_area <- ext(c(200000,430000,6400000,6637000))

# depth
depth.vect.small.area <- crop(depth.vect.3067.cropped, ext_vec_small_area)
plot(depth.vect.small.area, cex = 0.2)

coords_depth <- crds(depth.vect.small.area)
dist_depth <- as.matrix(dist(coords_depth))
diag(dist_depth) <- NA

# examine the distance to nearest grid point
depth_dist_nearest <- apply(dist_depth, 1, min, na.rm = TRUE)
hist(depth_dist_nearest) # seems to be the 1km as expected

# copernicus variables
cop.vect.small.area <- crop(grid_2021_july_vect_cropped, ext_vec_small_area)
plot(cop.vect.small.area, cex = 0.2)

coords_copernicus <- crds(cop.vect.small.area)
dist_copernicus <- as.matrix(dist(coords_copernicus))
diag(dist_copernicus) <- NA

# examine distance to nearest copernicus station
cop_dist_nearest <- apply(dist_copernicus,1,min,na.rm=TRUE)
hist(cop_dist_nearest[cop_dist_nearest < 1800], breaks = 40, probability = TRUE,
     main = "distances to nearest copernicus station", xlab = "distance (m)") # something like 1.6km


### CREATE 2 PREDICTIVE GRIDS ##
# 1) with deep sea, 1km x 1km
# 2) without deep sea, 500m x 500m DOES NOT MAKE SENSE!!! SINCE 1km and 1.6km are the resolutions for covariates

### MAKE AN 1km x 1km PREDICTION GRID
# parameters of the grid
spatial_grid <- vect("data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp")
ext_vec_grid <- ext(spatial_grid) # use the same extent than the spatial grid
length_grid_cell <- 1000 #1kmx1km grid
ncols <- ( ext_vec_grid[2] - ext_vec_grid[1] ) / length_grid_cell
nrows <- ( ext_vec_grid[4] - ext_vec_grid[3] ) / length_grid_cell

grid_raster <- rast(ext_vec_grid,ncols=ncols,nrows=nrows, crs = "EPSG:3067")

# erase land areas
grid_vector <- as.polygons(grid_raster)
cropped_grid <- erase(grid_vector, europe.vect.cut)

predictive_grid_1km <- cropped_grid

# save for later
writeVector(predictive_grid_1km, filename = "data/estonia_new/predictive_grid_1km_wo_variables/predictibe_grid_1km_wo_variables.shp")


### ADD NEAREST VALUES OF DEPTH AND COPERNICUS VARIABLES
predictive_grid_1km <- vect("data/estonia_new/predictive_grid_1km_wo_variables/predictibe_grid_1km_wo_variables.shp")
predictive_grid <- predictive_grid_1km

grid_2021_july_vect_cropped # here you have the copernicus values
depth.vect.3067.cropped # here you have the depth values

# Find the nearest grid points
nearest_depth <- nearest(predictive_grid, depth.vect.3067.cropped)
nearest_copernicus <- nearest(predictive_grid, grid_2021_july_vect_cropped)

# add nearest depths
predictive_grid$depth <- depth.vect.3067.cropped$depth.mean[nearest_depth$to_id]
predictive_grid$depth <- -1*predictive_grid$depth #make depths positive as in training data

# add nearest copernicus values
vars <- c("no3_bottom","o2_bottom","po4_bottom","zsd","bottomT","so_bottom","current_bottom","chl_bottom")
predictive_grid[,vars] <- grid_2021_july_vect_cropped[nearest_copernicus$to_id,vars]

# save the vector file
writeVector(predictive_grid,filename = "data/estonia_new/predictive_grid_1km_all_variables_2021_july/predictive_grid_1km_all_variables_2021_july.shp",
            overwrite = TRUE)


### PREPARE ALSO ONE PREDICTIVE GRID WITHOUT THE DEEP SEA
predictive_grid_wo_deep_sea <- erase(predictive_grid, deep.sea.vect)
writeVector(predictive_grid_wo_deep_sea,filename = "data/estonia_new/predictive_grid_1km_all_variables_2021_july_wo_deepsea/predictive_grid_1km_all_variables_2021_july_wo_deepsea.shp",
            overwrite = TRUE)




### VISUALIZE THE COVARIATES AS MAPS
vect_grid <- vect("data/estonia_new/predictive_grid_1km_all_variables_2021_july/predictive_grid_1km_all_variables_2021_july.shp")
#vect_grid$depth <- -1*vect_grid$depth
vect_grid$depth_to_secchi <- vect_grid$depth / vect_grid$zsd
vect_grid$light_bottom <- exp(-1.7*vect_grid$depth/vect_grid$zsd)
vect_rast <- rast(ext = ext(vect_grid), res = 1000, crs = "EPSG:3067") #1km x 1km resolution

#vars <- c("depth","no3_bottom","o2_bottom","po4_bottom","zsd","bottomT","so_bottom","current_bottom","chl_bottom")
vars <- c("depth","no3_bottom","o2_bottom","po4_bottom","zsd","bottomT","so_bottom","current_bottom","chl_bottom","depth_to_secchi","light_bottom")

for (var in vars) {
  print(var)
  r <- rasterize(vect_grid, vect_rast, field = var)
  png(paste0("plots/estonia_new/predictive_grid_1km_new/",var,".png"))
  plot(r)
  #plot(r.sub)
  dev.off()
}