library(terra)

grid_2021 <- read.csv("data/estonia_new/Copernicus_grid_2021.csv")
dim(grid_2021)
colnames(grid_2021)
months <- grid_2021$month
table(months)

### load in the training data
load("data/estonia_new/train_2020_2021.Rdata")
table(df_sub$month)

### use July = 7 values in predictive grid, from the latest year 2021
grid_2021 <- grid_2021[months == 7, ]

save(grid_2021,file="data/estonia_new/copernicus_grid_2021_July.Rdata")

### prepare grid by removing extra covariates and adding TM35FIN coordinates
load("data/estonia_new/copernicus_grid_2021_July.Rdata")

copernicus_grid_points <- readRDS("data/estonia_new/BGC_003_012 and PHY_003_011 Copernicus grid points.rds")

copernicus_grid_points.vect <- vect(copernicus_grid_points[,1:3], geom = c("longitude_wgs84","latitude_wgs84"), crs = "EPSG:4326")
copernicus_grid_points.vect <- project(copernicus_grid_points.vect, "EPSG:3067")

copernicus_grid_points.vect.sub <- terra::subset(copernicus_grid_points.vect, id_copernicus_phy_bgc %in% unique(grid_2021$id_copernicus_phy_bgc), NSE = T)


CGP <- as.data.frame(copernicus_grid_points.vect, geom = "XY")

### now go through the grid_points and add the corresponding TM35FIN (EPSG:3067) coordinates
grid_2021[,c("x","y")] <- 0
for (i in 1:nrow(grid_2021)) {
  cop_id <- grid_2021[i,"id_copernicus_phy_bgc"]
  coords <- CGP[cop_id,c("x","y")]
  grid_2021[i,c("x","y")] <- coords
}

### save the predictive grid
# take only interesting variables
vars <- c("no3_bottom","o2_bottom","po4_bottom","zsd","bottomT","so_bottom","current_bottom","chl_bottom")
grid_2021 <- grid_2021[,c("id_copernicus_phy_bgc","year","month","x","y",vars)]
save(grid_2021, file = "data/estonia_new/copernicus_grid_2021_July_3067_coords.Rdata")


### prepare the predictive grid for depth
depth <- readRDS("data/estonia_new/depth_SWM_SWM-Bekkby.rds")
depth.vect <- vect(depth, geom = c("x","y"), crs = "EPSG:3035")
depth.vect <- project(depth.vect, "EPSG:3067")

depth.df <- as.data.frame(depth.vect, geom = "XY")
depth_grid <- depth.df[,c("id_bs1km","x","y","depth.mean")]

save(depth_grid, file = "data/estonia_new/depth_grid_3067_coords.Rdata")
