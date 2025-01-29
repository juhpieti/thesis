### take a relevant subset of the huge estonia dataset

# load in packages
library(terra)

# load in the data
library(readxl)
estonia <- read_excel("data/estonia/2024_QuantitativeSamplesCoverKeySpecies.xlsx")
estonia <- as.data.frame(estonia)
head(estonia)

# change the CRS system to same as with Finnish data
estonia.vect <- vect(estonia, geom = c("longitude","latitude"), crs = "EPSG:4326") #original data is in WGS84
estonia.vect <- project(estonia.vect, "EPSG:3067")

# take a subset of the data
# first spatially
europe.vect <- vect("europe_coastline_shapefile/Europe_coastline_poly.shp")
europe.vect <- project(europe.vect, "EPSG:3067")

plot(estonia.vect, cex = 0.5, main = "Estonia observations")
ext_vec <- c(0,566613.8,6000000,6650000)
estonia.vect.cut <- crop(estonia.vect,ext_vec)
europe.vect.cut <- crop(europe.vect,ext_vec)
plot(estonia.vect.cut, cex = 0.5, main = "Estonia observations cropped", col = "red")
plot(europe.vect.cut, add = TRUE)

estonia <- as.data.frame(estonia.vect.cut)
estonia[,c("x","y")] <- crds(estonia.vect.cut)

# take only the interesting species
sp_list <- c("Amphibalanus improvisus","Chara aspera","Chara tomentosa",
             "Chorda filum","Cladophora glomerata","Fucus vesiculosus",
             "Furcellaria lumbricalis","Furcellaria lumbricalis loose form",
             "Myriophyllum spicatum","Mytilus trossulus","Potamogeton perfoliatus",
             "Pylaiella/Ectocarpus","Ruppia maritima","Stuckenia pectinata",
             "Tolypella nidifica","Ulva intestinalis","Vertebrata fucoides",
             "Zannichellia palustris","Zostera marina")

estonia <- estonia[,c(colnames(estonia)[1:16], sp_list, "x", "y")]

### visualize the full data
# distribution of depth
par(mfrow = c(1,1))
hist(estonia[,"depth"], probability = TRUE, xlab = "depth", main = "Full data distribution of depth", breaks = 40)

# distribution of bottom types
par(mfrow = c(3,4),
    mar = c(4,4,2,0))
for (col_name in colnames(estonia)[7:16]) {
  hist(estonia[,col_name], main = col_name, xlab = "%", probability = TRUE, breaks = 20)
}

# distribution of species coverages
par(mfrow = c(4,5),
    mar = c(4,4,2,0))

prevalences <- round(100*colMeans(estonia[,sp_list] > 0),2)

for (sp_name in sp_list) {
  hist(estonia[,sp_name], main = sp_name, xlab = "%", xlim = c(0,100), probability = TRUE, breaks = 20)
  text(50,0.1,paste0(prevalences[sp_name],"%"), col = "red")
}

### visualize the subset (by years) of the data
# first subset with year(s)
set.seed(123)
estonia_sub <- estonia[estonia[,"year"] %in% c(2022), ] #year 2022
estonia_sub <- estonia_sub[sample(1:nrow(estonia_sub), 500, replace = FALSE), ] #random sample with n = 500
estonia_sub.vect <- vect(estonia_sub, geom = c("x","y"), crs = "EPSG:3067") 

#without subsetting by years
# set.seed(123)
# estonia_sub <- estonia[sample(1:nrow(estonia), 500, replace = FALSE), ]
# estonia_sub.vect <- vect(estonia_sub, geom = c("x","y"), crs = "EPSG:3067")

# visualize along the whole data spatially, is well distributed across the area?
par(mfrow = c(1,1))
plot(estonia.vect.cut, cex = 0.5, main = "Estonia observations cropped", col = "red")
plot(europe.vect.cut, add = TRUE)
plot(estonia_sub.vect, cex = 1, col = "blue", add = TRUE)

# distribution of depth
hist(estonia_sub[,"depth"], breaks = 40, main = "Subset distribution of depth", xlab = "depth")

# distribution of bottom types
par(mfrow = c(3,4),
    mar = c(4,4,2,0))
for (col_name in colnames(estonia_sub)[7:16]) {
  hist(estonia_sub[,col_name], main = col_name, xlab = "%", probability = TRUE, breaks = 20)
}

# distribution of species coverages

par(mfrow = c(4,5),
    mar = c(4,4,2,0))

prevalences <- round(100*colMeans(estonia_sub[,sp_list] > 0),2)

for (sp_name in sp_list) {
  hist(estonia_sub[,sp_name], main = sp_name, xlab = "%", xlim = c(0,100), probability = TRUE)
  text(50,0.05,paste0(prevalences[sp_name],"%"), col = "red")
}

### save the data
save(estonia_sub, file = "data/estonia_sub/estonia_sub_df.RData")
writeVector(estonia_sub.vect, filename = "data/estonia_sub/estonia_sub_vect.shp")

