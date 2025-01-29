### script to visualize the data used for modeling
### this is a subset (n = 500) using year 2022 in the estonia data

#load in packages
library(terra)

#load("data/estonia_sub/estonia_sub_df.RData")
load("data/estonia_new/train_2020_2021.Rdata")

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
  hist(estonia_sub[,sp_name], main = sp_name, xlab = "%", xlim = c(0,100), probability = TRUE, breaks = 20)
  text(50,0.05,paste0(prevalences[sp_name],"%"), col = "red")
}

par(mfrow = c(4,4),
    mar = c(4,4,2,0))

### the same without zeros
for (sp_name in sp_list) {
  hist(estonia_sub[estonia_sub[,sp_name] > 0, sp_name], main = sp_name, xlab = "%", xlim = c(0,100), probability = TRUE, breaks = 20)
}

### pairplots for each species

for (i in 20:35) {
  y <- estonia_sub[,i]
  X <- estonia_sub[,c(11:19)]
  sp_name <- colnames(estonia_sub)[i]
  sp_name <- gsub(" ","_",sp_name)
  sp_name <- gsub("/","_",sp_name)
  
  png(paste0("plots/training_data/species_against_covariates/",sp_name))
  par(mfrow = c(3,3),
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


