library(terra)
library(rstan)
library(ggplot2)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("codes/helpers.R")

# load in the data
hanko2 <- read.csv("data/hanko2/velmu.to.Obama.2024-09-05.csv", header = TRUE)

# take only vegetation mapping (exact)
hanko2 <- hanko2[hanko2[,"MethodCategory"] == "Vegetation mapping", ]
# take only Dive or Transect as suggested by Louise
hanko2 <- hanko2[hanko2$Method %in% c("Dive","Transect"), ]
# drop the covers unidentified
hanko2 <- hanko2[hanko2$Cover != 9999, ]
# drop the ID columns
hanko2 <- hanko2[,3:ncol(hanko2)]

# turn data into wide format using ID:s (which ID? Transect?)
head(hanko2)

cover_vec <- c()
for (rec_ID in unique(hanko2$RecordID)) {
  sub_df <- hanko2[hanko2$RecordID == rec_ID, ]
  cover_vec <- c(cover_vec, sum(sub_df$Cover,na.rm = TRUE))
}

hist(cover_vec, breaks = 40)

# The RecordID seems to be the one telling which rows correspond to one dive
hanko_wide_X <- c()
hanko_wide_Y <- c()
n_species <- length(unique(hanko2$scientificName))
sp_names <- unique(hanko2$scientificName)

for (rec_ID in unique(hanko2$RecordID)) {
  Y_vec <- rep(0,n_species)
  names(Y_vec) <- sp_names
  sub_df <- hanko2[hanko2$RecordID == rec_ID, ]
  Y_vec[sub_df$scientificName] <- sub_df$Cover
  hanko_wide_Y <- rbind(hanko_wide_Y, Y_vec)
  
  hanko_wide_X <- rbind(hanko_wide_X, cbind(sub_df[1,c(1:6,9:14,16:20)]))
}

hanko_wide <- cbind(hanko_wide_X, hanko_wide_Y)

hist(rowSums(hanko_wide[,18:ncol(hanko_wide)] > 0))
hist(rowSums(hanko_wide[,18:ncol(hanko_wide)]), breaks = 40)

save(hanko_wide, file = "data/hanko2/hanko2_wide.Rdata")

### visualize on map!
### take a subset of n=500, randomly? any clever column to use for subsetting?

# first drop NA's
hanko_wide <- hanko_wide[complete.cases(hanko_wide), ]

#plot the data
hanko_wide.vect <- vect(hanko_wide, geom = c("X_coord","Y_coord"), crs = "EPSG:3067")
europe.vect <- vect("europe_coastline_shapefile/Europe_coastline_poly.shp")
europe.vect <- project(europe.vect, "EPSG:3067")
hanko.rast <- rast("data/hanko1/sal_rasters_Obama.2023-06-19.tif")
europe.vect <- crop(europe.vect, hanko.rast)

# histograms of variables in full data
par(mfrow = c(2,3))
for (i in c(1:5,11)) {
  hist(hanko_wide[,i], main = colnames(hanko_wide)[i], breaks = 5, col = "forestgreen",
       probability = TRUE, xlab = "value")
}

# histogram of prevalences in full data
# subset only interesting species
species_of_interest <- c("Fucus","Fucus vesiculosus", "Potamogeton", "Potamogeton perfoliatus",
                         "Zostera marina", "Chara", "Chara aspera", "Chara canescens")

prev_full <- 100*colMeans(hanko_wide[,species_of_interest] > 0)
par(mfrow = c(1,1))
hist(prev_full, main = "Prevalences in full data", breaks = 10, col = "forestgreen",
     xlab = "Prevalence (%)", probability = TRUE, xlim = c(0,25))


# could I subset by year?
dates <- hanko_wide$Date
years <- substr(dates,1,4)

set.seed(13)
#hanko_wide_sub <- hanko_wide[years %in% c("2014","2015"), ]
#hanko_wide_sub <- hanko_wide_sub[sample(1:nrow(hanko_wide_sub),500,FALSE), ]

# or just randomly?
set.seed(1234)
hanko_wide_sub <- hanko_wide[sample(1:nrow(hanko_wide),500,FALSE),]

hanko_wide_sub_vect <- vect(hanko_wide_sub, geom = c("X_coord","Y_coord"), crs = "EPSG:3067")

plot(hanko.rast, col = "grey")
plot(europe.vect, add = TRUE)
plot(hanko_wide.vect, col = "red", add = TRUE, cex = 1)
plot(hanko_wide_sub_vect, col = "blue", add = TRUE, cex = 0.8)

# histograms of covariates
par(mfrow = c(2,3))
for (i in c(1:5,11)) {
  hist(hanko_wide_sub[,i], main = colnames(hanko_wide_sub)[i], breaks = 5, col = "forestgreen",
       xlab = "value", probability = TRUE)
}

# histograms of prevalences
par(mfrow = c(1,1))
prev_sub <- 100*colMeans(hanko_wide_sub[,species_of_interest] > 0)
hist(prev_sub, main = "Prevalences in the subset (2014-2015)", breaks = 30, col = "forestgreen",
     xlab = "Prevalence (%)", probability = TRUE, xlim = c(0,25))

### start the analysis

X <- hanko_wide_sub[,c(1:5,11)]
X.scaled <- scale_covariates(X)

hanko_res_loo <- c()

for (sp_name in species_of_interest) {
  y <- hanko_wide_sub[,sp_name]
  y.bin <- as.numeric(y > 0)
  
  ### fit the model
  mod <- fit_logistic_regression(y.bin,X.scaled,4,1000,FALSE)
  #mod.iid <- fit_logistic_regression(y.bin,X.scaled,4,1000,TRUE)
  
  ### check convergence
  check_convergence(mod,FALSE)
  
  ### calculate loo
  mod.loo <- calc_loo(mod)
  # save to a table
  hanko_res_loo <- rbind(hanko_res_loo, c(mod.loo))
  
  ### examine the coefficient distributions
  sp_name_modified <- gsub(" ","_",sp_name)

  png(paste0("plots/hanko/regression/coefficient_distributions/",sp_name_modified,".png"))
  plot_distributions(mod,X.scaled,FALSE,3,4)
  dev.off()

  ### examine the responses
  png(paste0("plots/hanko/regression/response_curves/",sp_name_modified,".png"))
  plot_responses_hanko(mod,X.scaled,X,TRUE,FALSE,0,1,2,3)
  dev.off()
}
