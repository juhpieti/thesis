#########################################################
### Script to create different sizes of training data ###
#########################################################

# load in libraries
library(terra) # to handle spatial data (vector/raster)

# Load in the full data
df <- readRDS("data/estonia_new/2024_QuantitativeSamplesCoverKeySpeciesWithCopernicusDataAndDepth.rds")
table(df$year)

### Training data was decided to be taken from years 2020-2021
### Years 2022 onwards did not have the copernicus data
### Two years (instead of just the latest 2021) was taken to get more variation in the environmental covariates
df <- df[df[,"year"] <= 2021,]
df <- df[df[,"year"] >= 2020, ]

df <- as.data.frame(df)
table(df$year)
table(df$month)

### Take only summer months
df <- df[df[,"month"] %in% 5:9, ]

### change the coordinates to TM35FIN
df.vect <- vect(df, geom = c("longitude","latitude"), crs = "EPSG:4326") # original data is in WGS84
df.vect <- project(df.vect, "EPSG:3067") # reproject to TM35FIN

# back to dataframe
df <- as.data.frame(df.vect, geom = "XY")

# take the covariates of interest
vars <- c("depth","no3_bottom","o2_bottom","po4_bottom","zsd","bottomT","so_bottom","current_bottom","chl_bottom")

# and the species of interest
sp_list <- c("Amphibalanus improvisus","Chara aspera","Chara tomentosa",
             "Chorda filum","Cladophora glomerata","Fucus vesiculosus",
             "Furcellaria lumbricalis","Furcellaria lumbricalis loose form",
             "Myriophyllum spicatum","Mytilus trossulus","Potamogeton perfoliatus",
             "Pylaiella/Ectocarpus","Ruppia maritima","Stuckenia pectinata",
             "Tolypella nidifica","Ulva intestinalis","Vertebrata fucoides",
             "Zannichellia palustris","Zostera marina")

df <- df[,c(colnames(df)[1:8],"x","y",vars,sp_list)]

# save the full training data
train_full <- df
save(train_full, file="data/estonia_new/train/train_2020_2021_full.Rdata")

### take subsets of different sizes
set.seed(123)
n <- 500
rand_idx <- sample(1:nrow(train_full),n)
train_n500 <- train_full[rand_idx, ]
save(train_n500, file = "data/estonia_new/train/train_2020_2021_n500.Rdata")

set.seed(123)
n <- 1000
rand_idx <- sample(1:nrow(train_full),n)
train_n1000 <- train_full[rand_idx, ]
save(train_n1000, file = "data/estonia_new/train/train_2020_2021_n1000.Rdata")

set.seed(123)
n <- 2000
rand_idx <- sample(1:nrow(train_full),n)
train_n2000 <- train_full[rand_idx, ]
save(train_n2000, file = "data/estonia_new/train/train_2020_2021_n2000.Rdata")

set.seed(123)
n <- 100
rand_idx <- sample(1:nrow(train_full),n)
train_n100 <- train_full[rand_idx, ]
save(train_n100, file = "data/estonia_new/train/train_2020_2021_n100.Rdata")
