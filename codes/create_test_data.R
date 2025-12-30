######################################
### SCRIPT TO CREATE THE TEST DATA ###
######################################

#load in packages
library(terra)

### load in the full data
load("data/estonia_new/train/train_2020_2021_full_all_species.Rdata")
dim(df_all_species)

### load in the training data
load("data/estonia_new/train/train_2020_2021_all_species_n500.Rdata")
dim(train_n500_all_species)

### the training data are identified by row name
### this information can be used to exclude training data from the full data
length(unique(rownames(train_n500_all_species)))

### all training points are included in the full data set
sum(rownames(train_n500_all_species) %in% rownames(df_all_species))

### CREATE TEST SET OF SIZE m
m <- 500
train_idx <- rownames(train_n500_all_species)
test_full_all_species <- df_all_species[!(rownames(df_all_species) %in% train_idx), ] #drop the training points

set.seed(123)
test_idx <- sample(1:nrow(test_full_all_species), m)
test_n500_all_species <- test_full_all_species[test_idx, ]

### save the test sets
save(test_full_all_species, file = "data/estonia_new/test/test_2020_2021_all_species_full_except_n500.Rdata")
save(test_n500_all_species, file = "data/estonia_new/test/test_2020_2021_all_species_n500.Rdata")

