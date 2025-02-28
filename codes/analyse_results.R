### models fitted with M1-M4.R scripts
### this script loads fitted models and draws different maps and curves, as well as produces tables


### leave-one out CV

loo_table <- c()
loo(mod_amphi)$estimates[1,1]
loo(mod_amphi.ZI)$estimates[1,1]
loo(mod_amphi.spat)$estimates[1,1]
loo(mod_amphi.spat.ZI)$estimates[1,1]



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
