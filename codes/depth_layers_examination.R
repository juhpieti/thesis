### check how close my grid centers are from the depth layer
### how does this affect to the depth e.g. between my observed data and the layer I have

plot(depth.vect.3067.cropped, col = "red", cex = .1)
plot(europe.vect.cut, add = TRUE)
plot(estonia_sub.vect)


### zoom in, plot depth and grid centers?
ext_vec_small_area <- ext(c(200000,300000,6430000,6500000))
depth.vect.small.area <- crop(depth.vect.3067.cropped, ext_vec_small_area)

predictive_grid.small.area <- crop(predictive_grid, ext_vec_small_area)

copernicus_grid_points.vect.sub.cropped <- crop(copernicus_grid_points.vect.sub, ext(predictive_grid))
copernicus_grid_points.vect.small.area <- crop(copernicus_grid_points.vect.sub, ext_vec_small_area)

plot(copernicus_grid_points.vect.sub.cropped, cex = 0.1)

plot(depth.vect.small.area, cex = .3, col = "red")
plot(centroids(predictive_grid.small.area), add = TRUE, col = "blue", cex = 0.3)
#plot(copernicus_grid_points.vect.sub, add = TRUE, col = "forestgreen", cex=0.3)
#plot(europe.vect.cut, add = TRUE)

### see how close my predictive points are to a depth layer
depth.vect.small.area
predictive_grid.small.area
nearest_from_grid_to_depth <- nearest(predictive_grid.small.area,depth.vect.small.area)
hist(nearest_from_grid_to_depth$distance, breaks = 40)


### see how it affects for example the depth of my observed points
estonia_sub.vect
observed_depths <- estonia_sub.vect$depth
depths_from_grid <- predictive_grid$depth[nearest_from_obs_to_pred_grid$to_id]
depths_from_depth_grid <- depth.vect.3067.cropped$depth.mean[nearest_from_obs_to_depth_grid$to_id]

depth_df <- cbind(observed_depths,-1*depths_from_grid,-1*depths_from_depth_grid)

nearest_from_obs_to_pred_grid <- nearest(estonia_sub.vect,predictive_grid)
nearest_from_obs_to_depth_grid <- nearest(estonia_sub.vect,depth.vect.3067.cropped)

nearest_from_obs_to_pred_grid
