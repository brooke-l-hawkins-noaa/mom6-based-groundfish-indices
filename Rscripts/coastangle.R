# Install required packages if you don't have them:
# install.packages(c("terra", "sf", "tidyverse"))

library(terra)
library(sf)     
library(tidyverse)
devtools::load_all(here::here())

# 1. Load the Bathymetry NetCDF file
nc_path <- "STATIC/ocean_static.deptho.nc"

# Read the depth sub-dataset
depth_raster <- rast(nc_path, subds = "deptho")

# 2. Dynamically create a land/ocean mask from the depth profiles
# Deep ocean values are positive, while land is typically 0, negative, or NA.
# We create a binary mask: 1 for water (depth > 0), 0 for land.
mask_raster <- depth_raster > 0

# 3. Extract the Coastline Pixels using Vector Polygons
# Convert the binary mask to polygon features to grab clean border outlines
mask_poly <- as.polygons(mask_raster, values = TRUE)

# Use names() instead of colnames() for SpatVector objects
names(mask_poly) <- "is_water"

# Filter for the water polygon (where value is 1) to build a boundary line
water_poly <- mask_poly[mask_poly$is_water == 1, ] 

# Cast the outer edge of the water polygon into a line feature string
coastline_sf <- st_as_sf(water_poly) %>% 
  st_geometry() %>% 
  st_cast("MULTILINESTRING")


# 4. Extract Spatial Coordinates Along the Identified Coastline
coast_points <- st_coordinates(st_cast(coastline_sf, "POINT")) %>% 
  as.data.frame() %>% 
  rename(lon = X, lat = Y) %>% 
  select(lon, lat) # Clean up formatting artifacts

## plot to check how it looks
ggplot() +
  geom_sf(data = coastline_sf, color = "darkblue", size = 0.8) +
  theme_minimal() +
  labs(
    title = "Extracted Coastline Geometry",
    subtitle = "Generated from ocean static depth profile",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_line(color = "gray90")
  )

ggplot() + 
  geom_sf(data = coastline_sf, color = 'darkblue', size = 0.8) + 
  coord_sf(
    xlim = c(230,260), 
    ylim = c(20, 48), 
    expand = FALSE
  ) + 
  theme_minimal() + 
  labs( 
    title = "Extracted Coastline Geometry", 
    subtitle = "Generated from ocean static depth profile", 
    x = "Longitude", 
    y = "Latitude" 
  ) + 
  theme( 
    plot.title = element_text(face = "bold", size = 14), 
    panel.grid.major = element_line(color = "gray90") 
  )

# 5. Function to calculate local orientation angles via Moving Windows
# Apply across all valid localized shoreline coordinates
# 5.5 Filter data to be South of Canada (< 48°N) and North of Baja California (> 22°N)
coast_points_filtered <- coast_points %>%
  filter(lat >= 32 & lat <= 48 & lon>=230)

coast_points_filtered$coast_angle_rad <- sapply(1:nrow(coast_points_filtered), function(i) {
  get_coast_angle(i, coast_points_filtered, window_size = 15)
})

# Alternative: Map coastal points and color them by their orientation angle
ggplot(coast_points_filtered, aes(x = lon, y = lat, color = coast_angle_rad)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "plasma", name = "Angle (rad)") +
  theme_minimal() +
  labs(
    title = "Coastal Orientations Along Shoreline",
    x = "Longitude",
    y = "Latitude"
  )

# View the final processed dataframe containing coordinates and orientation metrics
head(coast_points_filtered)

# (Optional) If you also want to trim the line geometry for plotting:
coastline_sf_filtered <- st_crop(coastline_sf, xmin = -180, xmax = -100, ymin = 32, ymax = 48)

# 6. Export the dataframe as an R object (.rds)
# Replace "coastal_angles.rds" with your desired filename and path
saveRDS(coast_points_filtered, file = "STATIC/coastal_angles.rds")


