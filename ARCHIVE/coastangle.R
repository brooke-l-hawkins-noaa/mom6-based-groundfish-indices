# Install required packages if you don't have them:
# install.packages(c("terra", "sf", "tidyverse"))

library(terra)
library(sf)     
library(tidyverse)
devtools::load_all(here::here())

# 1. Load the Bathymetry NetCDF file
nc_path <- "STATIC/ocean_static.deptho.nc"

# Read the depth sub-dataset as a raster file
depth_raster <- rast(nc_path, subds = "deptho")

# plot the raster 
ocean_palette <- hcl.colors(100, palette = "Blues 3", rev = TRUE)
plot(depth_raster, 
     col = ocean_palette, 
     main = "Ocean Depth Bathymetry",
     xlab = "Longitude", 
     ylab = "Latitude")

# 2. Dynamically create a land/ocean mask from the depth profiles
# Deep ocean values are positive, while land is typically 0, negative, or NA.
# We create a binary mask: 1 for water (depth > 0), 0 for land.
mask_raster <- depth_raster > 0

# 3. Extract the coastline pixels using vector polygons
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

# 4. Extract spatial coordinates along the identified coastline
coast_points <- st_coordinates(st_cast(coastline_sf, "POINT")) %>% 
  as.data.frame() %>% 
  rename(lon = X, lat = Y) %>% 
  select(lon, lat) # Clean up formatting artifacts


# 5. Function to calculate local orientation angles via Moving Windows
# Apply across all valid localized shoreline coordinates
# 5.5 Filter data to be South of Canada (< 48°N) and North of Baja California (> 22°N)
coast_points_filtered <- coast_points %>%
  filter(lat >= 32 & lat <= 48 & lon>=230)%>%
  distinct(lat, .keep_all = TRUE) # Keeps only unique spatial pairs

coast_points_filtered$coast_angle_rad <- sapply(1:nrow(coast_points_filtered), function(i) {
  get_coast_angle(i, coast_points_filtered, window_size = 15)
})

# Alternative: Map coastal points and color them by their orientation angle
# note: some of these look like duplicated that 
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








# 1. Load Required Packages --------------------------------------------------
library(terra)
library(sf)
library(tidyverse)
devtools::load_all(here::here())

# 2. Load Bathymetry NetCDF --------------------------------------------------
nc_path      <- "STATIC/ocean_static.deptho.nc"
depth_raster <- rast(nc_path, subds = "deptho")

# 3. Process Mask & Detect Shoreline in the Raster Domain -------------------
# Create a binary water mask: 1 for water, 0 for land
mask_raster <- depth_raster > 0

# Detect the boundary pixels where land meets water
# This flags the exact raster cells that represent the coastline
coast_raster_boundary <- boundaries(mask_raster, inner = TRUE)

# Assign a temporary value to the identified coastline cells (1 = Coast)
coast_raster_boundary[coast_raster_boundary == 0] <- NA

# 4. Crop and Extract Pixels at the Native 1/10 Degree Grid -----------------
# Crop the raster to your US West Coast bounding box first
coast_raster_cropped <- crop(
  coast_raster_boundary, 
  ext(230, 260, 32, 48)
)

# Extract the exact center coordinates of the 1/10 degree boundary pixels
# as.points() converts each non-NA cell to its precise spatial grid center point
coast_points_vect <- as.points(coast_raster_cropped, values = FALSE)

coast_points_filtered <- as.data.frame(crds(coast_points_vect)) %>% 
  rename(lon = x, lat = y)

# 5. Calculate Moving-Window Shoreline Angles --------------------------------
# Standardizing your column names for get_coast_angle calculation
coast_points_filtered$coast_angle_rad <- sapply(
  1:nrow(coast_points_filtered), 
  function(i) { get_coast_angle(i, coast_points_filtered, window_size = 15) }
)

# 6. Data Visualization ------------------------------------------------------
# Visual check: The map points now align to a rigid, neat 1/10 degree spacing
ggplot(coast_points_filtered, aes(x = lon, y = lat, color = coast_angle_rad)) + 
  geom_point(size = 1.2, shape = 15) + # Uses square points to highlight grid alignment
  scale_color_viridis_c(option = "plasma", name = "Angle (rad)") + 
  theme_minimal() + 
  labs(
    title = "Coastal Orientations (Fixed 1/10° Resolution)",
    subtitle = "Points extracted directly from cell grid centers",
    x = "Longitude", y = "Latitude"
  ) +
  theme(plot.title = element_text(face = "bold", size = 14))

# 7. Export Processed Data ---------------------------------------------------
head(coast_points_filtered)
saveRDS(coast_points_filtered, file = "STATIC/coastal_angles.rds")
