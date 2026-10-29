

transport_2d <- function(u_dataset, v_dataset, long_dataset, coast_points_df, filter_dict, stage, precond, print_counts = TRUE) {
  start_time <- Sys.time()
  cat("--- Calculating 3D Coastal Transport (Alongshore & Cross-shore) ---\n")
  
  # Fetch native Python built-ins for slicing syntax
  py <- import_builtins()
  np <- import("numpy") 
  
  # =========================================================================
  # 1. TEMPORAL FILTERING
  # =========================================================================
  month_list <- seq(filter_dict[[stage]][["min_month"]], filter_dict[[stage]][["max_month"]])
  
  u_selected <- u_dataset[["uo_rotate"]]$sel(time = u_dataset$time$dt$month$isin(month_list))
  v_selected <- v_dataset[["vo_rotate"]]$sel(time = v_dataset$time$dt$month$isin(month_list))
  
  py_years <- u_selected$time$dt$year$values
  r_years <- py_to_r(py_years)
  
  # =========================================================================
  # 2. SPATIAL GEOGRAPHIC BOUNDARY SLICING
  # =========================================================================
  min_lat <- as.numeric(filter_dict[[stage]][["min_lat"]])
  max_lat <- as.numeric(filter_dict[[stage]][["max_lat"]])
  
  cat("Slicing spatial boundaries for U and V components...\n")
  u_selected <- u_selected$sel(lat = py$slice(min_lat, max_lat))
  v_selected <- v_selected$sel(lat = py$slice(min_lat, max_lat))
  
  # Detect the depth/vertical coordinate dimension name
  depth_name <- NULL
  available_dims <- py_to_r(py$list(u_selected$dims))
  for (d in c("depth", "lev", "plev", "vertical", "z_l")) {
    if (d %in% available_dims) {
      depth_name <- d
      break
    }
  }
  
  if (!is.null(depth_name) && "min_depth" %in% names(filter_dict[[stage]])) {
    min_depth <- as.numeric(filter_dict[[stage]][["min_depth"]])
    max_depth <- as.numeric(filter_dict[[stage]][["max_depth"]])
    cat(paste0("Slicing depth dimension (", depth_name, "): ", min_depth, " to ", max_depth, "\n"))
    
    depth_slice <- unname(list(py$slice(min_depth, max_depth))) |> setNames(depth_name)
    
    u_selected <- u_selected$sel(depth_slice)
    v_selected <- v_selected$sel(depth_slice)
  }
  
  # =========================================================================
  # 3. COORDINATE STANDARDIZATION & DATAFRAME ANGLE BROADCASTING
  # =========================================================================
  # In many global ocean models, U and V vectors are calculated on slightly staggered 
  # spatial grids (for MOM6 its the Arakawa C-grids). This section uses linear interpolation
  # to shift the V coordinate matrix points and the bathymetry masks so that they map 
  # exactly onto the coordinate points of the U grid matrix. If you are curious see
  # Kelly Kearneys excellent MOM6 grid explainer for more details 
  # https://beringnpz.github.io/ak_regional_mom6/cefi_dataset.html
  
  cat("Standardizing coordinate names and mapping dataframe angles...\n")
  dim_keys <- py_to_r(py$list(long_dataset$dims))
  rename_list <- list()
  if ("latitude" %in% dim_keys) rename_list$latitude <- "lat"
  if ("longitude" %in% dim_keys) rename_list$longitude <- "lon"
  if (!is.null(depth_name) && "deptho" %in% dim_keys) rename_list$deptho <- depth_name
  
  if (length(rename_list) > 0) {
    long_dataset <- long_dataset$rename(rename_list)
  }
  
  # Interpolate velocity matrices and masks onto matching grid structures using
  # the U grid matrix
  
  v_array_aligned <- v_selected$interp_like(u_selected, method = "linear")
  mask_sliced     <- long_dataset[["deptho"]]$interp_like(u_selected, method = "nearest")
  
  # Extract the mean angle value from your calculated flat R dataframe
  # and broadcast it uniformly across the active 3D array matrix dimensions
  # because everything is now aligned with the U grid, you can just use U here
  # Computes the mathematical average of the orientation angle values (in radians) 
  # stored inside coast_points_df table, ignoring missing data entries.
  mean_angle_rad <- mean(coast_points_df$coast_angle_rad, na.rm = TRUE)
  cangle_sliced <- u_selected * 0 + mean_angle_rad
  
  # =========================================================================
  # 4. MASK GENERATION
  # =========================================================================
  # we have a single mask generationg here, because again, we are working with everything
  # translated onto the U grid. Generates an object containing logical TRUE/FALSE 
  # states for every point on the grid. A coordinate point yields TRUE only if 
  # it matches your latitude box limits, falls within the depth constraints 
  # (min_iso to max_iso), and does not represent an empty/land data void.
  
  py_mask <- (mask_sliced$lat >= min_lat) & (mask_sliced$lat <= max_lat) & 
    (mask_sliced >= filter_dict[[stage]][['min_iso']]) & 
    (mask_sliced <= filter_dict[[stage]][['max_iso']]) & (!mask_sliced$isnull())
  
  # =========================================================================
  # 5. MATHEMATICAL ORIENTATION RESOLUTION & THICKNESS EXTRACTION
  # =========================================================================
  cat("Resolving vectors along high-res coastal orientation models...\n")
  
  # 5a. Standard coordinate rotation equations
  ua_array <- (v_array_aligned * np$sin(cangle_sliced)) - (u_selected * np$cos(cangle_sliced))
  uc_array <- (v_array_aligned * np$cos(cangle_sliced)) + (u_selected * np$sin(cangle_sliced))
  
  # 5b. Extract the layer thicknesses (dzh) matching your sliced domain
  # Note: Adjust "dzh" to your specific dataset variable name if it differs (e.g., "thkcello")
  dz_sliced <- long_dataset[["dzh"]]$interp_like(u_selected, method = "nearest")
  
  # 5c. Apply the spatial and depth bathymetry mask to both velocities and thicknesses
  ua_masked <- ua_array$where(py_mask)
  uc_masked <- uc_array$where(py_mask)
  dz_masked <- dz_sliced$where(py_mask) # Truncates layers at the true bottom depth
  
  # =========================================================================
  # 6. THICKNESS-WEIGHTED VERTICAL COLLAPSE & SPATIAL MEAN
  # =========================================================================
  cat("Collapsing water column via high-performance thickness weighting...\n")
  
  # 6a. Compute the integrated transport over the water column (Velocity * Layer Thickness)
  ua_integrated <- ua_masked * dz_masked
  uc_integrated <- uc_masked * dz_masked
  
  # 6b. Sum the integrated transport and divide by total unmasked water column thickness
  # This isolates the vertical dimension name dynamically detected in Section 2
  ua_weighted_vertical_mean <- ua_integrated$sum(dim = depth_name, skipna = TRUE) / dz_masked$sum(dim = depth_name, skipna = TRUE)
  uc_weighted_vertical_mean <- uc_integrated$sum(dim = depth_name, skipna = TRUE) / dz_masked$sum(dim = depth_name, skipna = TRUE)
  
  # 6c. Average the resulting 2D surface over the remaining horizontal spatial dimensions
  horizontal_dims <- py_to_r(py$list(ua_weighted_vertical_mean$dims))
  spatial_dims_to_drop <- horizontal_dims[horizontal_dims != "time"]
  
  cat("Collapsing final horizontal dimensions inside Python engine...\n")
  ua_series <- ua_weighted_vertical_mean$mean(dim = spatial_dims_to_drop, skipna = TRUE)
  uc_series <- uc_weighted_vertical_mean$mean(dim = spatial_dims_to_drop, skipna = TRUE)
  
  # =========================================================================
  # 7. PARSING DATA INTO R DATA FRAMES
  # =========================================================================
  cat("Extracting values free of memory overhead...\n")
  raw_ua <- py_to_r(ua_series$values)
  raw_uc <- py_to_r(uc_series$values)
  
  temp_df <- data.frame(
    year = r_years,
    u_along = raw_ua,
    u_cross = raw_uc
  )
  
  # =========================================================================
  # 8. R AGGREGATION & LAGGING WRAPPERS
  # =========================================================================
  cat("Aggregating annual timeline numbers in R...\n")
  ds_df_along <- aggregate(u_along ~ year, data = temp_df, FUN = mean, na.rm = TRUE)
  ds_df_cross <- aggregate(u_cross ~ year, data = temp_df, FUN = mean, na.rm = TRUE)
  
  ds_df <- merge(ds_df_along, ds_df_cross, by = "year")
  
  if (precond) {
    cat("Applying 1-year environmental precondition lag (Years + 1)...\n")
    ds_df$year <- ds_df$year + 1
  }
  
  rownames(ds_df) <- ds_df$year
  ds_df$year <- NULL
  
  rm(temp_df, raw_ua, raw_uc, r_years)
  gc()
  
  end_time <- Sys.time()
  runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
  cat(paste0("Runtime: ", round(runtime_minutes, 2), " minutes\n"))
  
  return(ds_df)
}

