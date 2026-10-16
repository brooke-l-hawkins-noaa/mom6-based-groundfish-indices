

transport_3d <- function(u_dataset, v_dataset, long_dataset, coast_points_df, filter_dict, stage, precond, print_counts = TRUE) {
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
  cat("Standardizing coordinate names and mapping dataframe angles...\n")
  dim_keys <- py_to_r(py$list(long_dataset$dims))
  rename_list <- list()
  if ("latitude" %in% dim_keys) rename_list$latitude <- "lat"
  if ("longitude" %in% dim_keys) rename_list$longitude <- "lon"
  if (!is.null(depth_name) && "deptho" %in% dim_keys) rename_list$deptho <- depth_name
  
  if (length(rename_list) > 0) {
    long_dataset <- long_dataset$rename(rename_list)
  }
  
  # Interpolate velocity matrices and masks onto matching grid structures
  v_array_aligned <- v_selected$interp_like(u_selected, method = "linear")
  mask_sliced     <- long_dataset[["deptho"]]$interp_like(u_selected, method = "nearest")
  
  # Extract the mean angle value from your calculated flat R dataframe
  # and broadcast it uniformly across the active 3D array matrix dimensions
  mean_angle_rad <- mean(coast_points_df$coast_angle_rad, na.rm = TRUE)
  cangle_sliced <- u_selected * 0 + mean_angle_rad
  
  # =========================================================================
  # 4. CROSS-SHORE MASK GENERATION
  # =========================================================================
  py_mask <- (mask_sliced$lat >= min_lat) & (mask_sliced$lat <= max_lat) & 
    (mask_sliced >= filter_dict[[stage]][['min_iso']]) & 
    (mask_sliced <= filter_dict[[stage]][['max_iso']]) & (!mask_sliced$isnull())
  
  # =========================================================================
  # 5. MATHEMATICAL ORIENTATION RESOLUTION (RADIAN EQUATIONS)
  # =========================================================================
  cat("Resolving vectors along high-res coastal orientation models...\n")
  ua_array <- v_array_aligned * np$sin(cangle_sliced) - u_selected * np$cos(cangle_sliced)
  uc_array <- v_array_aligned * np$cos(cangle_sliced) + u_selected * np$sin(cangle_sliced)
  
  ua_masked <- ua_array$where(py_mask)
  uc_masked <- uc_array$where(py_mask)
  
  # =========================================================================
  # 6. NATIVE DIMENSION COLLAPSING (XARRAY ENGINE)
  # =========================================================================
  all_dims <- py_to_r(py$list(ua_masked$dims))
  spatial_dims_to_drop <- all_dims[all_dims != "time"]
  
  cat("Collapsing spatial/vertical dimensions inside Python engine...\n")
  ua_series <- ua_masked$mean(dim = spatial_dims_to_drop, skipna = TRUE)
  uc_series <- uc_masked$mean(dim = spatial_dims_to_drop, skipna = TRUE)
  
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

