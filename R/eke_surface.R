#' Calculate Surface Eddy Kinetic Energy (EKE) using the first depth layer (0-5m)
#' @author Megan L. Feddern (Adapted)
#' 
surface_eke <- function(u_dataset, v_dataset, long_dataset, mask_type, filter_dict, stage, precond, print_counts = TRUE) {
  start_time <- Sys.time()
  cat("--- Calculating Surface EKE 0-5m (Velocity Anomalies) ---\n")
  
  # Fetch native Python built-ins
  py <- import_builtins()
  np <- import("numpy")
  
  # =========================================================================
  # 1. INITIAL TEMPORAL & LATITUDINAL SLICING UPFRONT
  # =========================================================================
  cat("Slicing spatial and temporal boundaries upfront...\n")
  month_list <- seq(filter_dict[[stage]][["min_month"]], filter_dict[[stage]][["max_month"]])
  min_lat <- as.numeric(filter_dict[[stage]][["min_lat"]])
  max_lat <- as.numeric(filter_dict[[stage]][["max_lat"]])
  
  u_filtered <- u_dataset[["uo_rotate"]]$sel(
    time = u_dataset$time$dt$month$isin(month_list),
    lat = py$slice(min_lat, max_lat)
  )
  v_filtered <- v_dataset[["vo_rotate"]]$sel(
    time = v_dataset$time$dt$month$isin(month_list),
    lat = py$slice(min_lat, max_lat)
  )
  
  r_years <- py_to_r(u_filtered$time$dt$year$values)
  r_months <- py_to_r(u_filtered$time$dt$month$values)
  
  # =========================================================================
  # 2.SHALLOWEST LAYER SELECTION
  # =========================================================================
  # 1. Create the named list tracking your depth parameter axis
  target_depth_slice <- unname(list(c(2.5))) |> setNames(depth_name)
  
  # 2. Convert to a native Python dictionary using reticulate's dict engine
  py_dict_slice <- reticulate::dict(target_depth_slice)
  
  # 3. Apply the slice safely via xarray
  u_surface <- u_filtered$sel(py_dict_slice, method = "nearest")
  v_surface <- v_filtered$sel(py_dict_slice, method = "nearest")

  #=========================================================================
  # 4. BROADCAST 1D DISTANCE METRIC ACROSS THE 2D LONGITUDE MESH (SEPARATE U & V MASKS)
  # =========================================================================
  if (!is.null(mask_type) && mask_type == "distance") {
    # Generate U-grid mask
    mask_sliced_u <- long_dataset[["dist"]]$interp_like(u_surface, method = "nearest")
    py_mask_u <- (mask_sliced_u$lat >= min_lat) & (mask_sliced_u$lat <= max_lat) &
      (mask_sliced_u >= filter_dict[[stage]][['min_dist']]) &
      (mask_sliced_u <= filter_dict[[stage]][['max_dist']]) &
      (mask_sliced_u$notnull())
    
    # Generate V-grid mask
    mask_sliced_v <- long_dataset[["dist"]]$interp_like(v_surface, method = "nearest")
    py_mask_v <- (mask_sliced_v$lat >= min_lat) & (mask_sliced_v$lat <= max_lat) &
      (mask_sliced_v >= filter_dict[[stage]][['min_dist']]) &
      (mask_sliced_v <= filter_dict[[stage]][['max_dist']]) &
      (mask_sliced_v$notnull())
    
  } else if (!is.null(mask_type) && mask_type == "isobath") {
    # Generate U-grid mask
    mask_sliced_u <- long_dataset[["deptho"]]$interp_like(u_surface, method = "nearest")
    py_mask_u <- (mask_sliced_u$lat >= min_lat) & (mask_sliced_u$lat <= max_lat) &
      (mask_sliced_u >= filter_dict[[stage]][['min_iso']]) &
      (mask_sliced_u <= filter_dict[[stage]][['max_iso']]) &
      (mask_sliced_u$notnull())
    
    # Generate V-grid mask
    mask_sliced_v <- long_dataset[["deptho"]]$interp_like(v_surface, method = "nearest")
    py_mask_v <- (mask_sliced_v$lat >= min_lat) & (mask_sliced_v$lat <= max_lat) &
      (mask_sliced_v >= filter_dict[[stage]][['min_iso']]) &
      (mask_sliced_v <= filter_dict[[stage]][['max_iso']]) &
      (mask_sliced_v$notnull())
    
  } else {
    # Default Longitude/Latitude bounds check
    py_mask_u <- (u_surface$lat >= min_lat) & (u_surface$lat <= max_lat) &
      (u_surface$lon >= filter_dict[[stage]][['min_long']]) &
      (u_surface$lon <= filter_dict[[stage]][['max_long']])
    
    py_mask_v <- (v_surface$lat >= min_lat) & (v_surface$lat <= max_lat) &
      (v_surface$lon >= filter_dict[[stage]][['min_long']]) &
      (v_surface$lon <= filter_dict[[stage]][['max_long']])
  }
  
  # Apply distinct masks completely lazily
  u_masked <- u_surface$where(py_mask_u)
  v_masked <- v_surface$where(py_mask_v)
  
  # =========================================================================
  # 5. CALCULATE SPATIAL ANOMALIES IN XARRAY (FULLY LAZY)
  # =========================================================================
  cat("Calculating monthly velocity anomalies across the spatial mesh inside Xarray...\n")
  
  # Group by month to calculate the climatological mean for each grid point
  u_climatology <- u_masked$groupby("time.month")$mean("time", skipna = TRUE)
  v_climatology <- v_masked$groupby("time.month")$mean("time", skipna = TRUE)
  
  # Subtract the monthly climatology from the original masked datasets to get anomalies
  u_anom_spatial <- u_masked$groupby("time.month") - u_climatology
  v_anom_spatial <- v_masked$groupby("time.month") - v_climatology
  
  # Calculate EKE on the full 3D spatial grid (0.5 * (u'^2 + v'^2))
  eke_spatial <- 0.5 * (u_anom_spatial^2 + v_anom_spatial^2)
  
  # =========================================================================
  # 6. SPATIAL MEAN (WITH AREA WEIGHTING)
  # =========================================================================
  cat("Calculating area weights and collapsing spatial dimensions...\n")
  
  # 1. Grab the latitude coordinate from your Python object
  # (Assumes your latitude dimension/coordinate is named "lat" or "latitude")
  lat <- eke_spatial$lat 
  
  # 2. Compute weights: cos(lat in radians)
  # np is assumed to be your imported numpy reference in R (e.g., np <- import("numpy"))
  weights <- np$cos(np$deg2rad(lat))
  
  # 3. Apply weights to the EKE field where data is not null
  # xarray handles broadcasting the 1D latitude weights across the 3D grid automatically
  weighted_eke <- eke_spatial$weighted(weights)
  
  # 4. Determine spatial dimensions to collapse
  all_remaining_dims <- py_to_r(py$list(eke_spatial$dims))
  spatial_dims_to_drop <- all_remaining_dims[all_remaining_dims != "time"]
  
  # 5. Compute the weighted mean across the spatial dimensions
  # xarray's weighted mean automatically accounts for NaNs/null masks in its denominator
  eke_time_series_lazy <- weighted_eke$mean(dim = spatial_dims_to_drop, skipna = TRUE)
  
  cat("Loading collapsed 1D weighted EKE time series vector into R memory...\n")
  # This triggers the computation graph
  raw_values <- as.numeric(py_to_r(eke_time_series_lazy$compute()$values))

  valid_count <- sum(!is.na(raw_values))
  cat(paste0(">>> Diagnostic: Out of ", length(raw_values), " time steps, ", valid_count, " contain valid numbers <<<\n"))
  
  if (valid_count == 0) {
    stop("All calculated EKE values are NaN! The spatial mesh grid masks did not find intersecting water elements.")
  }
  
  temp_df <- data.frame(year = r_years, value = raw_values)
  
  # =========================================================================
  # 7. AGGREGATE ANNUAL TIME SERIES
  # =========================================================================
  cat("Aggregating annual time series in R...\n")
  ds_df <- aggregate(value ~ year, data = temp_df, FUN = mean, na.rm = TRUE)
  
  # If preconditioning is TRUE, add a 1-year lag to the timeline data index
  if (precond==TRUE) {
    cat("Applying 1-year environmental precondition lag (Years + 1)...\n")
    # In a default pandas to_dataframe() conversion with groupby('time.year'),
    # the years are stored as the index (rownames) of the R data frame.
    current_years <- ds_df$year
    lagged_years <- current_years + 1
    ds_df$year <- lagged_years
  }
  # Format column and index headers
  rownames(ds_df) <- ds_df$year
  ds_df$year <- NULL
  names(ds_df) <- "EKE"
  
  eke_name <- paste0("EKE", stage)
  colnames(ds_df) <- c(eke_name)
  
  # Clean up only the variables used in the new Xarray workflow
  rm(temp_df, raw_values, r_years, r_months)
  gc()
  
  end_time <- Sys.time()
  cat(paste0("Calculation complete! Total execution time: ", round(difftime(end_time, start_time, units="secs"), 2), " seconds.\n"))
  return(ds_df)
}
