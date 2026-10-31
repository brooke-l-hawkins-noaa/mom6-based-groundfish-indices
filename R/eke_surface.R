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
    lat = py$slice(min_lat - 0.5, max_lat + 0.5)
  )
  
  r_years <- py_to_r(u_filtered$time$dt$year$values)
  r_months <- py_to_r(u_filtered$time$dt$month$values)
  
  # =========================================================================
  # 2. DYNAMIC SHALLOWEST LAYER SELECTION
  # =========================================================================
  depth_name <- NULL
  available_dims <- py_to_r(py$list(u_filtered$dims))
  for (d in c("depth", "lev", "plev", "vertical", "z_l")) {
    if (d %in% available_dims) { depth_name <- d; break }
  }
  
  if (!is.null(depth_name)) {
    nz <- length(py_to_r(u_filtered[[depth_name]]$values))
    u_surface <- u_filtered$isel(setNames(list(0L), depth_name))
    v_surface <- v_filtered$isel(setNames(list(0L), depth_name))
    if (nz > 1) {
      for (k in 2:nz) {
        layer_idx <- k - 1L
        u_surface <- u_surface$combine_first(u_filtered$isel(setNames(list(layer_idx), depth_name)))
        v_surface <- v_surface$combine_first(v_filtered$isel(setNames(list(layer_idx), depth_name)))
      }
    }
  } else {
    u_surface <- u_filtered; v_surface <- v_filtered
  }
  
  # =========================================================================
  # 3. ALIGN STAGGERED VELOCITY GRIDS LAZILY
  # =========================================================================
  cat("Aligning staggered velocity grids lazily...\n")
  v_array_aligned <- v_surface$interp_like(u_surface, method = "nearest")
  if (as.logical(v_array_aligned$isnull()$all()$values)) {
    v_array_aligned <- u_surface$copy(data = v_surface$values, deep = TRUE)
  }
  
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
    mask_sliced_v <- long_dataset[["dist"]]$interp_like(v_array_aligned, method = "nearest")
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
    mask_sliced_v <- long_dataset[["deptho"]]$interp_like(v_array_aligned, method = "nearest")
    py_mask_v <- (mask_sliced_v$lat >= min_lat) & (mask_sliced_v$lat <= max_lat) &
      (mask_sliced_v >= filter_dict[[stage]][['min_iso']]) &
      (mask_sliced_v <= filter_dict[[stage]][['max_iso']]) &
      (mask_sliced_v$notnull())
    
  } else {
    # Default Longitude/Latitude bounds check
    py_mask_u <- (u_surface$lat >= min_lat) & (u_surface$lat <= max_lat) &
      (u_surface$lon >= filter_dict[[stage]][['min_long']]) &
      (u_surface$lon <= filter_dict[[stage]][['max_long']])
    
    py_mask_v <- (v_array_aligned$lat >= min_lat) & (v_array_aligned$lat <= max_lat) &
      (v_array_aligned$lon >= filter_dict[[stage]][['min_long']]) &
      (v_array_aligned$lon <= filter_dict[[stage]][['max_long']])
  }
  
  # Apply distinct masks completely lazily
  u_masked <- u_surface$where(py_mask_u)
  v_masked <- v_array_aligned$where(py_mask_v)
  
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
  # 6. SPATIAL COLLAPSE WITH ZERO-COUNT SAFEGUARDS & MEMORY LOAD
  # =========================================================================
  cat("Collapsing horizontal spatial dimensions upfront to protect RAM...\n")
  all_remaining_dims <- py_to_r(py$list(eke_spatial$dims))
  spatial_dims_to_drop <- all_remaining_dims[all_remaining_dims != "time"]
  
  # Average EKE across space for each time step
  eke_sum <- eke_spatial$sum(dim = spatial_dims_to_drop, skipna = TRUE)
  eke_count <- eke_spatial$notnull()$sum(dim = spatial_dims_to_drop)
  eke_time_series_lazy <- eke_sum$where(eke_count > 0, np$nan) / eke_count$where(eke_count > 0, 1.0)
  
  cat("Loading collapsed 1D EKE time series vector into R memory...\n")
  # This is the single trigger point where Python executes the graph and loads the data
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
  
  if (precond) {
    cat("Applying 1-year environmental precondition lag (Years + 1)...\n")
    ds_df$year <- ds_df$year + 1
  }
  
  rownames(ds_df) <- as.character(ds_df$year)
  ds_df$year <- NULL
  names(ds_df) <- "EKE_surface"
  
  # Clean up only the variables used in the new Xarray workflow
  rm(temp_df, raw_values, r_years, r_months)
  gc()
  
  end_time <- Sys.time()
  cat(paste0("Calculation complete! Total execution time: ", round(difftime(end_time, start_time, units="secs"), 2), " seconds.\n"))
  return(ds_df)
}
