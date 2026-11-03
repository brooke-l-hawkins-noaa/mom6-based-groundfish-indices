#' Developing a function for 3-dimensional data that calculates transport along-shore
#' and across shelf simultaneously
#' using list()
#' @author Megan L. Feddern

transport_3d <- function(u_dataset, v_dataset, long_dataset, coast_ds, filter_dict, stage, precond, print_counts = TRUE) {
  start_time <- Sys.time()
  cat("--- Calculating 3D Coastal Transport (Alongshore & Cross-shelf) ---\n")
  
  # Fetch native Python built-ins for slicing syntax
  py <- import_builtins()
  np <- import("numpy")
  
  # =========================================================================
  # 1. SLICING TO THE DATA WE WANT
  # =========================================================================
  
  # 1A. temporal slicing of months
  month_list <- seq(filter_dict[[stage]][["min_month"]], filter_dict[[stage]][["max_month"]])
  
  u_selected <- u_dataset[["uo_rotate"]]$sel(time = u_dataset$time$dt$month$isin(month_list))
  v_selected <- v_dataset[["vo_rotate"]]$sel(time = v_dataset$time$dt$month$isin(month_list))
  
  py_years <- u_selected$time$dt$year$values
  r_years <- py_to_r(py_years)
  
  # 1b. latitudinal bounds
  min_lat <- as.numeric(filter_dict[[stage]][["min_lat"]])
  max_lat <- as.numeric(filter_dict[[stage]][["max_lat"]])
  
  cat("Slicing spatial boundaries for U and V components...\n")
  u_selected <- u_selected$sel(lat = py$slice(min_lat, max_lat))
  v_selected <- v_selected$sel(lat = py$slice(min_lat, max_lat))
  
  # 1c. depths bounds
  depth_name <- "z_l"
  min_depth <- as.numeric(filter_dict[[stage]][["min_depth"]])
  max_depth <- as.numeric(filter_dict[[stage]][["max_depth"]])
  
  
  depth_slice <- unname(list(py$slice(min_depth, max_depth))) |> setNames(depth_name)
  u_selected <- u_selected$sel(depth_slice)
  v_selected <- v_selected$sel(depth_slice)
  
  # from here on out we are working with u and v data that already has the lat, depth, 
  # and month domains that we want
  # =========================================================================
  # 2. ALIGNMENT VIA C-GRID INDEX SHIFTING (PSI GRID TRANSFORMATION)
  # =========================================================================
  cat("Transforming components to shared PSI grid via C-grid index shifts...\n")
  
  # 2a. Replicate MATLAB index-position averaging: u_psi = (u(:,1:end-1,:)+u(:,2:end,:))/2;
  # We extract underlying array values to prevent xarray from misaligning differing coordinate tags (i.e, stripping away the coordinate tracking of xarray)
  u_part1 <- u_selected$isel(lat = py$slice(0L, -1L))
  u_part2 <- u_selected$isel(lat = py$slice(1L, NULL))
  u_psi_values <- (u_part1$values + u_part2$values) / 2
  u_array_aligned <- u_part1$copy(data = u_psi_values)
  
  # 2b. Replicate MATLAB index-position averaging: v_psi = (v(1:end-1,:,:)+v(2:end,:,:))/2;
  v_part1 <- v_selected$isel(lon = py$slice(0L, -1L))
  v_part2 <- v_selected$isel(lon = py$slice(1L, NULL))
  v_psi_values <- (v_part1$values + v_part2$values) / 2
  v_array_aligned <- v_part1$copy(data = v_psi_values)
  
  # 2c. align the lontitudinal dataset with the psi grid 
  mask_sliced <- long_dataset[["deptho"]]$interp_like(u_array_aligned, method = "nearest")
  
  # 2e. Extract the raw angle and latitude coordinates
  raw_angles_deg <- py_to_r(coast_ds[["coast_angle"]]$values)
  coast_lats <- py_to_r(coast_ds[["lat"]]$values)
  
  # 2f. Convert degrees to radians for use with np$sin/np$cos
  # (Since Mike's MATLAB used 'sind/cosd', but Python's NumPy works exclusively in radians)
  angles_rad <- raw_angles_deg * (pi / 180)
  
  # 2g. Interpolate the 1D latitude angle vector to the U grid
  target_lats <- py_to_r(u_array_aligned$lat$values)
  valid_idx <- !is.na(coast_lats) & !is.na(angles_rad)
  clean_lats <- coast_lats[valid_idx]
  clean_angles <- angles_rad[valid_idx]
  
  if (length(unique(clean_lats)) >= 2) {
    interp_angles <- approx(x = clean_lats, y = clean_angles, xout = target_lats, method = "linear", rule = 2)$y
  } else if (length(clean_angles) > 0) {
    interp_angles <- rep(clean_angles, length(target_lats))
  } else {
    interp_angles <- rep(0, length(target_lats))
  }
  
  # We enforce matching sizes by assigning data evaluated against u_array_aligned$lat directly to u_array_aligned$lat coords
  cangle_sliced <- xr$DataArray(
    data = interp_angles,
    dims = list("lat"),
    coords = list(lat = u_array_aligned$lat)
  )
  
  # =========================================================================
  # 3. MASK GENERATION FOR LONGITUDINAL EXTENT
  # =========================================================================
  py_mask <- (mask_sliced$lat >= min_lat) & (mask_sliced$lat <= max_lat) & 
    (mask_sliced >= filter_dict[[stage]][['min_iso']]) & 
    (mask_sliced <= filter_dict[[stage]][['max_iso']]) & 
    (!mask_sliced$isnull())
  
  # =========================================================================
  # 4. VECTOR ROTATION & DYNAMIC THICKNESS CALCULATION
  # =========================================================================
  cat("Vector rotation and dynamic cell thickness calculation...\n")
  
  # 4a. Calculate alongshore and cross-shore wind stress
  
  ua_array <- (v_array_aligned * np$sin(cangle_sliced)) - (u_array_aligned * np$cos(cangle_sliced))
  uc_array <- (v_array_aligned * np$cos(cangle_sliced)) + (u_array_aligned * np$sin(cangle_sliced))
  
  # 4b. calculate depth profiles
  
  if ("z_l" %in% names(u_dataset)) {
    zl_vals <- py_to_r(u_dataset[["z_l"]]$values)
    nz <- length(zl_vals)
    zi_vals <- numeric(nz + 1)
    zi_vals[1] <- 0
    for (k in 1:nz) {
      zi_vals[k+1] <- 2 * zl_vals[k] - zi_vals[k]
    }
    dz_vals <- abs(diff(zi_vals))
    
    dz_full <- xr$DataArray(
      data = dz_vals,
      dims = list(depth_name),
      coords = setNames(list(u_dataset[[depth_name]]), depth_name)
    )
  }
  
  depth_selector <- unname(list(u_array_aligned[[depth_name]])) |> setNames(depth_name)
  dz_sliced <- dz_full$sel(depth_selector)
  
  ua_masked <- ua_array$where(py_mask)
  uc_masked <- uc_array$where(py_mask)
  dz_masked <- dz_sliced$where(py_mask) #the values in this are the depth bins
  
  # =========================================================================
  # 5. THICKNESS-INTEGRATED TRANSPORT & SPATIAL MEAN
  # =========================================================================
  cat("Collapsing water column via high-performance thickness weighting...\n")
  # 5a. integrating across depth (Interpolate to z_int in Mike's code)
  ua_integrated <- ua_masked * dz_masked
  uc_integrated <- uc_masked * dz_masked
  
  # 5b. calculating the mean of the vertical profile
  ua_weighted_vertical_mean <- ua_integrated$sum(dim = depth_name, skipna = TRUE) / dz_masked$sum(dim = depth_name, skipna = TRUE)
  uc_weighted_vertical_mean <- uc_integrated$sum(dim = depth_name, skipna = TRUE) / dz_masked$sum(dim = depth_name, skipna = TRUE)
  
  # 5C. calculating means across space 
  
  horizontal_dims <- py_to_r(py$list(ua_weighted_vertical_mean$dims))
  spatial_dims_to_drop <- horizontal_dims[horizontal_dims != "time"]
  cat("Aggregating means across space ...\n")
  ua_series <- ua_weighted_vertical_mean$mean(dim = spatial_dims_to_drop, skipna = TRUE)
  uc_series <- uc_weighted_vertical_mean$mean(dim = spatial_dims_to_drop, skipna = TRUE)
  
  # =========================================================================
  # 6. PARSING DATA INTO R DATA FRAMES
  # =========================================================================
  cat("Extracting values free of memory overhead...\n")
  raw_ua <- py_to_r(ua_series$values)
  raw_uc <- py_to_r(uc_series$values)
  
  # this data frame has monthly means
  temp_df <- data.frame(
    year = r_years,
    u_along = raw_ua,
    u_cross = raw_uc
  )
  
  # =========================================================================    
  # 7. MONTHLY MEANS AND TEMPORAL LAGS
  # =========================================================================
  cat("Aggregating annual timeline numbers in R...\n")
  # 7a. aggregating across months to take a monthly mean
  ds_df_along <- aggregate(u_along ~ year, data = temp_df, FUN = mean, na.rm = TRUE)
  ds_df_cross <- aggregate(u_cross ~ year, data = temp_df, FUN = mean, na.rm = TRUE)
  # 7b. combining u and v in one dataframe to export
  ds_df <- merge(ds_df_along, ds_df_cross, by = "year")
  # 7c. if preconditioning is specified, applying a lag
  if (precond) {
    cat("Applying 1-year environmental precondition lag (Years + 1)...\n")
    ds_df$year <- ds_df$year + 1
  }
  
  rownames(ds_df) <- ds_df$year
  ds_df$year <- NULL
  
  # 7d. Dynamically rename the columns based on the current active 'stage' string
  ast_name <- paste0("AST_", stage)
  cst_name <- paste0("CST_", stage)
  colnames(ds_df) <- c(ast_name, cst_name)
  
  rm(temp_df, raw_ua, raw_uc, r_years)
  gc()
  
  end_time <- Sys.time()
  runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
  cat(paste0("Runtime: ", round(runtime_minutes, 2), " minutes\n"))
  
  return(ds_df)
  
}

