#' Developing a function for 3-dimensional data that calculates transport along-shore
#' and across shelf simultaneously
#' using list() 
#' @author Megan L. Feddern
#' 

transport_3d <- function(u_dataset, v_dataset, long_dataset, coast_ds, filter_dict, stage, precond, print_counts = TRUE) {
  start_time <- Sys.time()
  cat("--- Calculating 3D Coastal Transport (Alongshore & Cross-shore) ---\n")
  
  # Fetch native Python built-ins for slicing syntax
  py <- import_builtins()
  np <- import("numpy") 
  
  # =========================================================================
  # 1. TEMPORAL FILTERING
  # =========================================================================
  #here we are just going to slice the months that we need first
  month_list <- seq(filter_dict[[stage]][["min_month"]], filter_dict[[stage]][["max_month"]])
  
  u_selected <- u_dataset[["uo_rotate"]]$sel(time = u_dataset$time$dt$month$isin(month_list))
  v_selected <- v_dataset[["vo_rotate"]]$sel(time = v_dataset$time$dt$month$isin(month_list))
  
  py_years <- u_selected$time$dt$year$values
  r_years <- py_to_r(py_years)
  
  # =========================================================================
  # 2. LATITUDINAL SLICING
  # =========================================================================
  # similarly we are going to slice the north / south bounds to decrease the 
  # dataset size we are working with
  min_lat <- as.numeric(filter_dict[[stage]][["min_lat"]])
  max_lat <- as.numeric(filter_dict[[stage]][["max_lat"]])
  
  cat("Slicing spatial boundaries for U and V components...\n")
  u_selected <- u_selected$sel(lat = py$slice(min_lat, max_lat))
  v_selected <- v_selected$sel(lat = py$slice(min_lat, max_lat))
  
  # Detect the depth/vertical coordinate dimension name
  # this is probably overkill but its built in flexibility for datasets that may have
  # different lames for depth variables. here it is "z_l" for MOM6 which is the 
  # center point of the depth grid
  
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
  # 3. COORDINATE STANDARDIZATION & DATAFRAME ANGLE MERGING
  # =========================================================================
  # In many global ocean models, U and V vectors are calculated on slightly staggered 
  # spatial grids (for MOM6 its the Arakawa C-grids). This section uses linear interpolation
  # to shift the V coordinate matrix points and the bathymetry masks so that they map 
  # exactly onto the coordinate points of the U grid matrix. If you are curious for more info see
  # Kelly Kearneys excellent MOM6 grid explainer for more details 
  # https://beringnpz.github.io/ak_regional_mom6/cefi_dataset.html
  
  cat("Standardizing coordinate names and mapping latitude-specific angles...\n")
  dim_keys <- py_to_r(py$list(long_dataset$dims))
  rename_list <- list()
  if ("latitude" %in% dim_keys) rename_list$latitude <- "lat"
  if ("longitude" %in% dim_keys) rename_list$longitude <- "lon"
  if (!is.null(depth_name) && "deptho" %in% dim_keys) rename_list$deptho <- depth_name
  
  if (length(rename_list) > 0) {
    long_dataset <- long_dataset$rename(rename_list)
  }
  
  # Interpolate velocity matrices and masks onto matching grid structures using
  # the U grid matrix. Velocity is linearly interpolated and the depth mask
  # is nearest neighbor
  
  v_array_aligned <- v_selected$interp_like(u_selected, method = "linear")
  mask_sliced     <- long_dataset[["deptho"]]$interp_like(u_selected, method = "nearest")
  
  # 3a. Extract the raw angle and latitude coordinates
  raw_angles_deg <- py_to_r(coast_ds[["coast_angle"]]$values)
  coast_lats     <- py_to_r(coast_ds[["lat"]]$values)
  
  # 3b. Convert degrees to radians for use with np$sin/np$cos
  # (Since Mike's MATLAB used 'sind/cosd', but Python's NumPy works exclusively in radians)
  angles_rad <- raw_angles_deg * (pi / 180)
  
  # 3c. Build a 1D xarray using standard coordinate names ('lat')
  cangle_1d <- xr$DataArray(
    data = angles_rad,
    dims = list("lat"),
    coords = list(lat = coast_lats)
  )
  
  # 3d. Interpolate the 1D latitude angle vector exactly to the U grid
  # xarray will merge this across time, longitude, and depth dimensions
  cangle_sliced <- cangle_1d$interp(lat = u_selected$lat, method = "linear")
  
  
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
  # 5. MATHEMATICAL ORIENTATION RESOLUTION & DYNAMIC THICKNESS CALCULATION
  # =========================================================================
  cat("Vector rotation and dynamic cell thickness calculation...\n")
  
  # 5a. Standard coordinate rotation equations
  ua_array <- (v_array_aligned * np$sin(cangle_sliced)) - (u_selected * np$cos(cangle_sliced))
  uc_array <- (v_array_aligned * np$cos(cangle_sliced)) + (u_selected * np$sin(cangle_sliced))
  
  #5b. dynamically reconstructing the depth interfaces z_i based on the depth midpoints z_l starting at the surface 0. use dz_val output to check that this
  # correctly aligns with MOM6 depth bins and isn't bissecting expected thicknesses (confirmation: it does)
  if ("z_l" %in% names(u_dataset)) {
    zl_vals <- py_to_r(u_dataset[["z_l"]]$values)
    nz <- length(zl_vals)
    
    # 1. Initialize the interface array (length nz + 1)
    zi_vals <- numeric(nz + 1)
    
    # 2. Assume the surface interface starts at 0 
    # (Works whether your grid goes positive or negative, as long as it starts at the surface)
    zi_vals[1] <- 0
    
    # 3. Progressively solve for each subsequent interface: z_i[k+1] = 2*z_l[k] - z_i[k]
    for (k in 1:nz) {
      zi_vals[k+1] <- 2 * zl_vals[k] - zi_vals[k]
    }
    
    # 4. The thickness of each layer is the absolute difference between its interfaces
    dz_vals <- abs(diff(zi_vals))
    
    # Construct a 1D xarray DataArray mapped to the original unsliced depth axis
    dz_full <- xr$DataArray(
      data = dz_vals,
      dims = list(depth_name),
      coords = setNames(list(u_dataset[[depth_name]]), depth_name)
    )
  }
  
  
  
  # Match the sliced depth footprint of u_selected matrix
  depth_selector <- unname(list(u_selected[[depth_name]])) |> setNames(depth_name)
  dz_sliced <- dz_full$sel(depth_selector)
  
  # 5c. Apply the spatial and depth bathymetry mask to both velocities and thicknesses
  ua_masked <- ua_array$where(py_mask)
  uc_masked <- uc_array$where(py_mask)
  dz_masked <- dz_sliced$where(py_mask) # Truncates cleanly at the seafloor bottom depth
  
  # =========================================================================
  # 6. THICKNESS-WEIGHTED VERTICAL COLLAPSE & SPATIAL MEAN
  # =========================================================================
  cat("Collapsing water column via high-performance thickness weighting...\n")
  
  # 6a. Compute integrated transport over the water column (Velocity * Layer Thickness)
  ua_integrated <- ua_masked * dz_masked
  uc_integrated <- uc_masked * dz_masked
  
  # 6b. Sum integrated transports and divide by the total water column thickness
  ua_weighted_vertical_mean <- ua_integrated$sum(dim = depth_name, skipna = TRUE) / dz_masked$sum(dim = depth_name, skipna = TRUE)
  uc_weighted_vertical_mean <- uc_integrated$sum(dim = depth_name, skipna = TRUE) / dz_masked$sum(dim = depth_name, skipna = TRUE)
  
  # 6c. Average the resulting 2D surface over the remaining horizontal spatial dimensions
  horizontal_dims <- py_to_r(py$list(ua_weighted_vertical_mean$dims))
  spatial_dims_to_drop <- horizontal_dims[horizontal_dims != "time"]
  
  cat("Aggregating means ...\n")
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
  
  # Dynamically rename the columns based on the current active 'stage' string
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
