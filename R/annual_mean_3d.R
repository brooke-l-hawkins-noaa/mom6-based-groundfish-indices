#' Developing a function for 3-dimensional data. This is relevant for MLD, SSH, and SST.
#' using list() 
#' @author Megan L. Feddern
#' 
annual_mean_3d <- function(var_dataset, long_dataset, distance, ds_var_name, df_var_name, filter_dict, stage, precond, print_counts = TRUE) {
  start_time <- Sys.time()
  cat(paste0("---", df_var_name, "---\n"))
  
  # Fetch native Python built-ins for slicing syntax
  py <- import_builtins()
  
  # 1. Temporal filter on the main array (1D operation)
  month_list <- seq(filter_dict[[stage]][["min_month"]], filter_dict[[stage]][["max_month"]])
  var_selected <- var_dataset$sel(time = var_dataset$time$dt$month$isin(month_list))
  
  # Extract year timelines safely
  py_years <- var_selected$time$dt$year$values
  r_years <- py_to_r(py_years)
  
  # Extract geographical coordinates cleanly as floating decimals
  min_lat <- as.numeric(filter_dict[[stage]][["min_lat"]])
  max_lat <- as.numeric(filter_dict[[stage]][["max_lat"]])
  
  # 2. Slice spatial coordinate boundaries (Lat and Depth if present)
  cat("Slicing spatial boundaries...\n")
  var_selected <- var_selected$sel(lat = py$slice(min_lat, max_lat))
  
  # Check for depth coordinates in filter_dict and slice if available
  # (Handles common depth names like 'depth', 'lev', 'plev')
  depth_name <- NULL
  available_dims <- py_to_r(py$list(var_selected$dims))
  for (d in c("depth", "lev", "plev", "vertical")) {
    if (d %in% available_dims) {
      depth_name <- d
      break
    }
  }
  
  if (!is.null(depth_name) && "min_depth" %in% names(filter_dict[[stage]])) {
    min_depth <- as.numeric(filter_dict[[stage]][["min_depth"]])
    max_depth <- as.numeric(filter_dict[[stage]][["max_depth"]])
    cat(paste0("Slicing depth dimension (", depth_name, "): ", min_depth, " to ", max_depth, "\n"))
    var_selected <- var_selected$sel(unname(list(py$slice(min_depth, max_depth))) |> setNames(depth_name))
  }
  
  var_array <- var_selected[[ds_var_name]]
  
  # 3. Dynamic Coordinate Name Standardization for long_dataset
  cat("Standardizing mask coordinate names...\n")
  dim_keys <- py_to_r(py$list(long_dataset$dims))
  rename_list <- list()
  if ("latitude" %in% dim_keys) rename_list$latitude <- "lat"
  if ("longitude" %in% dim_keys) rename_list$longitude <- "lon"
  if (!is.null(depth_name) && "deptho" %in% dim_keys) rename_list$deptho <- depth_name
  
  if (length(rename_list) > 0) {
    long_dataset <- long_dataset$rename(rename_list)
  }
  
  # 4. Use Xarray's interp_like to ensure alignment
  cat("Interpolating mask grid onto variable layout...\n")
  if (distance) {
    mask_sliced <- long_dataset[["dist"]]$interp_like(var_array, method = "nearest")
    py_mask <- (mask_sliced$lat >= min_lat) & (mask_sliced$lat <= max_lat) & 
      (mask_sliced >= filter_dict[[stage]][['min_dist']]) & 
      (mask_sliced <= filter_dict[[stage]][['max_dist']]) & (!mask_sliced$isnull())
  } else {
    mask_sliced <- long_dataset[["deptho"]]$interp_like(var_array, method = "nearest")
    py_mask <- (mask_sliced$lat >= min_lat) & (mask_sliced$lat <= max_lat) & 
      (mask_sliced >= filter_dict[[stage]][['min_iso']]) & 
      (mask_sliced <= filter_dict[[stage]][['max_iso']]) & (!mask_sliced$isnull())
  }
  
  # 5. Apply the mask layer
  cat("Applying mask layer...\n")
  var_masked <- var_array$where(py_mask)
  
  # Identify remaining spatial dimensions to collapse (lat, lon, and depth)
  all_dims <- py_to_r(py$list(var_masked$dims))
  spatial_dims_to_drop <- all_dims[all_dims != "time"]
  
  # 6. Collapse spatial and depth dimensions natively inside Python
  cat("Collapsing spatial/depth dimensions natively...\n")
  spatial_mean_series <- var_masked$mean(dim = spatial_dims_to_drop, skipna = TRUE)
  
  # 7. Pull the aggregated numeric timelines into R environment
  cat("Extracting values safe from memory overhead...\n")
  raw_values <- py_to_r(spatial_mean_series$values)
  
  # Manually construct a lightweight R data frame
  temp_df <- data.frame(
    year = r_years,
    value = raw_values
  )
  
  # 8. Aggregate annual timeline averages using R's engine
  cat("Aggregating annual timeline numbers in R...\n")
  ds_df <- aggregate(value ~ year, data = temp_df, FUN = mean, na.rm = TRUE)
  
  # --- ADDED PRECONDITIONING LAG LOGIC ---
  if (precond) {
    cat("Applying 1-year environmental precondition lag (Years + 1)...\n")
    current_years <- as.numeric(rownames(ds_df))
    lagged_years <- current_years + 1
    rownames(ds_df) <- as.character(lagged_years)
  }
  
  # Format column and index headers
  rownames(ds_df) <- ds_df$year
  ds_df$year <- NULL
  names(ds_df) <- df_var_name
  
  # Explicitly trigger system garbage collection
  rm(temp_df, raw_values, r_years)
  gc()
  
  # Print script runtime metrics
  end_time <- Sys.time()
  runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
  cat(paste0("Runtime: ", round(runtime_minutes, 2), " minutes\n"))
  
  return(ds_df)
}
