#' Developing a function for 3-dimensional data. This is relevant for MLD, SSH, and SST.
#' using list() 
#' @author Megan L. Feddern
#' 
annual_mean_3d <- function(var_dataset, long_dataset, distance, ds_var_name, df_var_name, filter_dict, stage, precond, print_counts = TRUE) {
  start_time <- Sys.time()
  cat(paste0("---", df_var_name, "---\n"))
  
  # Fetch native Python built-ins for slicing syntax
  py <- import_builtins()
  
  # 1. First, filter the variable called in "var_dataset" by month so you are only
  # interacting with the months for the lifestage of interest
  month_list <- seq(filter_dict[[stage]][["min_month"]], filter_dict[[stage]][["max_month"]])
  var_selected <- var_dataset$sel(time = var_dataset$time$dt$month$isin(month_list))
  
  # Extract the year information from the dataset - this will give you r_years that has
  # a year entry for each month defined in the lifestage dictionary you provide
  py_years <- var_selected$time$dt$year$values
  r_years <- py_to_r(py_years)
  
  # Extract latitudinal information as floating decimals from the lifestage list
  min_lat <- as.numeric(filter_dict[[stage]][["min_lat"]])
  max_lat <- as.numeric(filter_dict[[stage]][["max_lat"]])
  
  # 2. Now that we have the lat, we can slice the latitudinal spatial boundaries 
  # This isolates only the latitude window for the months of interest (above) into 'var_array'
  # so that the dataset is the "smallest" possible which helps with computation time
  cat("Slicing spatial boundaries...\n")
  var_selected <- var_selected$sel(lat = py$slice(min_lat, max_lat))
  
  # 3. since we are doing 3-d data we will now do the same thing with depth
  # this is the primary difference between the 3d and 2d functions
  min_depth <- as.numeric(filter_dict[[stage]][["min_depth"]])
  max_depth <- as.numeric(filter_dict[[stage]][["max_depth"]])
  cat(paste0("Slicing depth dimension (", depth_name, "): ", min_depth, " to ", max_depth, "\n"))
  
  var_selected <- var_selected$sel(z_l=py$slice(min_depth, max_depth))
  var_array <- var_selected[[ds_var_name]]
  
  # 4. We need the sliced Earthmover data to be able to be aligned with our masks
  # this will just dynamically rename "latitude" and "longitude" to "lat" "long" in the static datasets
  # which are defined as "long_dataset"
  # this step is written to be flexible in case the user uses masks with different labels
  # or a lifestage dictionary with different labels
  
  cat("Standardizing mask coordinate names...\n")
  dim_keys <- py_to_r(py$list(long_dataset$dims))
  rename_list <- list()
  if ("latitude" %in% dim_keys) rename_list$latitude <- "lat"
  if ("longitude" %in% dim_keys) rename_list$longitude <- "lon"
  if (!is.null(depth_name) && "deptho" %in% dim_keys) rename_list$deptho <- depth_name
  
  if (length(rename_list) > 0) {
    long_dataset <- long_dataset$rename(rename_list)
  }
  
  # now we have only the data we need AND appropriate naming convention for the mask dataset
  # to be able to apply a mask for longitudinal extent. This can be done 2 ways:
  # 1) applying a mask of bottom depth using the static mom6 bathymetry file 
  # 2) applying a distance from shore maske using the static mom6 file
  
  # 5. Use Xarray's interp_like to ensure alignment
  cat("Interpolating mask grid onto variable layout...\n")
  # this chunk does 2 things, ensures the mask data is aligned with the with the variable data
  # using nearest neighbor, and slices the masking file according to the parameters defined in 
  # the life history dictionary 
  
  # if distance is true, we take the dataset in "long_dataset" and filter it useing the 
  # 'min_dist' and "max_dist" definitions in the lifestage dictionary
  if (distance) {
    mask_sliced <- long_dataset[["dist"]]$interp_like(var_array, method = "nearest")
    py_mask <- (mask_sliced$lat >= min_lat) & (mask_sliced$lat <= max_lat) & 
      (mask_sliced >= filter_dict[[stage]][['min_dist']]) & 
      (mask_sliced <= filter_dict[[stage]][['max_dist']]) & (!mask_sliced$isnull())
  } else {#if distance is 'FALSE' then we filter using min_iso and max_iso
    mask_sliced <- long_dataset[["deptho"]]$interp_like(var_array, method = "nearest")
    py_mask <- (mask_sliced$lat >= min_lat) & (mask_sliced$lat <= max_lat) & 
      (mask_sliced >= filter_dict[[stage]][['min_iso']]) & 
      (mask_sliced <= filter_dict[[stage]][['max_iso']]) & (!mask_sliced$isnull())
  }
  
  # 6. Apply the mask layer
  cat("Applying mask layer...\n")
  var_masked <- var_array$where(py_mask)
  
  # This is just a check: Identify remaining spatial dimensions to summarise means over
  all_dims <- py_to_r(py$list(var_masked$dims))
  spatial_dims_to_drop <- all_dims[all_dims != "time"]
  
  #spatial_dims_to_drop are the spatial dimensions to take a mean over. It should
  # just be "lat" and "long"
  
  # 7. Taking the spatial mean 
  # this takes the mean over the defined latitude and the longitudinal extent
  # filtered by the mask.
  cat("Collapsing spatial/depth dimensions natively...\n")
  spatial_mean_series <- var_masked$mean(dim = spatial_dims_to_drop, skipna = TRUE)
  
  # 8. Pull the monthly means
  # you now have monthly means across the defined lat long. We need to pull them into
  # R and align them with the correct time
  cat("Extracting values safe from memory overhead...\n")
  raw_values <- py_to_r(spatial_mean_series$values)
  
  # Manually construct R data frame with the extracted r_years
  temp_df <- data.frame(
    year = r_years,
    value = raw_values
  )
  
# 9. Aggregate annual time series using r
  cat("Aggregating annual time series in R...\n")
  ds_df <- aggregate(value ~ year, data = temp_df, FUN = mean, na.rm = TRUE)
  
  # If preconditioning is TRUE, add a 1-year lag to the timeline data index
  if (precond) {
    cat("Applying 1-year environmental precondition lag (Years + 1)...\n")
    # In a default pandas to_dataframe() conversion with groupby('time.year'),
    # the years are stored as the index (rownames) of the R data frame.
     current_years <- as.numeric(rownames(ds_df))
    lagged_years <- current_years + 1
    rownames(ds_df) <- as.character(lagged_years)
  }
  
  # Format column and index headers
  rownames(ds_df) <- ds_df$year
  ds_df$year <- NULL
  names(ds_df) <- df_var_name
  
  # Explicitly remove what was generated with this function
  rm(temp_df, raw_values, r_years)
  gc()
  
  # Print script runtime metrics
  end_time <- Sys.time()
  runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
  cat(paste0("Runtime: ", round(runtime_minutes, 2), " minutes\n"))
  
  return(ds_df)
}
