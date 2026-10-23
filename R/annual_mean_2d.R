#' Developing a function for 2-dimensional data. This is relevant for MLD, SSH, and SST.
#' using list() 
#' @author Megan L. Feddern

annual_mean_2d <- function(var_dataset, long_dataset, distance, ds_var_name, df_var_name, filter_dict, stage, precond, print_counts = TRUE) {
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
  var_selected <- var_selected$sel(lat = py$slice(min_lat, max_lat))
  var_array <- var_selected[[ds_var_name]]
  
  # 3. We need the sliced Earthmover data to be able to be aligned with our masks
  # this will just dynamically rename "latitude" and "longitude" to "lat" "long" in the static datasets
  # which are defined as "long_dataset"
  # this step is written to be flexible in case the user uses masks with different labels
  # or a lifestage dictionary with different labels
  
  cat("Standardizing mask coordinate names...\n")
  dim_keys <- py_to_r(py$list(long_dataset$dims))
  
  rename_list <- list()
  if ("latitude" %in% dim_keys)  rename_list$latitude <- "lat"
  if ("longitude" %in% dim_keys) rename_list$longitude <- "lon"
  
  if (length(rename_list) > 0) {
    long_dataset <- long_dataset$rename(rename_list)
  }
  
  # now we have only the data we need AND appropriate naming convention for the mask dataset
  # to be able to apply a mask for longitudinal extent. This can be done 2 ways:
  # 1) applying a mask of bottom depth using the static mom6 bathymetry file 
  # 2) applying a distance from shore maske using the static mom6 file
  
  # 4. Use Xarray's interp_like to ensure alignment
  cat("Interpolating mask grid onto variable layout...\n")
  # this chunk does 2 things, ensures the mask data is aligned with the with the variable data
  # using nearest neighbor, and slices the masking file according to the parameters defined in 
  # the life history dictionary 
  
  # if distance is true, we take the dataset in "long_dataset" and filter it useing the 
  # 'min_dist' and "max_dist" definitions in the lifestage dictionary 
  if (distance) {
    mask_sliced <- long_dataset[["dist"]]$interp_like(var_array, method = "nearest")
    py_mask <- (mask_sliced$lat >= min_lat) & 
      (mask_sliced$lat <= max_lat) & 
      (mask_sliced >= filter_dict[[stage]][['min_dist']]) & 
      (mask_sliced <= filter_dict[[stage]][['max_dist']]) & 
      (!mask_sliced$isnull())
  } else { #if distance is 'FALSE' then we filter using min_iso and max_iso
    mask_sliced <- long_dataset[["deptho"]]$interp_like(var_array, method = "nearest")
    py_mask <- (mask_sliced$lat >= min_lat) & 
      (mask_sliced$lat <= max_lat) & 
      (mask_sliced >= filter_dict[[stage]][['min_iso']]) & 
      (mask_sliced <= filter_dict[[stage]][['max_iso']]) & 
      (!mask_sliced$isnull())
  }
  
  
  # 5. Apply the mask layer
  cat("Applying mask layer...\n")
  var_masked <- var_array$where(py_mask)
  
  # This is just a check: Identify remaining spatial dimensions to summarise means over
  all_dims <- var_masked$dims
  spatial_dims_to_drop <- all_dims[all_dims != "time"]
  
  #spatial_dims_to_drop are the spatial dimensions to take a mean over. It should
  # just be "lat" and "long"
  
  # 6. Taking the spatial mean 
  # this takes the mean over the defined latitude and the longitudinal extent
  # filtered by the mask.
  cat("Taking the spatial means...\n")
  spatial_mean_series <- var_masked$mean(dim = spatial_dims_to_drop, skipna = TRUE)
  spatial_mean_series$time
  # 7. Pull the monthly means
  # you now have monthly means across the defined lat long. We need to pull them into
  # R and align them with the correct time
  cat("Extracting values safe from memory overhead...\n")
  raw_values <- py_to_r(spatial_mean_series$values)
  
  # Manually construct R data frame with the extracted r_years
  temp_df <- data.frame(
    year = r_years,
    value = raw_values
  )
  
  # 8. Aggregate annual mean time series using R
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


# annual_mean_2d_fast <- function(var_dataset, long_dataset, ds_var_name, df_var_name, filter_dict,stage, precond, print_counts = TRUE) {
#   
#   # Start timer
#   start_time <- Sys.time()
#   cat(paste0("---", df_var_name, "---\n"))
#   
#   # Print initial info
#   print_ds_info(var_dataset, ds_var_name, print_counts)
#   
#   # Merge variable with mask file
#   ds_merged <- xr$merge(list(var_dataset, long_dataset))
#   print_ds_info(ds_merged, ds_var_name, print_counts)
#   
#   # Create list with months (R vectors auto-convert to Python lists)
#   month_list <- seq(filter_dict[[stage]][["min_month"]], filter_dict[[stage]][["max_month"]])
#   
#   # Select times and latitudes
#   # Using xr$slice for slicing dimensions via reticulate
#   ds_selected <- ds_merged$sel(
#     time = ds_merged$time$dt$month$isin(month_list), 
#     lat = py$slice(filter_dict[[stage]][["min_lat"]], filter_dict[[stage]][["max_lat"]])
#   )
#   print_ds_info(ds_selected, ds_var_name, print_counts)
#   
#   # Mask longs using bitwise AND (&) on the xarray DataArrays
#   ds_masked <- ds_selected$where(
#     (ds_selected$deptho >= filter_dict[[stage]][['min_iso']]) & (ds_selected$deptho <= filter_dict[[stage]][['max_iso']])
#   )
#   print_ds_info(ds_masked, ds_var_name, print_counts)
#   
#   # Calculate mean grouped by year
#   # Dim parameter accepts an R vector which automatically passes as a Python list
#   ds_mean <- ds_masked[[ds_var_name]]$groupby('time.year')$mean(dim = c('time', 'lon', 'lat'))
#   
#   # Create pandas dataframe and automatically convert it to an R data.frame
#   # Setting convert = TRUE ensures it returns a native R data frame
#   ds_df <- py_to_r(ds_mean$to_dataframe())



#   
#     # --- ADDED PRECOND LAG LOGIC ---
#   # If precond is TRUE, add a 1-year lag to the timeline data index
#   if (precond) {
#     cat("Applying 1-year environmental precondition lag (Years + 1)...\n")
#     # In a default pandas to_dataframe() conversion with groupby('time.year'),
#     # the years are stored as the index (rownames) of the R data frame.
#     current_years <- as.numeric(rownames(ds_df))
#     lagged_years <- current_years + 1
#     rownames(ds_df) <- as.character(lagged_years)
#   }
#   
#   # Rename column from dataset to dataframe variable name in R syntax
#   names(ds_df)[names(ds_df) == ds_var_name] <- df_var_name
#   
#   # End timer and print runtime
#   end_time <- Sys.time()
#   runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
#   cat(sprintf("Runtime: %.2f minutes\n", runtime_minutes))
#   
#   # Write output to CSV using standard R engine
#   #write.csv(ds_df, paste0("Output/hake_MOM6_", df_var_name, ".csv"), row.names = TRUE)
#   
#   return(ds_df)
# }
# 
# `