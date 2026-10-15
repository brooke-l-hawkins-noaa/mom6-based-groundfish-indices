#' Developing a function for 2-dimensional data. This is relevant for MLD, SSH, and SST.
#' using list() 
#' @author Megan L. Feddern

annual_mean_2d <- function(var_dataset, long_dataset, distance, ds_var_name, df_var_name, filter_dict, stage, precond, print_counts = TRUE) {
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
  
  # 2. Slice spatial coordinate boundaries on the main 1D variable array FIRST
  # This isolates only your 20-latitude window right out of the box
  var_selected <- var_selected$sel(lat = py$slice(min_lat, max_lat))
  var_array <- var_selected[[ds_var_name]]
  
  # 3. Dynamic Coordinate Name Standardization for long_dataset
  cat("Standardizing mask coordinate names...\n")
  dim_keys <- py_to_r(py$list(long_dataset$dims))
  
  rename_list <- list()
  if ("latitude" %in% dim_keys)  rename_list$latitude <- "lat"
  if ("longitude" %in% dim_keys) rename_list$longitude <- "lon"
  
  if (length(rename_list) > 0) {
    long_dataset <- long_dataset$rename(rename_list)
  }
  
  # 4. Use Xarray's interp_like to ensure alignment
  cat("Interpolating mask grid onto variable layout...\n")
  if (distance) {
    mask_sliced <- long_dataset[["dist"]]$interp_like(var_array, method = "nearest")
    py_mask <- (mask_sliced$lat >= min_lat) & 
      (mask_sliced$lat <= max_lat) & 
      (mask_sliced >= filter_dict[[stage]][['min_dist']]) & 
      (mask_sliced <= filter_dict[[stage]][['max_dist']]) & 
      (!mask_sliced$isnull())
  } else {
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
  
  # Identify remaining spatial dimensions to collapse
  all_dims <- var_masked$dims
  spatial_dims_to_drop <- all_dims[all_dims != "time"]
  
  # 6. Collapse spatial dimensions natively inside Python
  cat("Collapsing spatial dimensions natively...\n")
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
  
  # Explicitly trigger system garbage collection
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