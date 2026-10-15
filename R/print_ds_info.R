#' Print number of records and dimensions of Zarr object 
#' once converted to match Python's double bracket dataset subsetting
#' using list() 
#' @author Megan L. Feddern

print_ds_info <- function(dataset, var_name, nan_bool = FALSE) {
  
  # Access the variable array using standard R list indexing
  var_array <- dataset[[var_name]]
  
  # In Python, .size and .shape are attributes
  denominator <- py_to_r(var_array$size)
  shape <- py_to_r(var_array$shape)
  
  # Format shape into a clean string (e.g., "(12, 360, 720)")
  shape_str <- paste0("(", paste(shape, collapse = ", "), ")")
  
  cat(sprintf("%d records with dimensions %s\n", denominator, shape_str))
  
  if (nan_bool) {
    # Call NumPy's isfinite on the python data object
    # .sum() converts the boolean array to a count
    numerator <- py_to_r(np$isfinite(var_array$data)$sum())
    percent <- (numerator / denominator) * 100
    
    cat(sprintf("%d records (%.2f%%) are not np.nan\n", numerator, percent))
  }
  
  return(invisible(NULL))
}