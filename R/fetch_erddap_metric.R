#' Accessing erdapp this is relevant to basin-scale indices (PDO, NPGO, ONI)
#' @author Megan L. Feddern


fetch_erddap_metric <- function(dataset_id, variable, server_url = "https://oceanview.pfeg.noaa.gov/erddap/") {
  
  # 1. Define fields to return
  requested_fields <- c("time", variable)
  
  # 2. Build the string constraint exactly as a character string (e.g., 'index="NPGO"')
  constraint_string <- paste0('index="', variable, '"')
  
  # 3. Use do.call to bypass rerddap's non-standard evaluation issues.
  # This forces R to evaluate the text string before making the API request.
  args_list <- list(
    dataset_id,
    url = server_url,
    fields = requested_fields,
    constraint_string
  )
  
  raw_data <- do.call(rerddap::tabledap, args_list)
  
  # 4. Return as a standard data frame
  return(as.data.frame(raw_data))
}
