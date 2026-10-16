#install.packages(c("rerddap","httr"))
library(dplyr)
library(tidyr)
library(corrplot)
library(rerddap)
library(ggplot2)
library(rerddap)
library(ggplot2)
library(lubridate)
library(rerddap)
library(purrr)
library(stringr)
server_url <- "https://oceanview.pfeg.noaa.gov/erddap/"

# Creating a list for lifesatge (NOTE make this a R object so it can be accessed from multiple scripts once finalized)
# assign months, latitudes, and depths to filter
# all are integers, depths are in meters
# filters are inclusive (closed interval), e.g. depth 0-180 includes depths 0 and 180 and every depth between

lifestage_dict <- list(
  JA = list( # spawn
    min_month = 1,
    max_month = 4
  ), #long extent shelf break
  AS = list( # egg
    min_month = 4,
    max_month = 9
  )
)


# Define the clean rerddap function
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
  
  raw_data <- do.call(tabledap, args_list)
  
  # 4. Return as a standard data frame
  return(as.data.frame(raw_data))
}

climate_datasets <- data.frame(
  variable=c("NPGO", "PDO"),
  dataset_id=c( "cciea_OC_NPGO","cciea_OC_PDO") 
)

# 2. Fetch each dataset individually using the function
compiled_dataset<-NULL
for(i in 1:length(climate_datasets)){
    temp<-fetch_erddap_metric(dataset_id = climate_datasets$dataset_id[i], variable = climate_datasets$variable[i])
    # Fetch the raw data for the current row
    var_name <- climate_datasets$variable[i]
    # 3. Combine the data dynamically
    if (is.null(compiled_dataset)) {
      # First iteration: initialize the dataset with the first metric
      compiled_dataset <- temp
    } else {
      # Subsequent iterations: join side-by-side on the 'time' column
      compiled_dataset <- full_join(compiled_dataset, temp, by = "time")
    }
    # Clean columns inside the loop so joining works seamlessly
    
}
     
compiled_dataset$month <- month(compiled_dataset$time)
compiled_dataset$day   <- format(compiled_dataset$time, "%d")
compiled_dataset$year   <- year(compiled_dataset$time)

cleaned_dataset<-compiled_dataset

# 3. Summarize metrics across each element of lifestage_dict dynamically
lifestage_summaries <- lapply(names(lifestage_dict), function(stage_name) {
  
  mn <- lifestage_dict[[stage_name]]$min_month
  mx <- lifestage_dict[[stage_name]]$max_month
  

    stage_data <- cleaned_dataset %>% filter(month >= mn & month <= mx)
 
  
  summary_stats <- stage_data %>%
    group_by(year) %>%
    summarise(
      across(
        all_of(climate_datasets$variable), 
        ~ mean(.x, na.rm = TRUE), 
        .names = "{.col}"  # Creates: mean_NPGO, mean_PDO, etc.
      ),
      .groups = 'drop'
    ) %>%
    mutate(lifestage = stage_name)
  
  return(summary_stats)
}) %>% 
  bind_rows()

# 4. Reshape the dataset into a wide format grouped by Year
wide_lifestage_summary <- lifestage_summaries %>%
  pivot_wider(
    names_from = lifestage,
    values_from = c(-year, -lifestage),
    # This reshapes names from "mean_NPGO" + "JA" into "mean_NPGO_JA"
    names_glue = "{.value}{lifestage}"
  )

# View your tidy, single-row-per-year database structure
print(head(wide_lifestage_summary))

# 5. Calculate lagged variables
lag<-1
preconditioning <-wide_lifestage_summary%>%
  mutate(year=year+lag)%>%
  select(c(year,contains("AS")))%>%
  rename_with(~ str_replace(.x, "_AS$", "pre"), ends_with("_AS"))

# 6. Full dataset
full_dataset <- left_join(preconditioning,wide_lifestage_summary)
write.csv(full_dataset, "Output/ERRDAP_var.csv", row.names = FALSE)
