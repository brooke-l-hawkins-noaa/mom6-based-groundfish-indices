#' Calculate coast angle
#' @author Megan L. Feddern

get_coast_angle <- function(index, points_df, window_size = 15) {
  half_w <- floor(window_size / 2)
  start_idx <- max(1, index - half_w)
  end_idx <- min(nrow(points_df), index + half_w)
  
  # Subset local neighborhood coordinates
  subset_pts <- points_df[start_idx:end_idx, ]
  
  # Linear regression: lat ~ lon 
  # Note: For physical accuracy, projecting your data frame points to a 
  # metric CRS (like local UTM meters) before running lm() is highly recommended!
  fit <- lm(lat ~ lon, data = subset_pts)
  slope <- coef(fit)[2]
  
  # Calculate orientation trajectory angle counter-clockwise from East
  angle_rad <- atan(slope) 
  return(angle_rad)
}