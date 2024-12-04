rm(list = ls())
library(htmlwidgets)  # Ensure this is loaded for onRender
library(sf)
library(tigris)  # Provides US boundary data
# Load the leaflet library
library(leaflet)
library(tidyverse)
library(ggplot2)
load("WorkSpace.RDATA")

# set working directory
wd <- getwd()
load(paste0(wd,"/workingRdata/Station_Data.RDATA"))
load(paste0(wd,"/workingRdata/TMAX_df.RDATA"))
load(paste0(wd,"/workingRdata/TMIN_df.RDATA"))

## Pull all files in stationdir saved to hard drive
Station_Info_Files <- list.files(path = stationsdir,full.names = T)
### Separate the RData station info for the Temp and Precip Data
US_PRCP_Station_Info <- Station_Info_Files[4:8]
US_Temp_Station_Info <- Station_Info_Files[8:11]
### show all the paths to each saved RData info on Temp Station
US_Temp_Station_Info
### Load only the ones that has data since 1900
load(US_Temp_Station_Info[1])

### making an object to hold either all the stations since 1900 or all the ones I have previously saved
airports <- filter_station_data
# Fetch state boundaries
tigris::tigris_cache_dir(NULL)
us_states <- tigris::states(cb = TRUE)  # cb = TRUE for a generalization of the boundaries

#Tables + Code
TMAX_long <- merged_TMAX_df |>
  pivot_longer(cols = starts_with("Temp"),  # Select temperature columns
               names_to = "station_id", 
               values_to = "max_temp") |>
  mutate(station_id = sub("Temp ", "", station_id)) |>
  mutate(max_temp = ifelse(!is.na(max_temp), max_temp / 10, NA)) |>
  filter(abs(max_temp) <= 60 | is.na(max_temp))

TMIN_long <- merged_TMIN_df |>
  pivot_longer(cols = starts_with("Temp"),
               names_to = "station_id", 
               values_to = "min_temp") |>
  mutate(station_id = sub("Temp ", "", station_id)) |>
  mutate(min_temp = ifelse(!is.na(min_temp), min_temp / 10, NA)) |>
  filter(abs(min_temp) <= 60 | is.na(min_temp))
# Removes extreme values
# TMIN_long & TMAX_long are both showing a number of recorded tempuratures which go well beyond world records for low and high tempuratures, this normalizes to roughly the world record

TMAX_long_FEB = TMAX_long |>
  filter(month == 2)

TMIN_long_FEB = TMIN_long |>
  filter(month == 2)

FEB_avg_high_temp <- TMAX_long_FEB |>
  group_by(station_id, year) |>
  filter(!all(is.na(max_temp))) |>
  summarize(FEB_avg_high = mean(as.numeric(max_temp), na.rm = TRUE), .groups = "drop")

FEB_avg_low_temp <- TMIN_long_FEB |>
  group_by(station_id, year) |>
  filter(!all(is.na(min_temp))) |>
  summarize(FEB_avg_low = mean(as.numeric(min_temp), na.rm = TRUE), .groups = "drop")

FEB_avg_temps <- FEB_avg_low_temp |>
  left_join(FEB_avg_high_temp, by = c("station_id", "year")) |>
  arrange(station_id, year) |>
  group_by(station_id) |>
  mutate(
    year_gap = year - lag(year, default = first(year)),
    segment_id = cumsum(year_gap > 5)
  ) |>
  ungroup()

overall_FEB_avg_temp <- FEB_avg_temps |>
  group_by(year) |>
  summarize(
    overall_avg_high = mean(FEB_avg_high, na.rm = TRUE),
    overall_avg_low = mean(FEB_avg_low, na.rm = TRUE)
  )

# August Code
TMAX_long_AUG = TMAX_long |>
  filter(month == 8)

TMIN_long_AUG = TMIN_long |>
  filter(month == 8)

AUG_avg_high_temp <- TMAX_long_AUG |>
  group_by(station_id, year) |>
  filter(!all(is.na(max_temp))) |>
  summarize(AUG_avg_high = mean(as.numeric(max_temp), na.rm = TRUE), .groups = "drop")

AUG_avg_low_temp <- TMIN_long_AUG |>
  group_by(station_id, year) |>
  filter(!all(is.na(min_temp))) |>
  summarize(AUG_avg_low = mean(as.numeric(min_temp), na.rm = TRUE), .groups = "drop")

AUG_avg_temps <- AUG_avg_low_temp |>
  left_join(AUG_avg_high_temp, by = c("station_id", "year")) |>
  arrange(station_id, year) |>
  group_by(station_id) |>
  mutate(
    year_gap = year - lag(year, default = first(year)),
    segment_id = cumsum(year_gap > 5)
  ) |>
  ungroup()

overall_AUG_avg_temp <- AUG_avg_temps |>
  group_by(year) |>
  summarize(
    overall_avg_high = mean(AUG_avg_high, na.rm = TRUE),
    overall_avg_low = mean(AUG_avg_low, na.rm = TRUE)
  )

# Function to generate a temperature plot for an airport
generate_February_plot <- function(airport_id) {
  # Create a plot
  p <- ggplot(FEB_avg_temps, aes(x = year)) +
    # Other station lines (faint)
    geom_line(aes(
      y = FEB_avg_high,
      group = station_id
    ), alpha = 0.1, color = "red") +
    geom_line(aes(
      y = FEB_avg_low,
      group = station_id
    ), alpha = 0.1, color = "blue") +
    # Highlighted lines for the current station
    geom_line(data = FEB_avg_temps %>% filter(station_id == airport_id),
              aes(y = FEB_avg_high, group = station_id,
                  color = "Station High Temp"), size = 1) +
    geom_line(data = FEB_avg_temps %>% filter(station_id == airport_id),
              aes(y = FEB_avg_low, group = station_id,
                  color = "Station Low Temp"), size = 1) +
    # (same as in the static graph above)
    geom_line(data = overall_FEB_avg_temp, 
              aes(y = overall_avg_high, color = "Average High Temp"), 
              size = 1) +
    geom_line(data = overall_FEB_avg_temp, 
              aes(y = overall_avg_low, color = "Average Low Temp"), 
              size = 1) +
    scale_color_manual(
      values = c(
        "Average High Temp" = "darkred",
        "Average Low Temp" = "darkblue",
        "Station High Temp" = "black", 
        "Station Low Temp" = "darkgreen"
      )
    ) +
    labs(
      title = paste("Average February High and Low Temperatures at Station", airport_id),
      x = "Year",
      y = "Average Temperature (°C)",
      color = "Temperature Type"
    ) +
    theme_minimal()
  
  
  # Save the plot as a temporary image file
  temp_file <- tempfile(fileext = ".png")
  ggsave(temp_file, plot = p, width = 6, height = 4)
  
  # Convert the image to base64 for embedding in the map
  base64_image <- base64enc::dataURI(file = temp_file, mime = "image/png")
  return(base64_image)
}


# Function to generate a temperature plot for an airport
generate_August_plot <- function(airport_id) {
  # Create a plot
  p <- ggplot(AUG_avg_temps, aes(x = year)) +
    # Other station lines (faint)
    geom_line(aes(
      y = AUG_avg_high,
      group = station_id
    ), alpha = 0.1, color = "red") +
    geom_line(aes(
      y = AUG_avg_low,
      group = station_id
    ), alpha = 0.1, color = "blue") +
    # Highlighted lines for the current station
    geom_line(data = AUG_avg_temps %>% filter(station_id == airport_id),
              aes(y = AUG_avg_high, group = station_id,
                  color = "Station High Temp"), size = 1) +
    geom_line(data = AUG_avg_temps %>% filter(station_id == airport_id),
              aes(y = AUG_avg_low, group = station_id,
                  color = "Station Low Temp"), size = 1) +
    # (same as in the static graph above)
    geom_line(data = overall_AUG_avg_temp, 
              aes(y = overall_avg_high, color = "Average High Temp"), 
              size = 1) +
    geom_line(data = overall_AUG_avg_temp, 
              aes(y = overall_avg_low, color = "Average Low Temp"), 
              size = 1) +
    scale_color_manual(
      values = c(
        "Average High Temp" = "darkred",
        "Average Low Temp" = "darkblue",
        "Station High Temp" = "black", 
        "Station Low Temp" = "darkgreen"
      )
    ) +
    labs(
      title = paste("Average August High and Low Temperatures at Station", airport_id),
      x = "Year",
      y = "Average Temperature (°C)",
      color = "Temperature Type"
    ) +
    theme_minimal()
  
  # Save the plot as a temporary image file
  temp_file <- tempfile(fileext = ".png")
  ggsave(temp_file, plot = p, width = 6, height = 4)
  
  # Convert the image to base64 for embedding in the map
  base64_image <- base64enc::dataURI(file = temp_file, mime = "image/png")
  return(base64_image)
}



# Create a basic leaflet map
map <- leaflet() %>%
  setView(lng = -96, lat = 37.8, zoom = 4) %>%  # Centered on the US
  addTiles() %>%  # Add default OpenStreetMap tiles
  addPolygons(data = us_states, fillColor = "lightblue", weight = 1, opacity = 1, color = "black", 
              dashArray = "3", fillOpacity = 0.5)  # Adding U.S. states borders

# Add airport markers to the map
for (i in 1:nrow(airports)) {
  
  airport_id <- airports$id[i]
  February_plot <- generate_February_plot(airport_id)
  August_plot <- generate_August_plot(airport_id)
  popup_opts <- popupOptions(maxWidth = 500) # adjust maxWidth of popup
  print(i)
  map <- map %>%
    addMarkers(lng = airports$longitude[i], 
               lat = airports$latitude[i], 
               popup = paste("ID:", airports$id[i], "<br>", 
                             "Name:", airports$name[i], "<br>", 
                             "State:", airports$state[i], "<br>",
                            "<img src='", February_plot, "' width = '500'>",
                            "<img src='", August_plot, "' width = '500'>"),
               popupOptions = popup_opts)
}

map




# Save the map to an HTML file
saveWidget(map, "airport_map.html", selfcontained = TRUE)
browseURL("airport_map.html")
