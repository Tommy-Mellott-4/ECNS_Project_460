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

# Function to generate a temperature plot for an airport
generate_temp_plot <- function(airport_id) {
  # Filter temperature data for the given airport ID
  airport_temp__ = merged_TMAX_df[c(paste0("Temp ",airport_id),"jday")]
  colnames(airport_temp__) = c("temperature","date")
  airport_temp = na.omit(airport_temp__)
  # Create a plot
  p <- ggplot(airport_temp, aes(x = date, y = temperature)) +
    geom_line(color = "blue") +
    labs(title = paste("Temperature for Airport ID:", airport_id),
         x = "Date", y = "Temperature (°C)") +
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
  temp_plot <- generate_temp_plot(airport_id)
  print(i)
  map <- map %>%
    addMarkers(lng = airports$longitude[i], 
               lat = airports$latitude[i], 
               popup = paste("ID:", airports$id[i], "<br>", 
                             "Name:", airports$name[i], "<br>", 
                             "State:", airports$state[i], "<br>",
                            "<img src='", temp_plot, "' width = '300'>"))
}

map




# Save the map to an HTML file
saveWidget(map, "airport_map.html", selfcontained = TRUE)
browseURL("airport_map.html")
