rm(list = ls())
library(htmlwidgets)  # Ensure this is loaded for onRender
library(sf)
library(tigris)  # Provides US boundary data
# Load the leaflet library
library(leaflet)
library(tidyverse)

# set working directory
wd <- getwd()
load(paste0(wd,"/workingRdata/Station_Data.RDATA"))
load(paste0(wd,"/workingRdata/TMAX_df.RDATA"))
load(paste0(wd,"/workingRdata/TMIN_df.RDATA"))

## Pull all files in stationdir saved to harddrive
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
us_states <- tigris::states(cb = TRUE)  # cb = TRUE for a generalization of the boundaries
load("WorkSpace.RDATA")
# Create a basic leaflet map
map <- leaflet() %>%
  setView(lng = -96, lat = 37.8, zoom = 4) %>%  # Centered on the US
  addTiles() %>%  # Add default OpenStreetMap tiles
  addPolygons(data = us_states, fillColor = "lightblue", weight = 1, opacity = 1, color = "black", 
              dashArray = "3", fillOpacity = 0.5)  # Adding U.S. states borders

# Add airport markers to the map
for (i in 1:nrow(airports)) {
  map <- map %>%
    addMarkers(lng = airports$longitude[i], 
               lat = airports$latitude[i], 
               popup = paste("ID:", airports$id[i], "<br>", 
                             "Name:", airports$name[i], "<br>", 
                             "State:", airports$state[i]))
}

# Add a click event to find the nearest airport
map <- map %>%
  onRender("
    function(el, x) {
      var airports = %s;  // This will be replaced with R data
      this.on('click', function(e) {
        var lat = e.latlng.lat;
        var lng = e.latlng.lng;

        // Function to calculate the distance
        function getDistance(lat1, lon1, lat2, lon2) {
          var R = 6371; // Radius of the earth in km
          var dLat = (lat2 - lat1) * Math.PI / 180; 
          var dLon = (lon2 - lon1) * Math.PI / 180; 
          var a = 
            Math.sin(dLat/2) * Math.sin(dLat/2) +
            Math.cos(lat1 * Math.PI / 180) * Math.cos(lat2 * Math.PI / 180) * 
            Math.sin(dLon/2) * Math.sin(dLon/2); 
          var c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a)); 
          var distance = R * c; // Distance in km
          return distance;
        }

        // Find the nearest airport
        var nearestAirport = null;
        var nearestDistance = Infinity;

        for (var i = 0; i < airports.length; i++) {
          var distance = getDistance(lat, lng, airports[i].lat, airports[i].lng);
          if (distance < nearestDistance) {
            nearestDistance = distance;
            nearestAirport = airports[i].id;
          }
        }

        alert('Nearest Airport ID: ' + nearestAirport);
      });
    }
  ", 
           airports = jsonlite::toJSON(
             airports[, c("id", "latitude", "longitude")],
             dataframe = "rows",
             rownames = FALSE
           ))

# Save the map to an HTML file
saveWidget(map, "airport_map.html", selfcontained = TRUE)
browseURL("airport_map.html")
