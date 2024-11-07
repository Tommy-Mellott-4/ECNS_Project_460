###############################################################################
# Load the required library
library(digest)
library(tidyverse)
library(lubridate)
rm(list = ls())
load('WorkSpace.RData'); cleanup()
###############################################################################
my_require(c('sp'))
my_require('readr')
## See all airport's that we have saved data for
temp_files_downloaded <- list.files(RDatadir, full.names = T)
### Now form the remaining airport ID in a nice way to filter later
temp_files_noFN <-list.files(RDatadir)
values_to_keep <-gsub("\\.RData$","",temp_files_noFN)
## Pull all files in stationdir saved to harddrive
Station_Info_Files <- list.files(path = stationsdir,full.names = T)
### Separate the RData station info for the Temp and Precip Data
US_Temp_Station_Info <- Station_Info_Files[8:11]

### show all the paths to each saved RData info on Temp Station
load(US_Temp_Station_Info[1])
## Filter the station data to keep only the 100 randomly selected from earlier
filter_station_data <- US_stations |>
  filter(id %in% values_to_keep)

which.min(filter_station_data$first_year)

load(temp_files_downloaded[1])

wdata_filtered <- wdata|>
  select(1:4,contains("VALUE"))|>
  filter(element %in% c("TMAX","TMIN"))

id <- unique(wdata_filtered$id)
long_data <- wdata_filtered |>
  pivot_longer(
    cols = starts_with("VALUE"),
    names_to = "Day",
    values_to = paste0("Temp ",id)
  )|>
  filter(year>=1900)

df_TMAX <- long_data |>
  filter(element == "TMAX")|>
  mutate(Day = gsub("VALUE","",Day))|>
  mutate(jday = julian(make_date(year,month,Day)))|>
  select(year,month,Day,jday,everything())|>
  select(-id,-element)|>



df_TMIN <- long_data |>
  filter(element != "TMAX")|>
  mutate(Day = gsub("VALUE","",Day))|>
  mutate(jday = julian(make_date(year,month,Day)))|>
  select(year,month,Day,jday,-element,-id,everything())|>
  select(-id,-element)

# Define the start and end dates
start_date <- as.Date("1900-01-01")
end_date <- Sys.Date()

# Generate a sequence of all dates from start_date to end_date
dates <- seq.Date(from = start_date, to = end_date, by = "day")

# Create a dataframe with the dates
dates_df <- data.frame(date = dates)
tail(dates_df)
tail(merged_TMAX_df$jday)
# Extract year, month, and day into separate columns
dates_df <- dates_df %>%
  mutate(
    year = year(date),
    month = month(date),
    day = day(date)
  )|>
  mutate(jday = julian(make_date(year,month,day)))
  

# Display the first few rows
head(dates_df)
jday<-dates_df$jday
merged_TMAX_df <- dates_df
merged_TMIN_df <- dates_df
for(i in 1:100){
    load(temp_files_downloaded[i])
    
    wdata_filtered <- wdata|>
      select(1:4,contains("VALUE"))|>
      filter(element %in% c("TMAX","TMIN"))
    
    id <- unique(wdata_filtered$id)
    long_data <- wdata_filtered |>
      pivot_longer(
        cols = starts_with("VALUE"),
        names_to = "Day",
        values_to = paste0("Temp ",id)
      )|>
      filter(year>=1900)
    
    MAX <- long_data |>
      filter(element == "TMAX")|>
      mutate(Day = gsub("VALUE","",Day))|>
      mutate(jday = julian(make_date(year,month,Day)))|>
      select(year,month,Day,jday,everything())|>
      select(jday,paste0("Temp ",id))|>
      filter(!is.na(jday))
    
    
    
    MIN <- long_data |>
      filter(element != "TMAX")|>
      mutate(Day = gsub("VALUE","",Day))|>
      mutate(jday = julian(make_date(year,month,Day)))|>
      select(year,month,Day,jday,-element,-id,everything())|>
      select(jday,paste0("Temp ",id))|>
      filter(!is.na(jday))
    
    
    merged_TMAX_df <- merged_TMAX_df |>
      full_join(MAX,by = "jday")
      
    merged_TMIN_df <- merged_TMIN_df |>
      full_join(MIN,by = "jday")
    
  print(i)
}



save(merged_TMAX_df,file = paste0(getwd(),"/workingRdata/TMAX_df.RDATA"))
save(merged_TMIN_df,file = paste0(getwd(),"/workingRdata/TMIN_df.RDATA"))
save(filter_station_data, file = paste0(getwd(),"/workingRdata/Station_Data.RDATA"))

