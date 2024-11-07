###########
# Data Visualization & Figures
##########

library(ggplot2)
library(dplyr)
library(lubridate)
library(purrr)
library(tidyr)
library(RColorBrewer)

load(merged_TMAX_df,file = paste0(getwd(),"/workingRdata/TMAX_df.RDATA"))
load(merged_TMIN_df,file = paste0(getwd(),"/workingRdata/TMIN_df.RDATA"))
load(filter_station_data, file = paste0(getwd(),"/workingRdata/Station_Data.RDATA"))

##Reconvert to long format for graph making

TMAX_long <- merged_TMAX_df |>
  pivot_longer(cols = starts_with("Temp"),  # Select temperature columns
               names_to = "station_id", 
               values_to = "max_temp") |>
  mutate(station_id = sub("Temp ", "", station_id))

TMIN_long <- merged_TMIN_df |>
  pivot_longer(cols = starts_with("Temp"),
               names_to = "station_id", 
               values_to = "min_temp") |>
  mutate(station_id = sub("Temp ", "", station_id))

#A Few extra cleaning points, noticed while graphing

#Normalizing Tempurature values which previously were presented as "226 = 22.6 degrees"

TMAX_long <- TMAX_long |>
  mutate(max_temp = ifelse(!is.na(max_temp), max_temp / 10, NA))

TMIN_long <- TMIN_long |>
  mutate(min_temp = ifelse(!is.na(min_temp), min_temp / 10, NA))

#Removing extreme values

TMAX_long <- TMAX_long |>
  filter(abs(max_temp) <= 60 | is.na(max_temp))
#TMIN_long & TMAX_long are both showing a number of recorded tempuratures which go well beyond world records for low and high tempuratures, this normalizes to roughly the world record

TMIN_long <- TMIN_long |>
  filter(abs(min_temp) <= 60 | is.na(min_temp))

#Creating a monthly average for each station to handle missing values and reduce datapoints per graph

ann_avg_high_temp <- TMAX_long |>
  group_by(station_id, year) |>
  filter(!all(is.na(max_temp))) |>
  summarize(ann_avg_high = mean(as.numeric(max_temp), na.rm = TRUE), .groups = "drop")

ann_avg_low_temp <- TMIN_long |>
  group_by(station_id, year) |>
  filter(!all(is.na(min_temp))) |>
  summarize(ann_avg_low = mean(as.numeric(min_temp), na.rm = TRUE), .groups = "drop")

# Joining the annual average datasets

ann_avg_temps <- ann_avg_high_temp |>
  left_join(ann_avg_low_temp, by = c("station_id", "year"))

overall_avg_temp <- ann_avg_temps |>
  group_by(year) |>
  summarize(
    overall_avg_high = mean(ann_avg_high, na.rm = TRUE),
    overall_avg_low = mean(ann_avg_low, na.rm = TRUE)
  )

#Graph 1: A plot of monthly average high tempurature at our 100 stations from 1900 to present 

ggplot(ann_avg_temps, aes(x = year)) +
  geom_line(aes(y = ann_avg_high, group = station_id, color = "High Temp"), alpha = 0.2) +
  geom_line(aes(y = ann_avg_low, group = station_id, color = "Low Temp"), alpha = 0.2) +
  geom_line(data = overall_avg_temp, aes(y = overall_avg_high, color = "Overall High Temp"), size = 1.2) +
  geom_line(data = overall_avg_temp, aes(y = overall_avg_low, color = "Overall Low Temp"), size = 1.2) +
  geom_smooth(aes(y = ann_avg_high, group = 1, color = "High Temp Trend"), method = "lm", se = FALSE, size = 0.5) +
  geom_smooth(aes(y = ann_avg_low, group = 1, color = "Low Temp Trend"), method = "lm", se = FALSE, size = 0.5) +
  scale_color_manual(
    values = c("High Temp" = "red", "Low Temp" = "blue", "Overall High Temp" = "darkred", "Overall Low Temp" = "darkblue", 
               "High Temp Trend" = "white", "Low Temp Trend" = "white"),
    labels = c("High Temp", "Low Temp", "Overall High Temp", "Overall Low Temp", "High Temp Trend", "Low Temp Trend")
  ) +
  labs(
    title = "Annual Average Daily High and Low Temperatures\n at Sampled Stations (1900-2024)",
    x = "Year",
    y = "Average Temperature (°C)",
    color = "Temperature Type"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        plot.title.position = "plot",
        legend.position = "right"
  ) +
  scale_y_continuous(
    breaks = seq(
      floor(-15),
      ceiling(30),
      by = 5
    ),
    sec.axis = sec_axis(~. * 1, breaks = seq(floor(-15), ceiling(30), by = 5))
  )

# More data prep
# Creating dataframes to represent highest annual temp and lowest annual temp

ann_extr_max <- TMAX_long |>
  group_by(station_id, year) |>
  summarize(
    max_temp_reached = max(max_temp, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(is.finite(max_temp_reached))

ann_extr_min <- TMIN_long |>
  group_by(station_id, year) |>
  summarize(
    min_temp_reached = min(min_temp, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(is.finite(min_temp_reached))

ann_extr_temps <- ann_extr_min |>
  left_join(ann_extr_max, by = c("station_id", "year"))|>
  mutate(
    ann_extr_max = as.numeric(max_temp_reached),
    ann_extr_min = as.numeric(min_temp_reached)) |>
  arrange(station_id, year) |>
  group_by(station_id) |>
  mutate(
    year_gap = year - lag(year, default = first(year)),
    segment_id = cumsum(year_gap > 5)
  ) |>
  ungroup()

overall_extr_temp <- ann_extr_temps |>
  group_by(year) |>
  summarize(
    overall_extr_max = max(ann_extr_max, na.rm = TRUE),
    overall_extr_min = min(ann_extr_min, na.rm = TRUE),
    avg_extr_max = mean(ann_extr_max, na.rm = TRUE),
    avg_extr_min = mean(ann_extr_min, na.rm = TRUE)
  )

#Graph 2: A plot of the highest and lowest annual tempurature at sampled stations

ggplot(ann_extr_temps, aes(x = year)) +
  geom_line(aes(y = ann_extr_max, group = interaction(station_id, segment_id), color = "Max Temp"), alpha = 0.2) +
  geom_line(aes(y = ann_extr_min, group = interaction(station_id, segment_id), color = "Min Temp"), alpha = 0.2) +
  geom_line(data = overall_extr_temp, aes(y = overall_extr_max, color = "Overall Max Temp"), size = 0.6) +
  geom_line(data = overall_extr_temp, aes(y = overall_extr_min, color = "Overall Min Temp"), size = 0.6) +
  geom_line(data = overall_extr_temp, aes(y = avg_extr_max, color = "Average Max Temp"), size = 0.6) +
  geom_line(data = overall_extr_temp, aes(y = avg_extr_min, color = "Average Min Temp"), size = 0.6) +
  geom_smooth(aes(y = ann_extr_max, group = 1, color = "Max Temp Trend"), method = "lm", se = FALSE, size = 0.5) +
  geom_smooth(aes(y = ann_extr_min, group = 1, color = "Min Temp Trend"), method = "lm", se = FALSE, size = 0.5) +
  scale_color_manual(
    values = c("Max Temp" = "red", "Min Temp" = "blue", "Overall Max Temp" = "darkred", "Overall Min Temp" = "darkblue", 
               "Average Max Temp" = "black", "Average Min Temp" = "black", "Max Temp Trend" = "white", "Min Temp Trend" = "white"),
    labels = c("High Temp", "Low Temp", "Overall High Temp", "Overall Low Temp", "Average Max Temp", "Average Min Temp", 
               "Max Temp Trend", "Min Temp Trend")
  ) +
  labs(
    title = "Annual Maximum and Minimum Temperatures\n at Sampled Stations (1900-2024)",
    x = "Year",
    y = "Average Temperature (°C)",
    color = "Temperature Type"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        plot.title.position = "plot",
        legend.position = "right"
  ) +
  scale_y_continuous(
    breaks = seq(
      floor(-40),
      ceiling(40),
      by = 5
    ),
    sec.axis = sec_axis(~. * 1, breaks = seq(floor(-40), ceiling(40), by = 5))
  )

# More Data Work
# Creating an Average only for the month of August (to attempt and measure only a hot month, see below for same but for January)

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

#Graph 3: Average August Tempuratures at Sampled Stations (1900-2024)

ggplot(AUG_avg_temps, aes(x = year)) +
  geom_line(aes(y = AUG_avg_high, group = interaction(station_id, segment_id), color = "High Temp"), alpha = 0.2) +
  geom_line(aes(y = AUG_avg_low, group = interaction(station_id, segment_id), color = "Low Temp"), alpha = 0.2) +
  geom_line(data = overall_AUG_avg_temp, aes(y = overall_avg_high, color = "Average High Temp"), size = 0.8) +
  geom_line(data = overall_AUG_avg_temp, aes(y = overall_avg_low, color = "Average Low Temp"), size = 0.8) +
  geom_smooth(aes(y = AUG_avg_high, group = 1, color = "High Temp Trend"), method = "lm", se = FALSE, size = 0.5) +
  geom_smooth(aes(y = AUG_avg_low, group = 1, color = "Low Temp Trend"), method = "lm", se = FALSE, size = 0.5) +
   scale_color_manual(
    values = c("High Temp" = "red", "Low Temp" = "blue", "Average High Temp" = "darkred", "Average Low Temp" = "darkblue",
               "High Temp Trend" = "white", "Low Temp Trend" = "white"),
    labels = c("High Temp", "Low Temp", "Overall High Temp", "Overall Low Temp", "High Temp Trend", "Low Temp Trend")
  ) +
  labs(
    title = "Average High and Low Temperatures at Sampled\n Stations During the Month of August (1900-2024)",
    x = "Year",
    y = "Average Temperature (°C)",
    color = "Temperature Type"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        plot.title.position = "plot",
        legend.position = "right"
  )+
  scale_y_continuous(
    breaks = seq(
      floor(0),
      ceiling(40),
      by = 5 
    ),
    sec.axis = sec_axis(~. * 1, breaks = seq(floor(0), ceiling(40), by = 5))
  )

# More Data Work
# Creating an Average only for the month of February (to attempt and measure only a cold month, see below for same but for January)

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

#Graph 4: Average February Tempuratures at Sampled Stations (1900-2024)

ggplot(FEB_avg_temps, aes(x = year)) +
  geom_line(aes(y = FEB_avg_high, group = interaction(station_id, segment_id), color = "High Temp"), alpha = 0.2) +
  geom_line(aes(y = FEB_avg_low, group = interaction(station_id, segment_id), color = "Low Temp"), alpha = 0.2) +
  geom_line(data = overall_FEB_avg_temp, aes(y = overall_avg_high, color = "Average High Temp"), size = 1) +
  geom_line(data = overall_FEB_avg_temp, aes(y = overall_avg_low, color = "Average Low Temp"), size = 1) +
  geom_smooth(aes(y = FEB_avg_high, group = 1, color = "High Temp Trend"), method = "lm", se = FALSE, size = 0.5) +
  geom_smooth(aes(y = FEB_avg_low, group = 1, color = "Low Temp Trend"), method = "lm", se = FALSE, size = 0.5) +
  scale_color_manual(
    values = c("High Temp" = "red", "Low Temp" = "blue", "Average High Temp" = "darkred", "Average Low Temp" = "darkblue",
               "High Temp Trend" = "white", "Low Temp Trend" = "white"),
    labels = c("High Temp", "Low Temp", "Overall High Temp", "Overall Low Temp", "High Temp Trend", "Low Temp Trend")
  ) +
  labs(
    title = "Average High and Low Temperatures at Sampled\n Stations During the Month of February (1900-2024)",
    x = "Year",
    y = "Average Temperature (°C)",
    color = "Temperature Type"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        plot.title.position = "plot",
        legend.position = "right"
  ) +
  scale_y_continuous(
    breaks = seq(
      floor(-20),
      ceiling(20),
      by = 5
    ),
    sec.axis = sec_axis(~. * 1, breaks = seq(floor(-20), ceiling(20), by = 5))
)

# More Data Work
# Our final two graphs analyze change in temperature by elevation and by latitude. 
# For this we have to process and join our data with the filter_station_data

stations_clean = filter_station_data |>
  select(station_id = id, latitude, elevation)

temps_n_stations = stations_clean |>
  left_join(ann_avg_temps, by = c("station_id"))

temps_n_stations <- temps_n_stations |>
  mutate(elevation_group = cut(
    elevation,
    breaks = quantile(elevation, probs = seq(0, 1, by = 0.25)),
    include.lowest = TRUE,
    labels = FALSE 
  )) |>
  mutate(elevation_group = factor(
    elevation_group,
    levels = 1:4,
    labels = sapply(quantile(elevation, probs = seq(0, 1, by = 0.25))[-5], 
                    function(x) paste0(
                      round(x), " ", "–", " ", round(quantile(elevation, probs = seq(0, 1, by = 0.25))[
                        which(quantile(elevation, probs = seq(0, 1, by = 0.25)) == x) + 1])))))

# Table 5: Analyzing change in temperature by elevation  

ggplot(temps_n_stations, aes(x = year)) +
  geom_line(aes(y = ann_avg_high, color = elevation), size = 0.5) +
  geom_smooth(aes(y = ann_avg_high, linetype = elevation_group, group = elevation_group), 
              method = "lm", se = FALSE, size = 1, color = "red") +  
  scale_color_gradient(low = "yellow", high = "darkgreen", name = "Elevation") +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotted", "dotdash", "longdash"),
    name = "Elevation Group"
  ) +
  labs(
    title = "Change in Average High Tempurature\n by Elevation Group (1900-2024)",
    x = "Year",
    y = "Temperature (°C)",
    color = "Elevation"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right",
        plot.title.position = "plot") +
  scale_y_continuous(
    breaks = seq(
      floor(-10),
      ceiling(50),
      by = 5
    ),
    sec.axis = sec_axis(~. * 1, breaks = seq(floor(-10), ceiling(50), by = 5))
  )


# More Data Prep

temps_n_stations <- temps_n_stations |>
  mutate(latitude_group = cut(
    latitude,
    breaks = quantile(latitude, probs = seq(0, 1, by = 0.25)),
    include.lowest = TRUE,
    labels = FALSE 
  )) |>
  mutate(latitude_group = factor(
    latitude_group,
    levels = 1:4,
    labels = sapply(quantile(latitude, probs = seq(0, 1, by = 0.25))[-5], 
                    function(x) paste0(
                      round(x), " ", "–", " ", round(quantile(latitude, probs = seq(0, 1, by = 0.25))[
                        which(quantile(latitude, probs = seq(0, 1, by = 0.25)) == x) + 1])))))

# Table 6: Analyzing change in temperature by latitude  

ggplot(temps_n_stations, aes(x = year)) +
  geom_line(aes(y = ann_avg_high, color = latitude), size = 0.5) +
  geom_smooth(aes(y = ann_avg_high, linetype = latitude_group, group = latitude_group), 
              method = "lm", se = FALSE, size = 1, color = "red") +  
  scale_color_gradient(low = "yellow", high = "darkgreen", name = "Latitude") +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotted", "dotdash", "longdash"),
    name = "Latitude Group"
  ) +
  labs(
    title = "Change in Average High Tempurature\n by Latitude Group (1900-2024)",
    x = "Year",
    y = "Temperature (°C)",
    color = "Latitude"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right",
        plot.title.position = "plot") +
  scale_y_continuous(
    breaks = seq(
      floor(-10),
      ceiling(50),
      by = 5
    ),
    sec.axis = sec_axis(~. * 1, breaks = seq(floor(-10), ceiling(50), by = 5))
  )



