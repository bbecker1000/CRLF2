library(dplyr)
library(here)
library(tidyverse)
library(ggplot2)
library(lubridate)

setwd(here::here("code"))

# importing water temperature data
band_temps <- read_csv(here::here("data", "watertemp_banducci.csv")) %>% 
  mutate(timestamp = mdy_hm(timestamp),
         temperature = (value - 32) * (5/9)) %>% 
  select(timestamp, temperature)
rodeo_temps <- read_csv(here::here("data", "watertemp_rodeo.csv")) %>% 
  mutate(timestamp = mdy_hm(timestamp),
         temperature = if_else(value > 30, ((value - 32) * (5/9)), (value))) %>% 
  select(timestamp, temperature)

# importing survey data
band_surveys <- read_csv(here::here("data", "onset_of_breeding.csv")) %>% 
  filter(LocationID == "RC07") %>% 
  mutate(beginning_WY = paste0(BRDYEAR-1, "1001"),
         beginning_WY = ymd(beginning_WY),
         breeding_date = beginning_WY + ddays(first_breeding)) %>% 
  select(-MaxD, -MaxD_proportion, -MaxD_yearly, -WaterTemp, -AirTemp, -Watershed)
rodeo_surveys <- read_csv(here::here("data", "onset_of_breeding.csv")) %>% 
  filter(LocationID == "RL02") %>% 
  mutate(beginning_WY = paste0(BRDYEAR-1, "1001"),
         beginning_WY = ymd(beginning_WY),
         breeding_date = beginning_WY + ddays(first_breeding)) %>% 
  select(-MaxD, -MaxD_proportion, -MaxD_yearly, -WaterTemp, -AirTemp, -Watershed)


# getting daily average temps, adding degree days
# to calculate degree days, we want to add together all the days where
# the temperature was above a certain threshold temperature

temp_threshold = 10 # can change depending on threshold we choose

band_daily_temps <- band_temps %>% 
  mutate(date = date(timestamp)) %>% 
  group_by(date) %>% 
  summarize(mean_temp = mean(temperature)) %>% 
  ungroup() %>% 
  mutate(degree_days = if_else(mean_temp > temp_threshold, (mean_temp - temp_threshold), 0),
         BRDYEAR = if_else(month(date) > 9, year(date) + 1, year(date)),
         beginning_WY = ymd(paste0(BRDYEAR - 1, "1001"))) %>% 
  group_by(beginning_WY) %>% 
  mutate(cum_degree_days = cumsum(degree_days))

rodeo_daily_temps <- rodeo_temps %>% 
  mutate(date = date(timestamp)) %>% 
  group_by(date) %>% 
  summarize(mean_temp = mean(temperature)) %>% 
  ungroup() %>% 
  mutate(degree_days = if_else(mean_temp > temp_threshold, (mean_temp - temp_threshold), 0),
         BRDYEAR = if_else(month(date) > 9, year(date) + 1, year(date)),
         beginning_WY = ymd(paste0(BRDYEAR - 1, "1001"))) %>% 
  group_by(beginning_WY) %>% 
  mutate(cum_degree_days = cumsum(degree_days))

# joining temperature data to breeding data
band_breeding_degree_days <- band_daily_temps %>% 
  right_join(., band_surveys, by = join_by("date" == "breeding_date")) %>% 
  mutate(date = as.POSIXct(date))

write_csv(band_breeding_degree_days, here::here("data", "banducci_degree_days.csv"))

rodeo_breeding_degree_days <- rodeo_daily_temps %>% 
  right_join(., rodeo_surveys, by = join_by("date" == "breeding_date")) %>% 
  mutate(date = as.POSIXct(date))

write_csv(rodeo_breeding_degree_days, here::here("data", "rodeo_degree_days.csv"))

# vectors of first breeding for vlines in plots
band_first_breeding_dates <- as.vector(as.POSIXct(band_surveys$breeding_date))
rodeo_first_breeding_dates <- as.vector(as.POSIXct(rodeo_surveys$breeding_date))

# making plots
band_plot <- ggplot(data = band_temps_with_breeding, aes(x = date, y = mean_temp)) + 
  geom_line() +
  geom_vline(xintercept = as.POSIXct(band_first_breeding_dates), color = "blue") +
  scale_x_datetime(limits = c(as.POSIXct('2015-10-01'), as.POSIXct('2024-10-01')))
band_plot

rodeo_plot <- ggplot(data = rodeo_temps_with_breeding, aes(x = date, y = mean_temp)) + 
  geom_line() +
  geom_vline(xintercept = as.POSIXct(rodeo_first_breeding_dates), color = "blue") +
  scale_x_datetime(limits = c(as.POSIXct('2014-10-01'), as.POSIXct('2024-10-01')))
rodeo_plot
