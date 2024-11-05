library(tidyverse)
library(dplyr)
library(here)
library(lubridate)
library(ggplot2)

# reading in data from spreadsheets
setwd(here::here("code"))

sun_hours <- read_csv(here::here("data", "Sun_Hours_WY2024.csv")) %>% 
  mutate(str_time = substr(str_time, 1, 10),
         str_time = ymd(str_time))

sun_hours_manual <- read_csv(here::here("data", "sun_hours_manual.csv")) %>% 
  select(-StartLat, -StartLon, -LatLong) %>% 
  rename(day_length_20241221 = day_length_winter_solstice,
         day_length_20240620 = day_length_summer_solstice,
         day_length_20240319 = day_length_spring_equinox) %>% 
  pivot_longer(cols = starts_with("day_length_"), names_to = "date", values_to = "length") %>% 
  mutate(date = substr(date, 12, 20),
         date = ymd(date))

sun_hours_plot <- ggplot(sun_hours, aes(x = str_time, y = dir_dur)) +
  geom_line(aes(color = LocationID))

sun_hours_manual_plot <- ggplot(sun_hours_manual, aes(x = date, y = length)) +
  geom_smooth(aes(color = LocationID), linewidth = 0.5)
