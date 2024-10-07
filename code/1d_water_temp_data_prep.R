library(dplyr)
library(here)
library(tidyverse)
library(ggplot2)
library(lubridate)

setwd(here::here("code"))

band_temps <- read_csv(here::here("data", "watertemp_banducci.csv")) %>% 
  mutate(timestamp = mdy_hm(timestamp),
         temperature = (value - 32) * (5/9)) %>% 
  select(timestamp, temperature)
rodeo_temps <- read_csv(here::here("data", "watertemp_rodeo.csv")) %>% 
  mutate(timestamp = mdy_hm(timestamp),
         temperature = if_else(value > 30, ((value - 32) * (5/9)), (value))) %>% 
  select(timestamp, temperature)
band_surveys <- read_csv(here::here("data", "onset_of_breeding.csv")) %>% 
  filter(LocationID == "RC07") %>% 
  mutate(beginning_WY = paste0(BRDYEAR-1, "1001"),
         beginning_WY = ymd(beginning_WY),
         breeding_date = beginning_WY + ddays(first_breeding))
rodeo_surveys <- read_csv(here::here("data", "onset_of_breeding.csv")) %>% 
  filter(LocationID == "RL02") %>% 
  mutate(beginning_WY = paste0(BRDYEAR-1, "1001"),
         beginning_WY = ymd(beginning_WY),
         breeding_date = beginning_WY + ddays(first_breeding))

band_daily_temps <- band_temps %>% 
  mutate(date = date(timestamp)) %>% 
  group_by(date) %>% 
  summarize(mean_temp = mean(temperature)) %>% 
  ungroup()

band_temps_with_breeding <- band_daily_temps %>% 
  full_join(., band_surveys, by = join_by("date" == "breeding_date")) %>% 
  mutate(date = as.POSIXct(date))

band_first_breeding_dates <- as.vector(as.POSIXct(band_surveys$breeding_date))

band_plot <- ggplot(data = band_temps_with_breeding, aes(x = date, y = mean_temp)) + 
  geom_line() +
  geom_vline(xintercept = as.POSIXct(band_first_breeding_dates), color = "blue") +
  scale_x_datetime(limits = c(as.POSIXct('2015-10-01'), as.POSIXct('2024-10-01')))

band_plot

rodeo_daily_temps <- rodeo_temps %>% 
  mutate(date = date(timestamp)) %>% 
  group_by(date) %>% 
  summarize(mean_temp = mean(temperature)) %>% 
  ungroup()

rodeo_temps_with_breeding <- rodeo_daily_temps %>% 
  full_join(., rodeo_surveys, by = join_by("date" == "breeding_date")) %>% 
  mutate(date = as.POSIXct(date))

rodeo_first_breeding_dates <- as.vector(as.POSIXct(rodeo_surveys$breeding_date))

rodeo_plot <- ggplot(data = rodeo_temps_with_breeding, aes(x = date, y = mean_temp)) + 
  geom_line() +
  geom_vline(xintercept = as.POSIXct(rodeo_first_breeding_dates), color = "blue") +
  scale_x_datetime(limits = c(as.POSIXct('2014-10-01'), as.POSIXct('2024-10-01')))

rodeo_plot
