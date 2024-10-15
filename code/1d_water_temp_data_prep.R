library(dplyr)
library(here)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(stringr)

setwd(here::here("code"))

# importing water temperature data
band_temps <- read_csv(here::here("data", "watertemp_banducci.csv")) %>% 
  mutate(timestamp = mdy_hm(timestamp),
         temperature = (value - 32) * (5/9)) %>% 
  mutate(date = date(timestamp)) %>% 
  group_by(date) %>% 
  summarize(mean_temp = mean(temperature)) %>% 
  ungroup() %>% 
  select(date, mean_temp) %>% 
  mutate(BRDYEAR = if_else(month(date) > 9, year(date) + 1, year(date)),
         beginning_WY = ymd(paste0(BRDYEAR - 1, "1001")),
         day_number = as.numeric(date - beginning_WY))
rodeo_temps <- read_csv(here::here("data", "watertemp_rodeo.csv")) %>% 
  mutate(timestamp = mdy_hm(timestamp),
         temperature = if_else(value > 30, ((value - 32) * (5/9)), (value))) %>% 
  mutate(date = date(timestamp)) %>% 
  group_by(date) %>% 
  summarize(mean_temp = mean(temperature)) %>% 
  ungroup() %>% 
  select(date, mean_temp)


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

for(temp_threshold in 7:11) {
  band_temps <- band_temps %>% 
    mutate(!!paste0("degree_days_", temp_threshold) := if_else(mean_temp > temp_threshold, (mean_temp - temp_threshold), 0)) %>% 
    group_by(beginning_WY) %>% 
    mutate(!!paste0("cum_degree_days_", temp_threshold) := cumsum(!!sym(paste0("degree_days_", temp_threshold)))) %>% 
    ungroup()
}

for(temp_threshold in 7:11) {
  rodeo_temps <- rodeo_temps %>% 
    mutate(!!paste0("degree_days_", temp_threshold) := if_else(mean_temp > temp_threshold, (mean_temp - temp_threshold), 0),
           BRDYEAR = if_else(month(date) > 9, year(date) + 1, year(date)),
           beginning_WY = ymd(paste0(BRDYEAR - 1, "1001"))) %>% 
    group_by(beginning_WY) %>% 
    mutate(!!paste0("cum_degree_days_", temp_threshold) := cumsum(!!sym(paste0("degree_days_", temp_threshold)))) %>% 
    ungroup()
}

# joining temperature data to breeding data
band_breeding_degree_days <- band_temps %>% 
  pivot_longer(cols = starts_with("cum_degree_days"), names_to = "temperature_threshold", values_to = "cum_degree_days") %>%
  mutate(temperature_threshold = str_extract(temperature_threshold, "(?<=_)[^_]+$"),
         temperature_threshold = factor(temperature_threshold, levels = c('7', '8', '9', '10', '11'))) %>% 
  select(date, mean_temp, cum_degree_days, temperature_threshold) %>% 
  inner_join(., band_surveys, by = join_by("date" == "breeding_date")) %>% 
  mutate(date = as.POSIXct(date)) %>% 
  group_by(temperature_threshold) %>% 
  arrange(cum_degree_days) %>% 
  mutate(proportion_breeding = cumsum(NumberofEggMasses) / sum(NumberofEggMasses))

write_csv(band_breeding_degree_days, here::here("data", "banducci_degree_days.csv"))

rodeo_breeding_degree_days <- rodeo_temps %>% 
  pivot_longer(cols = starts_with("cum_degree_days"), names_to = "temperature_threshold", values_to = "cum_degree_days") %>%
  mutate(temperature_threshold = str_extract(temperature_threshold, "(?<=_)[^_]+$"),
         temperature_threshold = factor(temperature_threshold, levels = c('7', '8', '9', '10', '11'))) %>% 
  select(date, mean_temp, cum_degree_days, temperature_threshold) %>% 
  inner_join(., band_surveys, by = join_by("date" == "breeding_date")) %>% 
  mutate(date = as.POSIXct(date)) %>% 
  group_by(temperature_threshold) %>% 
  arrange(cum_degree_days) %>% 
  mutate(proportion_breeding = cumsum(NumberofEggMasses) / sum(NumberofEggMasses))

write_csv(rodeo_breeding_degree_days, here::here("data", "rodeo_degree_days.csv"))

# plots disregarding year -- proportion of breeding for verious thresholds
band_degree_day_plot <- ggplot(data = band_breeding_degree_days, aes(x = cum_degree_days, y=  proportion_breeding)) +
  geom_step(aes(colour = temperature_threshold))

rodeo_degree_day_plot <- ggplot(data = rodeo_breeding_degree_days, aes(x = cum_degree_days, y=  proportion_breeding)) +
  geom_step(aes(colour = temperature_threshold))

# plots for degree days vs. day of breeding year, colored by year
# TODO: make this plot but for rodeo. also maybe figure out how to display this in a way that's easier to interpret...
band_temp_plot <- ggplot(data = band_temps %>% filter(BRDYEAR > 2015), aes(y = cum_degree_days_9, x = day_number)) +
  geom_step(aes(colour = as.factor(BRDYEAR))) +
  geom_vline(data = band_surveys %>% filter(BRDYEAR > 2015), aes(xintercept = first_breeding, color = as.factor(BRDYEAR)))

# vectors of first breeding for vlines in plots
band_first_breeding_dates <- as.vector(as.POSIXct(band_surveys$breeding_date))
rodeo_first_breeding_dates <- as.vector(as.POSIXct(rodeo_surveys$breeding_date))

# making plots
band_plot <- ggplot(data = band_breeding_degree_days, aes(x = date, y = cum_degree_days)) + 
  geom_line() +
  geom_vline(xintercept = as.POSIXct(band_first_breeding_dates), color = "blue") +
  scale_x_datetime(limits = c(as.POSIXct('2015-10-01'), as.POSIXct('2024-10-01')))
band_plot

rodeo_plot <- ggplot(data = rodeo_breeding_degree_days, aes(x = date, y = cum_degree_days)) + 
  geom_line() +
  geom_vline(xintercept = as.POSIXct(rodeo_first_breeding_dates), color = "blue") +
  scale_x_datetime(limits = c(as.POSIXct('2014-10-01'), as.POSIXct('2024-10-01')))
rodeo_plot
