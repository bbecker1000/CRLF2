library(dplyr)
library(here)
library(tidyverse)
library(ggplot2)
library(lubridate)

setwd(here::here("code"))

band_temps <- read_csv(here::here("data", "watertemp_banducci.csv")) %>% 
  mutate(timestamp = mdy_hm(timestamp))
rodeo_temps <- read_csv(here::here("data", "watertemp_rodeo.csv")) %>% 
  mutate(timestamp = mdy_hm(timestamp))

band_temps_at_noon <- band_temps %>% 
  filter(hour(timestamp) == 12 & minute(timestamp) == 00)

rodeo_temps_at_noon <- rodeo_temps %>% 
  filter(hour(timestamp) == 12 & minute(timestamp) == 00)

band_plot <- ggplot(band_temps_at_noon, aes(x = timestamp, y = value)) + geom_point()

rodeo_plot <- ggplot(rodeo_temps_at_noon, aes(x = timestamp, y = value)) + geom_point()

