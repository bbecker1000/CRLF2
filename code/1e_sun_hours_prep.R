library(tidyverse)
library(dplyr)
library(here)
library(lubridate)
library(ggplot2)

# reading in data from spreadsheets
setwd(here::here("code"))

sun_hours <- read_csv(here::here("data", "SunHours.csv")) %>% 
  mutate(str_time = substr(str_time, 1, 10),
         str_time = ymd(str_time))

sun_hours_plot <- ggplot(sun_hours, aes(x = str_time, y = dir_dur)) +
  geom_line(aes(color = LocationID))
