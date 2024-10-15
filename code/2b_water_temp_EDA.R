library(dplyr)
library(here)
library(tidyverse)
library(ggplot2)
library(lubridate)

setwd(here::here("code"))

# reading in water temp data
band_temps <- read_csv(here::here("data", "banducci_degree_days.csv"))
rodeo_temps <- read_csv(here::here("data", "rodeo_degree_days.csv"))