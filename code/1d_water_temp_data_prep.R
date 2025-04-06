library(dplyr)
library(here)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(stringr)
library(mgcv)
library(lme4)

setwd(here::here("code"))

temp_threshold <- 9

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
         day_number = as.numeric(date - beginning_WY),
         degree_days = if_else(mean_temp > temp_threshold, (mean_temp - temp_threshold), 0)) %>% 
  group_by(BRDYEAR) %>% 
  mutate(cum_degree_days = cumsum(degree_days)) %>% 
  ungroup()
rodeo_temps <- read_csv(here::here("data", "watertemp_rodeo.csv")) %>% 
  mutate(timestamp = mdy_hm(timestamp),
         temperature = if_else(value > 30, ((value - 32) * (5/9)), (value))) %>% 
  mutate(date = date(timestamp)) %>% 
  group_by(date) %>% 
  summarize(mean_temp = mean(temperature)) %>% 
  ungroup() %>% 
  select(date, mean_temp) %>% 
  mutate(BRDYEAR = if_else(month(date) > 9, year(date) + 1, year(date)),
         beginning_WY = ymd(paste0(BRDYEAR - 1, "1001")),
         day_number = as.numeric(date - beginning_WY),
         degree_days = if_else(mean_temp > temp_threshold, (mean_temp - temp_threshold), 0)) %>% 
  group_by(BRDYEAR) %>% 
  mutate(cum_degree_days = cumsum(degree_days)) %>% 
  ungroup()


# importing survey data
band_surveys <- read_csv(here::here("data", "onset_of_breeding.csv")) %>% 
  filter(LocationID == "RC07") %>% 
  mutate(beginning_WY = paste0(BRDYEAR-1, "1001"),
         beginning_WY = ymd(beginning_WY))
rodeo_surveys <- read_csv(here::here("data", "onset_of_breeding.csv")) %>% 
  filter(LocationID == "RL02") %>% 
  mutate(beginning_WY = paste0(BRDYEAR-1, "1001"),
         beginning_WY = ymd(beginning_WY))


# joining temperature data to breeding data
band_breeding_degree_days <- band_temps %>% 
  select(date, mean_temp, cum_degree_days) %>% 
  inner_join(., band_surveys, by = c("date" = "breeding_date")) %>% 
  mutate(date = as.POSIXct(date))
write_csv(band_breeding_degree_days, here::here("data", "banducci_degree_days.csv"))

rodeo_breeding_degree_days <- rodeo_temps %>% 
  select(date, mean_temp, cum_degree_days) %>% 
  inner_join(., rodeo_surveys, by = c("date" = "breeding_date")) %>% 
  mutate(date = as.POSIXct(date))
write_csv(rodeo_breeding_degree_days, here::here("data", "rodeo_degree_days.csv"))

breeding_degree_days_combined <- rbind(band_breeding_degree_days, rodeo_breeding_degree_days)

#### adding GAM model here because I don't wanna make a separate file for it ####

# scaling covariates

scaled_within_year_dd <- breeding_degree_days_combined %>% 
  mutate(
    BRDYEAR_scaled = as.vector(scale(BRDYEAR)),
    yearly_rain_scaled = as.vector(scale(yearly_rain)),
    rain_to_date_scaled = as.vector(scale(rain_to_date)),
    water_flow = as.factor(water_flow),
    water_regime = as.factor(water_regime), 
    cum_degree_days_scaled = as.vector(scale(cum_degree_days)),
    mean_temp_scaled = as.vector(scale(mean_temp)))

# glm model
within_year_glm <- lmer(first_breeding ~ 
                         rain_to_date_scaled +
                         BRDYEAR_scaled + 
                         cum_degree_days_scaled +
                         (1 | LocationID),
                       data = scaled_within_year_dd)
summary(within_year_glm)
plot(within_year_glm)

ggplot(scaled_within_year_band, aes(x = first_breeding, y = cum_degree_days_scaled)) + geom_point()


ggplot(scaled_within_year_band, aes(x = first_breeding, y = cum_degree_days_scaled)) + geom_point()

# plots disregarding year -- proportion of breeding for various thresholds
band_degree_day_plot <- ggplot(data = band_breeding_degree_days, aes(x = cum_degree_days, y=  proportion_breeding)) +
  geom_step(aes(colour = temperature_threshold))

rodeo_degree_day_plot <- ggplot(data = rodeo_breeding_degree_days, aes(x = cum_degree_days, y=  proportion_breeding)) +
  geom_step(aes(colour = temperature_threshold))

# plots for degree days vs. day of breeding year
band_temp_plot <- ggplot(data = band_breeding_degree_days %>% filter(BRDYEAR > 2015, temperature_threshold == 9), aes(y = cum_degree_days, x = first_breeding)) +
  geom_point(aes(colour = as.factor(BRDYEAR))) +
  geom_smooth(method = "lm", color = "black", alpha = 0.2)
band_temp_plot

rodeo_temp_plot <- ggplot(data = rodeo_breeding_degree_days %>% filter(BRDYEAR > 2015, temperature_threshold == 9), aes(y = cum_degree_days, x = first_breeding)) +
  geom_point(aes(colour = as.factor(BRDYEAR))) +
  geom_smooth(method = "lm", color = "black", alpha = 0.2)
rodeo_temp_plot

# comparison plot for water temp (on the day of breeding) vs. breeding date
water_temp_plot <- ggplot(data = rbind(band_breeding_degree_days, rodeo_breeding_degree_days), aes(x = first_breeding, y = mean_temp)) +
  geom_point(aes(color = LocationID))

# comparison plot for degree days
degree_day_plot <- ggplot(data = rbind(band_breeding_degree_days, rodeo_breeding_degree_days) %>% filter(temperature_threshold == 9), aes(x = first_breeding, y = cum_degree_days)) +
  geom_point(aes(color = LocationID))

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
