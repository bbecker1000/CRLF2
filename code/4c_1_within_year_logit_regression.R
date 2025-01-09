library(tidyverse)
library(mgcv)
library(here)
library(nlme) 
library(gratia)
library(ggplot2)
library(cowplot)

#### prepping data for analysis ####
setwd(here::here("code"))
#rename file
onset_of_breeding <- read_csv(here::here("data", "onset_of_breeding.csv"))

# scaling covariates
scaled_within_year <- onset_of_breeding %>% 
  mutate(
    BRDYEAR_scaled = as.vector(scale(BRDYEAR)),
    yearly_rain_scaled = as.vector(scale(yearly_rain)),
    rain_to_date_scaled = as.vector(scale(rain_to_date)),
    water_flow = as.factor(water_flow),
    water_regime = as.factor(water_regime), 
    Watershed = as.factor(Watershed),
    LocationID = as.factor(LocationID),
    cum_sun_hours_scaled = as.vector(scale(cum_sun_hours)),
    dir_dur_scaled = as.vector(scale(dir_dur))) %>% 
  select(-NumberofEggMasses)

# creating a "complete case" column
scaled_within_year$complete_case <- complete.cases(scaled_within_year)
complete_onset <- scaled_within_year %>% filter(complete_case == TRUE) %>% select(-complete_case)

#### LOGISTIC REGRESSION ####

logit_onset <- glm(breeding_status ~ 
                     rain_to_date_scaled +
                     cum_sun_hours_scaled,
                   data = complete_onset,
                   family = "binomial")
summary(logit_onset)
plot(logit_onset)

ggplot(complete_onset, aes(x = rain_to_date_scaled, y = breeding_status)) +
  geom_point() +
  stat_smooth(method = "glm", se = FALSE, method.args = list(family = binomial)) +
  theme_bw()

ggplot(complete_onset, aes(x = cum_sun_hours_scaled, y = breeding_status)) +
  geom_point() +
  stat_smooth(method = "glm", se = FALSE, method.args = list(family = binomial)) +
  theme_bw()
