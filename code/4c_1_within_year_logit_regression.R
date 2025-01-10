library(tidyverse)
library(here)
library(lme4) 
library(gratia)
library(ggplot2)
library(cowplot)
library(mgcv)

#### prepping data for analysis ####
setwd(here::here("code"))
#rename file
onset_of_breeding <- read_csv(here::here("data", "onset_of_breeding.csv"))
breeding_timing <- read_csv(here::here("data", "breeding_timing.csv"))

# scaling covariates
scaled_onset <- onset_of_breeding %>% 
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
  select(-NumberofEggMasses) %>% 
  mutate(complete_case = complete.cases(.)) %>% 
  filter(complete_case == TRUE)

scaled_timing <- breeding_timing %>% 
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
  mutate(complete_case = complete.cases(.)) %>% 
  filter(complete_case == TRUE)

### CHOOSE ONSET OR TIMING ###
within_year <- scaled_onset
# OR
within_year <- scaled_timing

#### LOGISTIC REGRESSION ####

logit_within_year <- glmer(breeding_status ~ 
                     rain_to_date_scaled +
                     cum_sun_hours_scaled +
                     (1 | LocationID) +
                     (1 | BRDYEAR),
                   data = within_year,
                   family = "binomial")
summary(logit_within_year)
plot(logit_within_year)

#model diagnostics
residuals_df <- data.frame(
  fitted = predict(logit_within_year, type = "response"),
  residuals = residuals(logit_within_year, type = "pearson")
)
qqnorm(residuals(logit_within_year, type = "pearson"))


# trying a GAM for logistic regression because the glmer didn't have good model diagnostics
logit_within_year_gam <- gam(breeding_status ~
                         s(rain_to_date_scaled) +
                         s(cum_sun_hours_scaled) +
                         s(LocationID, bs = "re"),
                       data = within_year,
                       family=  "binomial")
summary(logit_within_year_gam)
plot(logit_within_year_gam)

appraise(logit_within_year_gam)
gam.check(logit_within_year_gam)

# plots
ggplot(within_year, aes(x = rain_to_date_scaled, y = breeding_status)) +
  geom_point() +
  stat_smooth(method = "glm", se = FALSE, method.args = list(family = binomial)) +
  theme_bw()

ggplot(within_year, aes(x = cum_sun_hours_scaled, y = breeding_status)) +
  geom_point() +
  stat_smooth(method = "glm", se = FALSE, method.args = list(family = binomial)) +
  theme_bw()
