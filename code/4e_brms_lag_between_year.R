library(brms)
library(tidyverse)

##### prep lag data set ####
# see EDA.rmd for exploration and selection of continuous sites
lag_between_year <- read_csv(here::here("data", "lag_between_year_data.csv")) %>% filter(LocationID=='KC01' | LocationID=='KC02' |LocationID=='KC03' | LocationID=='LS01' | LocationID=='LS04' | LocationID=='LS05' | LocationID=='LS06' | LocationID=='LS07' | LocationID=='LS08' | LocationID=='LS09' | LocationID=='LS11' | LocationID=='LS12' | LocationID=='MC01' | LocationID=='RC01' | LocationID=='RC02' | LocationID=='RC03' | LocationID=='RC07' | LocationID=='RC10' | LocationID=='RC11' | LocationID=='RC13' | LocationID=='RC14' | LocationID=='RC15' | LocationID=='RC17' | LocationID=='RC18' | LocationID=='RC20' | LocationID=='RC24' | LocationID=='RC25' | LocationID=='RC26' | LocationID=='RL02' | LocationID=='RL04' | LocationID=='RL07' | LocationID=='TV02' | LocationID=='TV03' | LocationID=='TV06' | LocationID=='WG01')

# complete cases
## creating a "complete case" column
lag_between_year$complete_case <- complete.cases(lag_between_year)
complete_lag_btw_data <- lag_between_year %>% filter(complete_case == TRUE)

## write to CSV
write_csv(complete_lag_btw_data, here::here("data", "complete_lag_btw_data.csv"))

## scaling covariates
scaled_lag_between_year <- complete_lag_btw_data %>% 
  mutate(
    BRDYEAR_scaled = as.vector(scale(BRDYEAR)),
    mean_percent_sub_scaled = as.vector(scale(mean_percent_sub)),
    mean_percent_emerg_scaled = as.vector(scale(mean_percent_emerg)),
    mean_percent_water_scaled = as.vector(scale(mean_percent_water)),
    interpolated_canopy_scaled = as.vector(scale(interpolated_canopy)),
    yearly_rain_scaled = as.vector(scale(yearly_rain)),
    mean_max_depth_scaled = as.vector(scale(mean_max_depth)),
    max_depth_scaled = as.vector(scale(max_depth)),
    AirTemp_scaled = as.vector(scale(AirTemp)),
    WaterTemp_scaled = as.vector(scale(WaterTemp)),
    yearly_rain_lag_scaled = as.vector(scale(yearly_rain_lag)),
    num_egg_masses_lag_scaled = as.vector(scale(num_egg_masses_lag)))

