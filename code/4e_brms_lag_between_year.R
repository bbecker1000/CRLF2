library(brms)
library(tidyverse)

##### prep lagged data set ####
# see EDA.rmd for exploration and selection of continuous sites
lag_between_year <- read_csv(here::here("data", "lag_between_year_data.csv")) %>% 
  filter(LocationID=='KC01' | LocationID=='KC02' |LocationID=='KC03' | LocationID=='LS01' | LocationID=='LS04' | LocationID=='LS05' | LocationID=='LS06' | LocationID=='LS07' | LocationID=='LS08' | LocationID=='LS09' | LocationID=='LS11' | LocationID=='LS12' | LocationID=='MC01' | LocationID=='RC01' | LocationID=='RC02' | LocationID=='RC03' | LocationID=='RC07' | LocationID=='RC10' | LocationID=='RC11' | LocationID=='RC13' | LocationID=='RC14' | LocationID=='RC15' | LocationID=='RC17' | LocationID=='RC18' | LocationID=='RC20' | LocationID=='RC24' | LocationID=='RC25' | LocationID=='RC26' | LocationID=='RL02' | LocationID=='RL04' | LocationID=='RL07' | LocationID=='TV02' | LocationID=='TV03' | LocationID=='TV06' | LocationID=='WG01')

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

#### priors ####
# taken 10-16-2024 from bprior.no.sal.linear.zi priors (4d)
lag.priors <- c(
  prior(normal(0, 0.5), coef = BRDYEAR_scaled),
  prior(normal(-0.5, 0.5), coef = interpolated_canopy_scaled),
  prior(normal(0.25, 0.5), coef = mean_percent_sub_scaled),
  prior(normal(0.0, 0.5), coef = mean_percent_water_scaled), # add squared term
  prior(normal(0.5, 1), coef = water_flowlentic),
  prior(normal(-0.5, 1), coef = water_flowlotic),
  prior(normal(0.5, 1), coef = water_regimeseasonal),
  prior(normal(0.5, 1), coef = yearly_rain_lag_scaled),
  prior(normal(0.25, 1), coef = WaterTemp_scaled), # add squared term
  prior(normal(0.25, 1), coef = yearly_rain_scaled),
  prior(normal(0.25, 1), coef = yearly_rain_scaled:water_regimeseasonal),
  prior(normal(0.25, 1), coef = yearly_rain_lag_scaled:water_regimeseasonal),
  prior(normal(0.5, 1), coef = yearly_rain_lag_scaled:num_egg_masses_lag_scaled)
)

#### model ####
lag.zi.linear <- brm(
  num_egg_masses ~ 
    BRDYEAR_scaled + 
    mean_percent_water_scaled + 
    interpolated_canopy_scaled +
    WaterTemp_scaled +  
    mean_percent_sub_scaled +
    yearly_rain_scaled +
    yearly_rain_scaled : water_regime +
    yearly_rain_lag_scaled +
    yearly_rain_lag_scaled:num_egg_masses_lag_scaled +
    num_egg_masses_lag_scaled +
    water_regime +
    yearly_rain_lag_scaled : water_regime +
    water_flow +
    proportion_high_water_vis +
    proportion_na_water_vis +
    (1 | Watershed/LocationID),
  data = scaled_lag_between_year,  # run 4a to prep data file)
  family = zero_inflated_negbinomial(),
  prior = lag.priors,
  chains = 3, cores = 3,
  iter = 5500, # 11500, # only need about 500 for inference
  warmup = 4500, #11000, 
  control = list(adapt_delta = 0.97))

summary(lag.zi.linear, prob = 0.89)

#### plots ####