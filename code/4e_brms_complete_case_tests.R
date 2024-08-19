# comparing removing covariates from complete case filter
# testing model with added depth covariate

library(brms)
library(tidyverse)

#### complete cases ####
btw_data <- read_csv(here::here("data", "between_year_data.csv"))  %>% 
  mutate(LocationID = as.factor(LocationID),
         Watershed = as.factor(Watershed), 
         LocationInWatershed = interaction(Watershed, LocationID))

scaled_btw <- btw_data %>% 
  mutate(
    BRDYEAR_scaled = as.vector(scale(BRDYEAR)),
    mean_percent_sub_scaled = as.vector(scale(mean_percent_sub)),
    mean_percent_emerg_scaled = as.vector(scale(mean_percent_emerg)),
    mean_percent_water_scaled = as.vector(scale(mean_percent_water)),
    interpolated_canopy_scaled = as.vector(scale(interpolated_canopy)),
    yearly_rain_scaled = as.vector(scale(yearly_rain)),
    mean_max_depth_scaled = as.vector(scale(mean_max_depth)),
    max_depth_scaled = as.vector(scale(max_depth)),
    # AirTemp_scaled = as.vector(scale(AirTemp)),
    WaterTemp_scaled = as.vector(scale(WaterTemp))
    # mean_salinity_scaled = as.vector(scale(mean_salinity)),
    # max_salinity_scaled = as.vector(scale(max_salinity)),
  ) %>% 
  select(-BRDYEAR,-mean_max_depth, -max_depth, -AirTemp, -WaterTemp, -mean_percent_emerg, -mean_percent_sub, -mean_percent_water, -interpolated_canopy, -yearly_rain, -CoastalSite) #remove non-scaled variables

# creating a "complete case" column
scaled_btw$complete_case <- complete.cases(scaled_btw)
complete_scaled_btw <- scaled_btw %>% 
  filter(complete_case == TRUE) %>% 
  select(-complete_case)

#### priors ####
bprior.no.sal <- c(
  # yearly rain (positive)
  # prior(student_t(1, 0.5, 0.5), class = b, coef = syearly_rain_scaled_1),
  
  # yearly rain interactions (positive but more so for seasonal)
  prior(student_t(1, 0, 0.5), coef = syearly_rain_scaled:water_regimeperennial_1),
  prior(student_t(1, 0.5, 0.5), coef =  syearly_rain_scaled:water_regimeseasonal_1),
  
  # yearly rain interactions -- hurdle. not working and idk how to run get_prior() properly for them
  # I think we want these to be negative since we're measuring the probability of zero eggs
  # prior(student_t(1, -0.25, 0.5), coef = syearly_rain_scaled:water_regimeperennial_1, dpar = hu),
  # prior(student_t(1, -0.5, 0.5), coef = syearly_rain_scaled:water_regimeseasonal_1, dpar = hu),
  
  # other covariates (feel free to change as you see fit)
  # prior(student_t(1, 0.5, 0.5), coef =  water_regimeseasonal), # slightly positive based on hypotheses
  prior(student_t(1, 0.5, 0.5), coef =  water_flowlentic), # slightly positive based on hypotheses
  prior(student_t(1, -0.5, 0.5), coef =  water_flowlotic), # slightly negtive based on hypotheses
  prior(student_t(1, -0.5, 0.5), coef =  sinterpolated_canopy_scaled_1), # slightly negative based on hypotheses
  prior(student_t(1, 0.25, 0.5), coef =  smean_percent_sub_scaled_1), # slightly positive? based on hypotheses
  prior(student_t(1, 0, 0.5), coef =  sBRDYEAR_scaled_1 ),
  prior(student_t(1, -0.25, 0.5), coef =  smean_percent_water_scaled_1),
  prior(student_t(1, 0, 0.5), coef =  sWaterTemp_scaled_1)
)


#### hurdle model ####
##### without Airtemp #####
t0 <- Sys.time()

hurdle_AirTempdata <- brm(
  bf(num_egg_masses ~ 
       s(BRDYEAR_scaled) + 
       s(mean_percent_water_scaled) + 
       s(interpolated_canopy_scaled) +
       s(WaterTemp_scaled) +  
       s(mean_percent_sub_scaled) +
       s(yearly_rain_scaled, by =water_regime) +
       (water_flow) +
       (1 | Watershed/LocationID),
     hu ~ 
       s(yearly_rain_scaled, by = water_regime) +      # inflated model for zeros
       (1|Watershed/LocationID)),
  data = complete_scaled_btw,
  family = hurdle_negbinomial(),
  prior = bprior.no.sal,
  chains = 3, cores = 3,
  iter = 11500, # only need about 1000 for inference (3500-2500 warmup = 1000)
  warmup = 10500, 
  control = list(adapt_delta = 0.99)
)

##### without air tmp, with mean_max_depth #####
t0 <- Sys.time()
hurdle_depth <- brm(
  bf(num_egg_masses ~ 
       s(BRDYEAR_scaled) + 
       s(mean_percent_water_scaled) + 
       s(interpolated_canopy_scaled) +
       s(WaterTemp_scaled) +  
       s(mean_percent_sub_scaled) +
       s(mean_max_depth_scaled) +
       s(yearly_rain_scaled, by =water_regime) +
       (water_flow) +
       (1 | Watershed/LocationID),
     hu ~ 
       s(yearly_rain_scaled, by = water_regime) +      # inflated model for zeros
       (1|Watershed/LocationID)),
  data = complete_scaled_btw,
  family = hurdle_negbinomial(),
  prior = bprior.no.sal,
  chains = 3, cores = 3,
  iter = 15500, # only need about 1000 for inference (3500-2500 warmup = 1000)
  warmup = 14500, 
  control = list(adapt_delta = 0.99)
)



#### plots ####
t1 <- Sys.time()
t1-t0
beepr::beep(0)

hurdle_test <- hurdle_AirTempdata 
hurdle_test <- hurdle_depth

summary(hurdle_test)
conditional_effects(hurdle_test, surface = FALSE, prob = 0.8)
conditional_effects(hurdle_test, dpar = "hu")

