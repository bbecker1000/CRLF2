bprior <- c(
  # yearly rain (positive)
  # prior(student_t(1, 0.5, 0.5), class = b, coef = syearly_rain_scaled_1),
  
  # yearly rain interactions (positive but more so for seasonal)
  #prior(student_t(1, 0, 0.5), coef = yearly_rain_scaled:water_regime),
  prior(student_t(1, 0.5, 0.5), coef =  yearly_rain_scaled),
  
  # yearly rain interactions -- hurdle. not working and idk how to run get_prior() properly for them
  # I think we want these to be negative since we're measuring the probability of zero eggs
  # prior(student_t(1, -0.25, 0.5), coef = syearly_rain_scaled:water_regimeperennial_1, dpar = hu),
  # prior(student_t(1, -0.5, 0.5), coef = syearly_rain_scaled:water_regimeseasonal_1, dpar = hu),
  
  # other covariates (feel free to change as you see fit)
  # prior(student_t(1, 0.5, 0.5), coef =  water_regimeseasonal), # slightly positive based on hypotheses
  prior(student_t(1, 0.5, 0.5), coef =  water_flowlentic), # slightly positive based on hypotheses
  prior(student_t(1, -0.5, 0.5), coef =  water_flowlotic), # slightly negtive based on hypotheses
  prior(student_t(1, -0.5, 0.5), coef =  interpolated_canopy_scaled), # slightly negative based on hypotheses
  prior(student_t(1, 0.25, 0.5), coef =  mean_percent_sub_scaled), # slightly positive? based on hypotheses
  prior(student_t(1, 0, 0.5), coef =  BRDYEAR_scaled),
  #prior(student_t(1, -0.25, 0.5), coef =  mean_percent_water_scaled),
  prior(student_t(1, 0, 0.5), coef =  WaterTemp_scaled)
)






t0 <- Sys.time()

hurdle_AirTempdata.linear <- brm(
  bf(num_egg_masses ~ 
       BRDYEAR_scaled + 
       mean_percent_water_scaled + 
       interpolated_canopy_scaled +
       WaterTemp_scaled +  
       mean_percent_sub_scaled +
       yearly_rain_scaled*water_regime +
       water_flow +
       #(1 | Watershed/LocationID),
     (BRDYEAR_scaled | Watershed/LocationID),
     hu ~ 
       yearly_rain_scaled*water_regime +      # inflated model for zeros
       #(1|Watershed/LocationID)
     (yearly_rain_scaled:water_regime | Watershed/LocationID)
     ),
  data = complete_scaled_btw,
  family = hurdle_negbinomial(),
  prior = bprior,
  chains = 3, cores = 3,
  iter = 2000, # only need about 1000 for inference (3500-2500 warmup = 1000)
  warmup = 1500, 
  control = list(adapt_delta = 0.98)
)

t1 <- Sys.time()
t1-t0
beepr::beep(0)


summary(hurdle_AirTempdata.linear)
conditional_effects(hurdle_AirTempdata.linear, surface = FALSE, prob = 0.8)
conditional_effects(hurdle_AirTempdata.linear, dpar = "hu")



bprior <- c(
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

#random slopes
#### hurdle model ####
##### without Airtemp #####
t0 <- Sys.time()

hurdle_AirTempdata.slopes <- brm(
  bf(num_egg_masses ~ 
       s(BRDYEAR_scaled) + 
       s(mean_percent_water_scaled) + 
       s(interpolated_canopy_scaled) +
       s(WaterTemp_scaled) +  
       s(mean_percent_sub_scaled) +
       s(yearly_rain_scaled, by =water_regime) +
       (water_flow) +
       Watershed +
       (yearly_rain_scaled:water_regime | Watershed/LocationID),
     hu ~ 
       s(yearly_rain_scaled, by = water_regime) +    
       Watershed + # inflated model for zeros
       (yearly_rain_scaled:water_regime | Watershed/LocationID)
     ),
  data = complete_scaled_btw,
  family = hurdle_negbinomial(),
  prior = bprior,
  chains = 3, cores = 3,
  iter = 5000, # only need about 1000 for inference (3500-2500 warmup = 1000)
  warmup = 500, 
  control = list(adapt_delta = 0.98)
)
t1 <- Sys.time()
t1-t0
beepr::beep(0)


summary(hurdle_AirTempdata.slopes)
conditional_effects(hurdle_AirTempdata.slopes, surface = FALSE, prob = 0.8)
conditional_effects(hurdle_AirTempdata.slopes, dpar = "hu")





