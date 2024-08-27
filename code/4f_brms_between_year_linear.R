

t0 <- Sys.time()

hurdle_AirTempdata.linear <- brm(
  bf(num_egg_masses ~ 
       BRDYEAR_scaled + 
       mean_percent_water_scaled + 
       interpolated_canopy_scaled +
       WaterTemp_scaled +  
       mean_percent_sub_scaled +
       yearly_rain_scaled:water_regime +
       water_flow +
       (1 | Watershed/LocationID),
     hu ~ 
       yearly_rain_scaled:water_regime +      # inflated model for zeros
       (1|Watershed/LocationID)),
  data = complete_scaled_btw,
  family = hurdle_negbinomial(),
  #prior = bprior.no.sal,
  chains = 3, cores = 3,
  iter = 2000, # only need about 1000 for inference (3500-2500 warmup = 1000)
  warmup = 1500, 
  control = list(adapt_delta = 0.95)
)

t1 <- Sys.time()
t1-t0
beepr::beep(0)