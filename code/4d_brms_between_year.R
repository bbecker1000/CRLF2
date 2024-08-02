library(brms)

t0 <- Sys.time()
t0<-Sys.time()

#### BRM model (no longer using) ####
mod.brm <- brm(bf(num_egg_masses ~  #bf creates a model statement for compilation
                      s(BRDYEAR_scaled) + 
                      # s(yearly_rain_scaled) + ## commented out because water_regime added
                      s(mean_percent_water_scaled) + 
                      s(interpolated_canopy_scaled) +
                      s(WaterTemp_scaled) +  
                      max_depth_scaled +
                      max_salinity_scaled * CoastalSite + # not smoothed
                      s(yearly_rain_scaled * water_regime) + # added new co-variate, smoothed?
                      (1|Watershed/LocationID),
               zi ~ s(yearly_rain_scaled) +      # inflated model for zeros
                 (1|Watershed/LocationID)),   # added random effects 
               data = scaled_between_year,
               family = zero_inflated_negbinomial(),
               chains = 2, cores = 2,
               iter = 10000, # needs more iterations with added covariates
               control = list(adapt_delta = 0.98)) #reduce divergences

save(mod.brm, file = "Output/mod.brm.RData")
#load("Output/mod.brm.RData")
t1<-Sys.time()
t1-t0

beepr::beep(0)
  
summary(mod.brm)

#### priors ####

# # getting the names of all the priors -- useful for setting priors, but not needed to run the model
# get_prior(num_egg_masses ~ 
#             s(BRDYEAR_scaled) + 
#             s(mean_percent_water_scaled) + 
#             s(interpolated_canopy_scaled) +
#             s(WaterTemp_scaled) +  
#             s(max_salinity_scaled, by = CoastalSite) + 
#             s(yearly_rain_scaled, by = water_regime) +
#             (water_flow) +
#             (1 | Watershed/LocationID),
#           hu ~ 
#             s(yearly_rain_scaled, by = water_regime) +
#             (1|Watershed/LocationID),
#         data = scaled_between_year,
#         family = hurdle_negbinomial())

# could add more priors if helpful
bprior <- c(
  # yearly rain (positive)
  # prior(student_t(1, 0.5, 0.5), class = b, coef = syearly_rain_scaled_1),

  # yearly rain interactions (positive but more so for seasonal)
  # prior(student_t(1, 0.25, 0.5), coef = syearly_rain_scaled:water_regimeperennial_1),
  # prior(student_t(1, 0.5, 0.5), coef =  syearly_rain_scaled:water_regimeseasonal_1),

  # yearly rain interactions -- hurdle. not working and idk how to run get_prior() properly for them
  # I think we want these to be negative since we're measuring the probability of zero eggs
  # prior(student_t(1, -0.25, 0.5), coef = syearly_rain_scaled:water_regimeperennial_1, dpar = hu),
  # prior(student_t(1, -0.5, 0.5), coef = syearly_rain_scaled:water_regimeseasonal_1, dpar = hu),

  # salinity / coastal site interactions
  # no effect for non-coastal sites since they're all 0, slightly negative for coastal sites
  # prior(student_t(1, 0, 0.5), coef =  smax_salinity_scaled:CoastalSiteFALSE_1),
  # prior(student_t(1, -0.25, 0.5), coef =  smax_salinity_scaled:CoastalSiteTRUE_1),

  # other covariates (feel free to change as you see fit)
  # prior(student_t(1, 0.5, 0.5), coef =  water_regimeseasonal), # slightly positive based on hypotheses
  prior(student_t(1, 0.25, 0.5), coef =  water_flowlentic), # slightly positive based on hypotheses
  prior(student_t(1, -0.25, 0.5), coef =  water_flowlotic), # slightly negtive based on hypotheses
  prior(student_t(1, -0.25, 0.5), coef =  sinterpolated_canopy_scaled_1), # slightly negative based on hypotheses
  prior(student_t(1, 0, 0.5), coef =  sBRDYEAR_scaled_1 ),
  prior(student_t(1, 0, 0.5), coef =  smean_percent_water_scaled_1),
  prior(student_t(1, 0, 0.5), coef =  sWaterTemp_scaled_1)
)


#### hurdle model ####
mod.hurdle <- brm(
  bf(num_egg_masses ~ 
       s(BRDYEAR_scaled) + 
       s(mean_percent_water_scaled) + 
       s(interpolated_canopy_scaled) +
       s(WaterTemp_scaled) +  
       # s(max_salinity_scaled, by = CoastalSite) + 
       # s(yearly_rain_scaled, by =water_regime) +
       (water_flow) +
       (1 | Watershed/LocationID),
     hu ~ 
       # s(yearly_rain_scaled, by = water_regime) +      # inflated model for zeros
       (1|Watershed/LocationID)),
  data = scaled_between_year,
  family = hurdle_negbinomial(),
  prior = bprior,
  chains = 3, cores = 3,
  iter = 40000, # needs more iterations with added covariates
  control = list(adapt_delta = 0.99999, max_treedepth = 20)
)

#### hurdle model Remove salinity ####
bprior.no.sal <- c(
  # yearly rain (positive)
  # prior(student_t(1, 0.5, 0.5), class = b, coef = syearly_rain_scaled_1),
  
  # yearly rain interactions (positive but more so for seasonal)
  prior(student_t(1, 0.25, 0.5), coef = syearly_rain_scaled:water_regimeperennial_1),
  prior(student_t(1, 0.5, 0.5), coef =  syearly_rain_scaled:water_regimeseasonal_1),
  
  # yearly rain interactions -- hurdle. not working and idk how to run get_prior() properly for them
  # I think we want these to be negative since we're measuring the probability of zero eggs
  # prior(student_t(1, -0.25, 0.5), coef = syearly_rain_scaled:water_regimeperennial_1, dpar = hu),
  # prior(student_t(1, -0.5, 0.5), coef = syearly_rain_scaled:water_regimeseasonal_1, dpar = hu),
 
  # other covariates (feel free to change as you see fit)
  # prior(student_t(1, 0.5, 0.5), coef =  water_regimeseasonal), # slightly positive based on hypotheses
  prior(student_t(1, 0.25, 0.5), coef =  water_flowlentic), # slightly positive based on hypotheses
  prior(student_t(1, -0.25, 0.5), coef =  water_flowlotic), # slightly negtive based on hypotheses
  prior(student_t(1, -0.25, 0.5), coef =  sinterpolated_canopy_scaled_1), # slightly negative based on hypotheses
  prior(student_t(1, 0, 0.5), coef =  sBRDYEAR_scaled_1 ),
  prior(student_t(1, 0, 0.5), coef =  smean_percent_water_scaled_1),
  prior(student_t(1, 0, 0.5), coef =  sWaterTemp_scaled_1)
)

t0 <- Sys.time()

mod.hurdle.no.salinity <- brm(
  bf(num_egg_masses ~ 
       s(BRDYEAR_scaled) + 
       s(mean_percent_water_scaled) + 
       s(interpolated_canopy_scaled) +
       s(WaterTemp_scaled) +  
       #s(max_salinity_scaled, by = CoastalSite) + #Dropped
       #CoastalSite + #needed?
       s(yearly_rain_scaled, by =water_regime) +
       (water_flow) +
       (1 | Watershed/LocationID),
     hu ~ 
       s(yearly_rain_scaled, by = water_regime) +      # inflated model for zeros
       (1|Watershed/LocationID)),
  data = scaled_between_year,
  family = hurdle_negbinomial(),
  prior = bprior.no.sal,
  chains = 3, cores = 3,
  iter = 7500, # only need about 1000 for inference (3500-2500 warmup = 1000)
  warmup = 6500, 
  control = list(adapt_delta = 0.99)
)

t1<-Sys.time()
t1- t0

beepr::beep(0)

save(mod.hurdle.no.salinity, file = "Output/mod.hurdle.no.salinity.RData")
#load("Output/mod.hurdle.no.salinity.RData")

summary(mod.hurdle.no.salinity)
mod.hurdle <- mod.hurdle.no.salinity
summary(mod.hurdle)

#### hurdle model plots ####

#pairs(mod.brm)
pairs(mod.hurdle)
conditional_effects(mod.brm, surface = FALSE, prob = 0.8)
conditional_effects(mod.hurdle, surface = FALSE, prob = 0.8)
#from Mark
conditional_effects(mod.brm)|>
  plot(points = TRUE, theme = theme_classic())

#plot the hurdle effect by adding dpar
conditional_effects(mod.hurdle, dpar = "hu")

#another plot type that shows good interaction with rain and water_regime
conditional_smooths(mod.hurdle.no.salinity, prob = 0.9)


#try some plots with 
library(tidybayes)
#search tidybayes + brms for vignette

library(brms)
library(tidybayes)
library(ggplot2)
library(tidyverse)
# Tried something suggested by Chatgpt...

# Posterior predictive distribution plot
mod.brm|>
  tidybayes::add_predicted_draws(newdata = scaled_between_year, allow_new_levels = TRUE)|>
  ggplot(aes(x = .prediction, y = num_egg_masses)) +
  geom_point(alpha = 0.5) +
  geom_density_2d() +
  labs(x = "Predicted Number of Egg Masses", y = "Observed Number of Egg Masses") +
  theme_classic()

# Posterior intervals for fixed effects
mod.brm|>
  tidybayes::spread_draws(b_Intercept, b_CoastalSiteTRUE, b_zi_yearly_rain_scaled)|>
  tidyr::pivot_longer(cols = starts_with("b_"), names_to = "Variable", values_to = "Estimate") |>
  ggplot(aes(x = Estimate, y = Variable)) +
  tidybayes::geom_halfeyeh() +
  labs(x = "Estimate", y = "Variable") +
  theme_classic()

# Smooth terms visualization
mod.brm|>
  tidybayes::add_epred_draws(newdata = scaled_between_year, allow_new_levels = TRUE)|>
  ggplot(aes(x = BRDYEAR_scaled, y = .epred)) +
  tidybayes::stat_lineribbon(aes(y = num_egg_masses), .width = c(0.95, 0.8, 0.5)) +
  geom_point(data = scaled_between_year, aes(x = BRDYEAR_scaled, y = num_egg_masses), alpha = 0.5) +
  scale_fill_brewer() +
  labs(x = "BRDYEAR_scaled", y = "Number of Egg Masses") +
  theme_classic()

# Conditional effects plot
conditional_effects(mod.brm)|>
  plot(points = TRUE, theme = theme_classic())


# to get fitted values, newdata must match fitted data. 'Audubon Canyon', 'Easkoot Creek', 
# 'Garden Club Canyon', 'Olema Creek' were dropped from the model due to lack of data, so I am 
# dropping them here too
newdata <- scaled_between_year %>% 
  filter(Watershed != ('Audubon Canyon') & 
           Watershed != ('Easkoot Creek') & 
           Watershed != ('Garden Club Canyon') & 
           Watershed != ('Olema Creek')) %>% 
  droplevels()
  
fv <- add_epred_draws(mod.hurdle, newdata = newdata)

# yearly rain
rain_plot <- ggplot(fv, aes(x = yearly_rain)) + 
  stat_lineribbon(aes(y = .epred), .width = 0.89, size = 0.7, alpha = 0.5) +
  geom_point(aes(y = num_egg_masses), alpha = 0.007, size = 0.75, color = "darkblue") +
  labs(x = "Yearly Rain (inches)", y = "Predicted number of egg masses") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()
ggsave("rain_plot.jpg", width = 7, height = 6)
rain_plot

#Mark update: does this look better? use the geom_line() function with the stat_summary() function 
fv_summary <- fv %>%
  group_by(BRDYEAR) %>%
  summarize(mean_epred = mean(.epred), 
            lower = quantile(.epred, 0.05), 
            upper = quantile(.epred, 0.95))

rain_plot_modified2 <- ggplot(fv_summary, aes(x = BRDYEAR, y = mean_epred)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  labs(x = "Yearly Rain (inches)", y = "Predicted number of egg masses") +
  scale_fill_brewer(palette = "Greens") +
  theme_minimal()
rain_plot_modified2

#Mark: revised the codes a little.
#year
year_plot <- ggplot(fv, aes(x = BRDYEAR)) + 
  stat_lineribbon(aes(y = .epred), .width = 0.89, size = 0.7, alpha = 0.2, color = "blue", fill = "lightblue") +
  geom_point(aes(y = num_egg_masses), alpha = 0.1, size = 1, color = "darkblue") +
  labs(title = "Predicted Number of Egg Masses Over Years",
       x = "Year",
       y = "Predicted Number of Egg Masses",
       fill = "Confidence Level") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(year_plot)


# percent water (revised the codes a little for visualization)
percent_water_plot <- ggplot(fv, aes(x = mean_percent_water, y = .epred)) + 
  stat_lineribbon(aes(y = .epred), .width = 0.89, size = 0.7, alpha = 0.5) +
  geom_point(aes(y = num_egg_masses), alpha = 0.5, size = 1, color = "darkblue") +
  labs(x = "Percent Water Cover", y = "Predicted Number of Egg Masses") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
print(percent_water_plot)


# percent canopy
percent_canopy_plot <- ggplot(fv, aes(x = interpolated_canopy)) + 
  stat_lineribbon(aes(y = .epred), .width = 0.89, size = 0.7, alpha = 0.5) +
  geom_point(aes(y = num_egg_masses), alpha = 0.007, size = 0.75, color = "darkblue") +
  labs(x = "Percent canopy cover", y = "Predicted number of egg masses") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()
ggsave("canopy_cover_plot.jpg", width = 7, height = 6)
print(percent_canopy_plot)

# water temperature
water_temp_plot <- ggplot(fv, aes(x = WaterTemp)) + 
  stat_lineribbon(aes(y = .epred), .width = 0.89, size = 0.7, alpha = 0.5) +
  geom_point(aes(y = num_egg_masses), alpha = 0.007, size = 0.75, color = "darkblue") +
  labs(x = "Water Temperature (C)", y = "Predicted number of egg masses") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()
ggsave("water_temp_plot.jpg", width = 7, height = 6)
print(water_temp_plot)

# salinity for coastal sites
salinity_plot <- ggplot(fv, aes(x = max_salinity)) + 
  stat_lineribbon(aes(y = .epred), .width = 0.89, size = 0.7, alpha = 0.5) +
  geom_point(aes(y = num_egg_masses), alpha = 0.007, size = 0.75, color = "darkblue") +
  labs(x = "Salinity", y = "Predicted number of egg masses") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()
ggsave("salinity_plot.jpg", width = 7, height = 6)
print(salinity_plot)
