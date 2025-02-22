library(brms)
library(lme4)
library(marginaleffects)
library(cowplot)

scaled_between_year <- read_csv(here::here("data", "scaled_between_year.csv"))

# ZI linear model
##### priors: bprior.no.sal.linear.zi ####
bprior.no.sal.linear.zi <- c(
  #counts
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
  prior(normal(0.25, 1), coef = yearly_rain_lag_scaled:water_regimeseasonal)
  # prior(normal(0.25, 1), coef = proportion_high_water_vis),
  # prior(normal(0, 1), coef = proportion_na_water_vis)
  # prior(normal(0, 1), coef = water_visunknown)
)


##### model: mod.zi.no.salinity.linear ####
mod.zi.no.salinity.linear <- brm(
  num_egg_masses ~ 
       BRDYEAR_scaled + 
       mean_percent_water_scaled + 
       interpolated_canopy_scaled +
       WaterTemp_scaled +  
       mean_percent_sub_scaled +
       yearly_rain_scaled +
       yearly_rain_scaled : water_regime +
       yearly_rain_lag_scaled +
       water_regime +
       yearly_rain_lag_scaled : water_regime +
       water_flow +
       # proportion_high_water_vis +
       # proportion_na_water_vis +
       (1 | Watershed/LocationID),
  data = scaled_between_year, 
  family = zero_inflated_negbinomial(),
  prior = bprior.no.sal.linear.zi,
  chains = 3, cores = 3,
  iter = 13000,
  warmup = 12000,
  control = list(adapt_delta = 0.99))

summary(mod.zi.no.salinity.linear, prob = 0.89)


##### plots: mod.zi.no.salinity.linear #####
# using marginaleffects
plot_predictions(mod.zi.no.salinity.linear, by = "BRDYEAR_scaled", conf_level = 0.89)

pred <- predictions(mod.zi.no.salinity.linear, conf_level = 0.89, type = "prediction", ndraws = 10, re_formula = NA)
pred <- get_draws(pred)

write_csv(pred, here::here("data", "pred.csv"))

# unscaling response variables for plotting
col_means <- read_csv(here::here("data", "between_year_col_means.csv"))
col_sd <- read_csv(here::here("data", "between_year_col_sd.csv"))
pred_unscaled <- pred %>% 
  mutate(
    interpolated_canopy_unscaled = (interpolated_canopy_scaled * col_sd$interpolated_canopy) + col_means$interpolated_canopy,
    BRDYEAR_unscaled = (BRDYEAR_scaled * col_sd$BRDYEAR) + col_means$BRDYEAR,
    mean_percent_water_unscaled = (mean_percent_water_scaled * col_sd$mean_percent_water) + col_means$mean_percent_water
  )

# color palette because i want the plots to look pretty
# main_color <- "#CC5803"
# background <- "#FF9505"
# background2 <- "#FFB627"

# OR


main_color <- "#49741A"
background <- "#7DC82D"
background2 <- "#91D548"

# canopy -- significant
canopy_plot <- ggplot(pred_unscaled, aes(x = interpolated_canopy_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(aes(y = num_egg_masses), color = background, alpha = 0.035) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color, linewidth = 1.5) +
  labs(x = "Percent canopy cover", y = "Number of egg masses")

# percent open water -- significant
water_plot <- ggplot(pred_unscaled, aes(x = mean_percent_water_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(aes(y = num_egg_masses), color = background, alpha = 0.035) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color, linewidth = 1.5) +
  labs(x = "Percent open water cover", y = " ")

# BRDYEAR -- almost significant, nice to see trends over time
year_plot <- ggplot(pred_unscaled, aes(x = BRDYEAR_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(aes(y = num_egg_masses), color = background, alpha = 0.035) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color, linewidth = 1.5) +
  labs(x = "Water year", y = " ")

plot_grid(canopy_plot, water_plot, year_plot, nrow = 1)


# geom_smooth()# using sjPlot
conditional_effects(mod.zi.no.salinity.linear, surface = FALSE, prob = 0.89)

library(bayesplot)
posterior <- as.matrix(mod.zi.no.salinity.linear)

# scheme_green <- c(main_color, background, background2, "black", background, background2)

palette_green <- c(
  "#A5DD69",
  "#87D237",
  "#68A626",
  "#49741A",
  "#2A430F",
   "black")


color_scheme_set(palette_green)

# forest plot
mcmc_intervals(posterior, point_est = "mean", prob = 0.89, prob_outer = 0.89,
               inner_size = 1, 
               point_size = 3,
               point_args = list(shape = 21, fill = "white", color = "black", stroke = 1),
               pars = c(
                 "b_water_flowlentic",
                 "b_yearly_rain_scaled:water_regimeseasonal",
                 "b_mean_percent_sub_scaled",
                 "b_yearly_rain_lag_scaled",
                 "b_water_regimeseasonal",
                 "b_yearly_rain_scaled",
                 "b_WaterTemp_scaled",
                 "b_yearly_rain_lag_scaled:water_regimeseasonal",
                 "b_BRDYEAR_scaled",  
                 "b_mean_percent_water_scaled",
                 "b_interpolated_canopy_scaled",
                 "b_water_flowlotic"                                
                         
               )) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_y_discrete(labels = c(
    "b_BRDYEAR_scaled" = "Breeding year",
    "b_mean_percent_water_scaled" = "Mean percent water",
    "b_interpolated_canopy_scaled" = "Interpolated canopy cover",
    "b_WaterTemp_scaled" = "Water temperature",
    "b_mean_percent_sub_scaled" = "Mean percent submergent veg.",
    "b_yearly_rain_scaled" = "Yearly rain",
    "b_yearly_rain_lag_scaled" = "Lagged yearly rain",
    "b_water_regimeseasonal" = "Seasonal water regime",
    "b_water_flowlentic" = "Lentic flow",
    "b_water_flowlotic" = "Lotic flow",
    "b_yearly_rain_scaled:water_regimeseasonal" = "Yearly rain × Seasonal water",
    "b_yearly_rain_lag_scaled:water_regimeseasonal" = "Lagged rain × Seasonal water"
  )) +
  labs(x = "Estimate") +
  theme_minimal()

plot_model(mod.zi.no.salinity.linear, type = "pred", terms = c("yearly_rain_lag_scaled"))






##### OLD MODELS

# hurdle GAM model no salinity
##### priors: bprior.no.sal #######
bprior.no.sal <- c(
  #counts
  prior(normal(0, 0.5), coef = sBRDYEAR_scaled_1), 
  prior(normal(-0.5, 0.5), coef = sinterpolated_canopy_scaled_1), 
  prior(normal(0.25, 0.5), coef = smean_percent_sub_scaled_1), 
  prior(normal(0.5, 0.5), coef = smean_percent_water_scaled_1), 
  prior(normal(0.0, 1), coef = sWaterTemp_scaled_1), 
  prior(normal(0.0, 1), coef = syearly_rain_lag_scaled:water_regimeperennial_1), 
  prior(normal(0.0, 1), coef = syearly_rain_lag_scaled:water_regimeseasonal_1), 
  prior(normal(0.0, 1), coef = syearly_rain_scaled:water_regimeperennial_1), 
  prior(normal(0.0, 1), coef = syearly_rain_scaled:water_regimeseasonal_1), 
  prior(normal(0.0, 1), coef = water_flowlentic), 
  prior(normal(-0.5, 0.5), coef = water_flowlotic),
  #hu priors 
  prior(normal(0, 0.5), dpar = "hu", coef = sBRDYEAR_scaled_1), 
  prior(normal(-0.5, 0.5), dpar = "hu", coef = sinterpolated_canopy_scaled_1), 
  prior(normal(0.25, 0.5), dpar = "hu", coef = smean_percent_sub_scaled_1), 
  prior(normal(0.5, 0.5), dpar = "hu", coef = smean_percent_water_scaled_1), 
  prior(normal(0.0, 1), dpar = "hu", coef = sWaterTemp_scaled_1), 
  prior(normal(0.0, 1), dpar = "hu", coef = syearly_rain_lag_scaled:water_regimeperennial_1), 
  prior(normal(0.0, 1), dpar = "hu", coef = syearly_rain_lag_scaled:water_regimeseasonal_1), 
  prior(normal(0.0, 1), dpar = "hu", coef = syearly_rain_scaled:water_regimeperennial_1), 
  prior(normal(0.0, 1), dpar = "hu", coef = syearly_rain_scaled:water_regimeseasonal_1), 
  prior(normal(0.0, 1), dpar = "hu", coef = water_flowlentic), 
  prior(normal(-0.5, 0.5), dpar = "hu", coef = water_flowlotic)
)

##### model: mod.hurdle.no.salinity.gam ####  
t0 <- Sys.time()

# 2024-10-07 adding all terms to hurdle model 
mod.hurdle.no.salinity.gam <- brm(
  bf(num_egg_masses ~ 
       s(BRDYEAR_scaled) + 
       s(mean_percent_water_scaled) + 
       s(interpolated_canopy_scaled) +
       s(WaterTemp_scaled) +  
       s(mean_percent_sub_scaled) +
       s(yearly_rain_scaled) +
       s(yearly_rain_scaled, by = water_regime) +
       s(yearly_rain_lag_scaled) +
       s(yearly_rain_lag_scaled, by = water_regime) +
       water_flow +
       (1 | Watershed/LocationID),
     hu ~ 
       s(BRDYEAR_scaled) + 
       s(mean_percent_water_scaled) + 
       s(interpolated_canopy_scaled) +
       s(WaterTemp_scaled) +  
       s(mean_percent_sub_scaled) +
       s(yearly_rain_scaled) +
       s(yearly_rain_scaled, by = water_regime) +
       s(yearly_rain_lag_scaled) +
       s(yearly_rain_lag_scaled, by = water_regime) +
       water_flow +     
       (1|Watershed/LocationID)),
  data = scaled_between_year,  # run 4a to prep data file)
  family = hurdle_negbinomial(),
  prior = bprior.no.sal,
  chains = 3, cores = 3,
  iter = 3000, # 11500, # only need about 500 for inference
  warmup = 2800, #11000, 
  control = list(adapt_delta = 0.97)
)

t1<-Sys.time()
t1- t0

beepr::beep(0)

save(mod.hurdle.no.salinity.gam, file = "Output/mod.hurdle.no.salinity.gam.RData")
#load("Output/mod.hurdle.no.salinity.RData")

summary(mod.hurdle.no.salinity.gam)
mod.hurdle <- mod.hurdle.no.salinity.gam
summary(mod.hurdle)

##### plots: hurdle no salinity #####

pairs(mod.hurdle)
# conditional_effects(mod.brm, surface = FALSE, prob = 0.8)
conditional_effects(mod.hurdle, surface = FALSE, prob = 0.89)
#from Mark
conditional_effects(mod.hurdle)|>
  plot(points = TRUE, theme = theme_classic())

#plot the hurdle effect by adding dpar
conditional_effects(mod.hurdle, dpar = "hu")

#another plot type that shows good interaction with rain and water_regime
conditional_smooths(mod.hurdle, prob = 0.89)


## zi.gam

# ZI GAM model
##### priors: bprior.no.sal.zi.gam #####
bprior.no.sal.zi.gam <- c(
  #counts
  prior(normal(0, 0.5), coef = sBRDYEAR_scaled_1), 
  prior(normal(-0.5, 0.5), coef = sinterpolated_canopy_scaled_1), 
  prior(normal(0.25, 0.5), coef = smean_percent_sub_scaled_1), 
  prior(normal(0.5, 0.5), coef = smean_percent_water_scaled_1), 
  prior(normal(0.0, 1), coef = sWaterTemp_scaled_1), 
  prior(normal(0.0, 1), coef = syearly_rain_lag_scaled:water_regimeperennial_1), 
  prior(normal(0.0, 1), coef = syearly_rain_lag_scaled:water_regimeseasonal_1), 
  prior(normal(0.0, 1), coef = syearly_rain_scaled:water_regimeperennial_1), 
  prior(normal(0.0, 1), coef = syearly_rain_scaled:water_regimeseasonal_1), 
  prior(normal(0.0, 1), coef = water_flowlentic), 
  prior(normal(-0.5, 0.5), coef = water_flowlotic)
)


##### model: mod.zi.no.salinity.gam  ####
t0 <- Sys.time()

# 2024-10-07 adding all terms to hurdle model
mod.zi.no.salinity.gam <- brm(
  num_egg_masses ~ 
    s(BRDYEAR_scaled) + 
    s(mean_percent_water_scaled) + 
    s(interpolated_canopy_scaled) +
    s(WaterTemp_scaled) +  
    s(mean_percent_sub_scaled) +
    s(yearly_rain_scaled) +
    s(yearly_rain_scaled, by = water_regime) +
    s(yearly_rain_lag_scaled) +
    s(yearly_rain_lag_scaled, by = water_regime) +
    water_flow +
    (1 | Watershed/LocationID),
  data = scaled_between_year,  # run 4a to prep data file)
  family = zero_inflated_negbinomial(),
  prior = bprior.no.sal.zi.gam,
  chains = 3, cores = 3,
  iter = 3000, # 11500, # only need about 500 for inference
  warmup = 2800, #11000, 
  control = list(adapt_delta = 0.97)
)


summary(mod.zi.no.salinity.gam)


##### plots: mod.zi.no.salinity.gam #####

#pairs(mod.brm)
pairs(mod.zi.no.salinity.gam)
# conditional_effects(mod.brm, surface = FALSE, prob = 0.8)
conditional_effects(mod.zi.no.salinity.gam, surface = FALSE, prob = 0.89)
#from Mark
conditional_effects(mod.zi.no.salinity.gam)|>
  plot(points = TRUE, theme = theme_classic())

#another plot type that shows good interaction with rain and water_regime
conditional_smooths(mod.zi.no.salinity.gam, prob = 0.89)



# hurdle LINEAR model
##### lmer test #####
mod.1 <- glmer(
  num_egg_masses ~ 
    BRDYEAR_scaled + 
    mean_percent_water_scaled + 
    interpolated_canopy_scaled +
    WaterTemp_scaled +  
    mean_percent_sub_scaled +
    yearly_rain_scaled +
    yearly_rain_scaled : water_regime +
    yearly_rain_lag_scaled +
    water_regime +
    yearly_rain_lag_scaled : water_regime +
    water_flow +
    (water_flow + water_regime | Watershed/LocationID),
  data = scaled_between_year,  # run 4a to prep data file)
  family = negative.binomial(1)
)
summary(mod.1)
plot_model(mod.1, type = "diag", ci.lvl = 0.89)
plot_model(mod.1, type = "resid", ci.lvl = 0.89)
plot_model(mod.1)
plot_model(mod.1, type = "eff", ci.lvl = 0.89)
plot_model(mod.1, type = "int", ci.lvl = 0.89)

##### priors: bprior.no.sal.linear #####
get_prior(
  bf(num_egg_masses ~ 
       BRDYEAR_scaled + 
       mean_percent_water_scaled + 
       interpolated_canopy_scaled +
       WaterTemp_scaled +  
       mean_percent_sub_scaled +
       yearly_rain_scaled +
       yearly_rain_scaled : water_regime +
       yearly_rain_lag_scaled +
       water_regime +
       yearly_rain_lag_scaled : water_regime +
       water_flow +
       (water_flow + water_regime | Watershed/LocationID),
     hu ~ 
       yearly_rain_scaled * water_regime +      # inflated model for zeros
       (water_flow + water_regime | Watershed/LocationID)),
  data = scaled_between_year,  # run 4a to prep data file)
  family = hurdle_negbinomial())




bprior.no.sal.linear <- c(
  #counts
  prior(normal(0, 0.5), coef = BRDYEAR_scaled), 
  prior(normal(-0.5, 0.5), coef = interpolated_canopy_scaled), 
  prior(normal(0.25, 0.5), coef = mean_percent_sub_scaled), 
  prior(normal(0.5, 0.5), coef = mean_percent_water_scaled), 
  prior(normal(0.0, 1), coef = water_flowlentic), 
  prior(normal(0.0, 1), coef = water_flowlotic), 
  prior(normal(0.0, 1), coef = water_regimeseasonal), 
  prior(normal(0.0, 1), coef = yearly_rain_lag_scaled), 
  prior(normal(0.25, 1), coef = WaterTemp_scaled), 
  prior(normal(0.25, 1), coef = yearly_rain_scaled), 
  prior(normal(0.0, 1), coef = yearly_rain_scaled:water_regimeseasonal),
  prior(normal(0.0, 1), coef = yearly_rain_lag_scaled:water_regimeseasonal),
  #note hu parameters are p(0), so opposite normal intuition
  prior(normal(0, 0.5), dpar = "hu", coef = BRDYEAR_scaled), 
  prior(normal(0.25, 0.5), dpar = "hu", coef = interpolated_canopy_scaled), 
  prior(normal(0.25, 0.5), dpar = "hu", coef = mean_percent_sub_scaled), 
  prior(normal(-0.25, 0.5), dpar = "hu", coef = mean_percent_water_scaled), 
  prior(normal(0.0, 1), dpar = "hu", coef = water_flowlentic), 
  prior(normal(0.0, 1), dpar = "hu", coef = water_flowlotic), 
  prior(normal(0.0, 1), dpar = "hu", coef = water_regimeseasonal), 
  prior(normal(-0.0, 1), dpar = "hu", coef = yearly_rain_lag_scaled), 
  prior(normal(-0.25, 1), dpar = "hu", coef = WaterTemp_scaled), 
  prior(normal(-0.25, 1), dpar = "hu", coef = yearly_rain_scaled), 
  prior(normal(0.0, 1), dpar = "hu", coef = yearly_rain_scaled:water_regimeseasonal),
  prior(normal(0.0, 1), dpar = "hu", coef = yearly_rain_lag_scaled:water_regimeseasonal)
)




##### model: mod.hurdle.no.salinity.linear #####
t0 <- Sys.time()

mod.hurdle.no.salinity.linear <- brm(
  bf(num_egg_masses ~ 
       BRDYEAR_scaled + 
       mean_percent_water_scaled + 
       interpolated_canopy_scaled +
       WaterTemp_scaled +  
       mean_percent_sub_scaled +
       yearly_rain_scaled +
       yearly_rain_scaled : water_regime +
       yearly_rain_lag_scaled +
       water_regime +
       yearly_rain_lag_scaled : water_regime +
       water_flow +
       (1 | Watershed/LocationID),
     hu ~ 
       BRDYEAR_scaled + 
       mean_percent_water_scaled + 
       interpolated_canopy_scaled +
       WaterTemp_scaled +  
       mean_percent_sub_scaled +
       yearly_rain_scaled +
       yearly_rain_scaled : water_regime +
       yearly_rain_lag_scaled +
       water_regime +
       yearly_rain_lag_scaled : water_regime +
       water_flow +    # inflated model for zeros
       (1 | Watershed/LocationID)),
  data = scaled_between_year,  # run 4a to prep data file)
  family = hurdle_negbinomial(),
  prior = bprior.no.sal.linear,
  chains = 3, cores = 3,
  iter = 3000, # 11500, # only need about 500 for inference
  warmup = 2800, #11000, 
  control = list(adapt_delta = 0.97)
)

t1<-Sys.time()
t1- t0

beepr::beep(0)


mod.hurdle <- mod.hurdle.no.salinity.linear
summary(mod.hurdle, prob = 0.89)

##### plots: mod.hurdle.no.salinity.linear #####

#pairs(mod.brm)
# conditional_effects(mod.brm, surface = FALSE, prob = 0.8)
conditional_effects(mod.hurdle, surface = FALSE, prob = 0.89)
#from Mark
conditional_effects(mod.hurdle)|>
  plot(points = TRUE, theme = theme_classic())

#plot the hurdle effect by adding dpar
conditional_effects(mod.hurdle, dpar = "hu", prob = 0.89)

library(bayesplot)
posterior <- as.matrix(mod.hurdle)
mcmc_areas(posterior, 
           pars = c("b_BRDYEAR_scaled",  
                    "b_mean_percent_water_scaled",
                    "b_interpolated_canopy_scaled",
                    "b_WaterTemp_scaled",
                    "b_mean_percent_sub_scaled",
                    "b_yearly_rain_scaled",                             
                    "b_yearly_rain_lag_scaled",                              
                    "b_water_regimeseasonal",                                
                    "b_water_flowlentic",                                    
                    "b_water_flowlotic",
                    "b_yearly_rain_scaled:water_regimeseasonal",
                    "b_yearly_rain_lag_scaled:water_regimeseasonal",        
                    "b_hu_BRDYEAR_scaled",                               
                    "b_hu_mean_percent_water_scaled",                    
                    "b_hu_interpolated_canopy_scaled",                   
                    "b_hu_WaterTemp_scaled",                             
                    "b_hu_mean_percent_sub_scaled",                      
                    "b_hu_yearly_rain_scaled",                           
                    "b_hu_yearly_rain_lag_scaled",                       
                    "b_hu_water_regimeseasonal",                         
                    "b_hu_water_flowlentic",                             
                    "b_hu_water_flowlotic",                              
                    "b_hu_yearly_rain_scaled:water_regimeseasonal",      
                    "b_hu_yearly_rain_lag_scaled:water_regimeseasonal"   
           ), 
           prob = 0.5) +
  geom_vline(xintercept = 0, linetype = 2)


mcmc_intervals(posterior, point_est = "mean", prob = 0.89, prob_outer = 0.89,
               inner_size = 1,
               pars = c("b_BRDYEAR_scaled",  
                        "b_mean_percent_water_scaled",
                        "b_interpolated_canopy_scaled",
                        "b_WaterTemp_scaled",
                        "b_mean_percent_sub_scaled",
                        "b_yearly_rain_scaled",                             
                        "b_yearly_rain_lag_scaled",                              
                        "b_water_regimeseasonal",                                
                        "b_water_flowlentic",                                    
                        "b_water_flowlotic",
                        "b_yearly_rain_scaled:water_regimeseasonal",
                        "b_yearly_rain_lag_scaled:water_regimeseasonal",        
                        "b_hu_BRDYEAR_scaled",                               
                        "b_hu_mean_percent_water_scaled",                    
                        "b_hu_interpolated_canopy_scaled",                   
                        "b_hu_WaterTemp_scaled",                             
                        "b_hu_mean_percent_sub_scaled",                      
                        "b_hu_yearly_rain_scaled",                           
                        "b_hu_yearly_rain_lag_scaled",                       
                        "b_hu_water_regimeseasonal",                         
                        "b_hu_water_flowlentic",                             
                        "b_hu_water_flowlotic",                              
                        "b_hu_yearly_rain_scaled:water_regimeseasonal",      
                        "b_hu_yearly_rain_lag_scaled:water_regimeseasonal"   
               )) +
  geom_vline(xintercept = 0, linetype = 2)



