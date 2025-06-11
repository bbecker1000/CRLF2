library(brms)
# library(lme4)
library(marginaleffects)
library(cowplot)
library(bayesplot)
library(sjPlot)
library(priorsense)
# library(tidyverse)
library(ggridges)
library(rstan)
library(tidybayes)

scaled_between_year <- read_csv(here::here("data", "scaled_between_year.csv")) %>% 
  mutate(water_flow = as.factor(water_flow),
         water_regime = as.factor(water_regime),
         water_flow = fct_infreq(water_flow))

set.seed(42) # so the model will give us the same results each time

# to help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ZI linear model
##### priors: bprior.no.sal.linear.zi ####
bprior.no.sal.linear.zi <- c(
  #counts
  prior(normal(0, 0.5), coef = BRDYEAR_scaled),
  prior(normal(-0.5, 0.5), coef = interpolated_canopy_scaled), 
  prior(normal(0.25, 0.5), coef = mean_percent_sub_scaled), 
  prior(normal(0.0, 0.5), coef = mean_percent_water_scaled),
  prior(normal(0, 1), coef = water_flowboth),
  prior(normal(-0.5, 1), coef = water_flowlotic), 
  prior(normal(0.5, 1), coef = water_regimeseasonal), 
  prior(normal(0.5, 1), coef = yearly_rain_lag_scaled), 
  prior(normal(0.25, 1), coef = WaterTemp_scaled),
  prior(normal(0.25, 1), coef = yearly_rain_scaled), 
  prior(normal(0.25, 1), coef = yearly_rain_scaled:water_regimeseasonal),
  prior(normal(0.25, 1), coef = yearly_rain_lag_scaled:water_regimeseasonal)
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
    (1 | Watershed/LocationID),
  data = scaled_between_year, 
  family = zero_inflated_negbinomial(),
  prior = bprior.no.sal.linear.zi,
  chains = 3, cores = 3,
  iter = 13000,
  warmup = 10000,
  sample_prior = TRUE,
  control = list(adapt_delta = 0.99))

summary(mod.zi.no.salinity.linear, prob = 0.89)
powerscale_sensitivity(mod.zi.no.salinity.linear)
powerscale_plot_dens(mod.zi.no.salinity.linear)

# plotting prior and posterior distributions

prior_dist <- prior_draws(mod.zi.no.salinity.linear,
                          variable = c("b_BRDYEAR_scaled", 
                                       "b_mean_percent_water_scaled",
                                       "b_interpolated_canopy_scaled",
                                       "b_WaterTemp_scaled",
                                       "b_mean_percent_sub_scaled",
                                       "b_yearly_rain_scaled",
                                       "b_yearly_rain_lag_scaled",
                                       "b_water_regimeseasonal",
                                       "b_water_flowlotic",
                                       "b_yearly_rain_scaled:water_regimeseasonal",
                                       "b_yearly_rain_lag_scaled:water_regimeseasonal"
                          )) %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  mutate(parameter = case_match(
    parameter,
    "b_BRDYEAR_scaled" ~ "Breeding Year",
    "b_mean_percent_water_scaled" ~ "Mean Percent Water",
    "b_interpolated_canopy_scaled" ~ "Mean Percent Canopy",
    "b_WaterTemp_scaled" ~ "Water Temperature",
    "b_mean_percent_sub_scaled" ~ "Percent Submergent Vegetation",
    "b_yearly_rain_scaled" ~ "Yearly Rain",
    "b_yearly_rain_lag_scaled" ~ "Lagged Yearly Rain",
    "b_water_regimeseasonal" ~ "Seasonal Water Regime",
    "b_water_flowlotic" ~ "Lotic Water Flow",
    "b_yearly_rain_scaled:water_regimeseasonal" ~ "Yearly Rain Scaled x Seasonal Water Regime",
    "b_yearly_rain_lag_scaled:water_regimeseasonal" ~ "Lagged Yearly Rain Scaled x Seasonal Water Regime"
  ),
  distribution = "Prior")


posterior_dist <- as.matrix(mod.zi.no.salinity.linear) %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  filter(parameter %in% c("b_BRDYEAR_scaled", 
                          "b_mean_percent_water_scaled",
                          "b_interpolated_canopy_scaled",
                          "b_WaterTemp_scaled",
                          "b_mean_percent_sub_scaled",
                          "b_yearly_rain_scaled",
                          "b_yearly_rain_lag_scaled",
                          "b_water_regimeseasonal",
                          "b_water_flowlotic",
                          "b_yearly_rain_scaled:water_regimeseasonal",
                          "b_yearly_rain_lag_scaled:water_regimeseasonal"
  )) %>% 
  mutate(parameter = case_match(
    parameter,
    "b_BRDYEAR_scaled" ~ "Breeding Year",
    "b_mean_percent_water_scaled" ~ "Mean Percent Water",
    "b_interpolated_canopy_scaled" ~ "Mean Percent Canopy",
    "b_WaterTemp_scaled" ~ "Water Temperature",
    "b_mean_percent_sub_scaled" ~ "Percent Submergent Vegetation",
    "b_yearly_rain_scaled" ~ "Yearly Rain",
    "b_yearly_rain_lag_scaled" ~ "Lagged Yearly Rain",
    "b_water_regimeseasonal" ~ "Seasonal Water Regime",
    "b_water_flowlotic" ~ "Lotic Water Flow",
    "b_yearly_rain_scaled:water_regimeseasonal" ~ "Yearly Rain Scaled x Seasonal Water Regime",
    "b_yearly_rain_lag_scaled:water_regimeseasonal" ~ "Lagged Yearly Rain Scaled x Seasonal Water Regime"
  ),
  distribution = "Posterior")

prior_post <- rbind(prior_dist, posterior_dist) %>% 
  group_by(parameter, distribution) %>% 
  mutate(mean = mean(value)) %>% 
  ungroup() %>% 
  mutate(parameter = fct_reorder(parameter, desc(mean)),
         distribution = fct_inorder(distribution))

# color palette
prior_color <- "#54494B"
posterior_color <- "#7DC82D"


prior_post_plot <- ggplot(data = prior_post, aes(x = value, y = parameter, fill = distribution)) +
  geom_density_ridges(alpha = 0.5, rel_min_height = 0.01, scale = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.06))) +
  theme_ridges(center_axis_labels = TRUE) +
  scale_fill_manual(values = c("Prior" = prior_color, "Posterior" = posterior_color))
prior_post_plot

##### plots: mod.zi.no.salinity.linear #####
# using marginaleffects

pred <- predictions(mod.zi.no.salinity.linear, conf_level = 0.89, type = "prediction", ndraws = 10, re_formula = NA)
pred <- get_draws(pred)

# unscaling response variables for plotting
col_means <- read_csv(here::here("data", "between_year_col_means.csv"))
col_sd <- read_csv(here::here("data", "between_year_col_sd.csv"))
pred_unscaled <- pred %>% 
  mutate(
    interpolated_canopy_unscaled = (interpolated_canopy_scaled * col_sd$interpolated_canopy) + col_means$interpolated_canopy,
    BRDYEAR_unscaled = (BRDYEAR_scaled * col_sd$BRDYEAR) + col_means$BRDYEAR,
    mean_percent_water_unscaled = (mean_percent_water_scaled * col_sd$mean_percent_water) + col_means$mean_percent_water,
    lagged_rain_unscaled = (yearly_rain_lag_scaled * col_sd$yearly_rain_lag) + col_means$yearly_rain_lag,
    mean_percent_sub_unscaled = (mean_percent_sub_scaled * col_sd$mean_percent_sub) + col_means$mean_percent_sub,
    rain_unscaled = (yearly_rain_scaled * col_sd$yearly_rain) + col_means$yearly_rain,
    water_temp_unscaled = (WaterTemp_scaled * col_sd$WaterTemp) + col_means$WaterTemp,
  )

# color palettes
main_color <- "#49741A"
background <- "#7DC82D"

main_color_2 <- "#CC5803"
background_2 <- "#FF9505"

main_color_3 <- "#0070CC"
background_3 <- "#47ACFF"


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
  geom_point(aes(y = num_egg_masses), color = background_2, alpha = 0.035) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color_2, linewidth = 1.5) +
  labs(x = "Water year", y = " ")

# lagged yearly rain -- almost significant
lag_rain_plot <- ggplot(pred_unscaled, aes(x = lagged_rain_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(aes(y = num_egg_masses), color = background_3, alpha = 0.035) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color_3, linewidth = 1.5) +
  labs(x = "Lagged yearly rainfall (cm)", y = " ")

# percent submergent vegetation -- almost significant
sub_veg_plot <- ggplot(pred_unscaled, aes(x = mean_percent_sub_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(aes(y = num_egg_masses), color = background_3, alpha = 0.035) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color_3, linewidth = 1.5) +
  labs(x = "Percent submergent vegetation", y = "Number of egg masses")

# water temp
water_temp_plot <- ggplot(pred_unscaled, aes(x = water_temp_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(aes(y = num_egg_masses), color = background_3, alpha = 0.035) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", linewidth = 1.5, color = main_color_3) +
  labs(x = "Water temperature (°C)", y = " ")

cowplot::plot_grid(canopy_plot, water_plot, year_plot, sub_veg_plot, lag_rain_plot, water_temp_plot, nrow = 2, align = "hv")

palette_green <- c(
  "#A5DD69",
  "#87D237",
  "#68A626",
  "#49741A",
  "#2A430F",
  "black")

color_scheme_set(palette_green)

# forest plot in appendix
posterior <- as.matrix(mod.zi.no.salinity.linear)
mcmc_intervals(posterior, point_est = "mean", prob = 0.89, prob_outer = 0.89,
               inner_size = 1, 
               point_size = 2,
               pars = c(
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
    "b_water_flowlotic" = "Lotic flow",
    "b_yearly_rain_scaled:water_regimeseasonal" = "Yearly rain × Seasonal water",
    "b_yearly_rain_lag_scaled:water_regimeseasonal" = "Lagged rain × Seasonal water"
  )) +
  labs(x = "Estimate") +
  theme_minimal()
