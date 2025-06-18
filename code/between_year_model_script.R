library(brms)
library(marginaleffects)
library(cowplot)
library(bayesplot)
library(sjPlot)
library(priorsense)
library(ggridges)
library(rstan)
library(tidybayes)
library(ggeffects)

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
# color palettes
main_color <- "#49741A"
background <- "#7DC82D"

main_color_2 <- "#CC5803"
background_2 <- "#FF9505"

main_color_3 <- "#0070CC"
background_3 <- "#47ACFF"

# effects plots using sjPlot
scaled_between_year_round <- scaled_between_year %>% 
  mutate(across(where(is.double), ~ round(., digits = 2)))

sjPlot_effects <- function(term, xlab, color, ylab = "Number of egg masses") {
  min_val <- min(scaled_between_year_round[[paste0(term, "_scaled")]])
  max_val <- max(scaled_between_year_round[[paste0(term, "_scaled")]])
  
  switch(color,
         green = {
           mc <- main_color
           bg <- background
         },
         orange = {
           mc <- main_color_2
           bg <- background_2
         },
         blue = {
           mc <- main_color_3
           bg <- background_3
         })
  plot_data <- as.data.frame(get_model_data(mod.zi.no.salinity.linear, type = "pred", terms = paste0(term, "_scaled [" , min_val, ":", max_val, ", by = 0.01]"), interval = "confidence")) %>% 
    select(-group, -group_col) %>% 
    mutate(unscaled = (x * col_sd[[term]]) + col_means[[term]])
  
  xlim <- c(round(min(plot_data$unscaled)), round(max(plot_data$unscaled)))
  
  ggplot(plot_data, aes(x = unscaled, y = predicted)) +
    # point data for appendix plot
    # geom_point(data = scaled_between_year, aes(x = .data[[term]], y = num_egg_masses), alpha = 0.5, color = bg) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, fill = bg) +
    geom_line(linewidth = 1, color = mc) +
    theme_bw() +
    labs(x = xlab, y = ylab) +
    scale_x_continuous(limits = xlim) +
  scale_y_continuous(limits = c(-1, 55))
}
canopy_plot <- sjPlot_effects("interpolated_canopy", "Percent canopy cover", "green")
water_plot <- sjPlot_effects("mean_percent_water", "Percent open water", "green", " ")
year_plot <- sjPlot_effects("BRDYEAR", "Breeding year", "orange", " ")
sub_veg_plot <- sjPlot_effects("mean_percent_sub", "Percent submergent vegetation", "blue")
lag_rain_plot <- sjPlot_effects("yearly_rain_lag", "Lagged yearly rain (cm)", "blue", " ")
water_temp_plot <- sjPlot_effects("WaterTemp", "Water temperature (°C)", "blue", " ")
yearly_rain_plot <- sjPlot_effects("yearly_rain", "Yearly rain (cm)", "green", " ")

# interaction plot for water regime and yearly rain
# directly calling ggpredict because sjPlot's wrapper (get_model_data) won't let me get interaction data at intervals of 0.01
plot_data <- ggpredict(mod.zi.no.salinity.linear, 
                       terms = c("yearly_rain_scaled [-1.58:1.58 by=0.01]", "water_regime"),
                       ci.level = 0.89,
                       interval = "confidence") %>% 
  as.data.frame() %>% 
  mutate(unscaled = (x * col_sd$yearly_rain) + col_means$yearly_rain)

interaction_plot <- ggplot(plot_data, aes(x = unscaled, y = predicted)) +
  # geom_point(data = scaled_between_year, aes(x = yearly_rain, y = num_egg_masses, color = water_regime), alpha = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(linewidth = 1, aes(color = group)) +
  theme_bw() +
  labs(x = "Yearly rain (cm)", y = "Number of egg masses", color = "Water regime", fill = "Water regime") +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = c("black", "#49741A")) +
  scale_fill_manual(values = c("#54494B", "#7DC82D"))

# 3 cols, 2 rows
effects_plots <- cowplot::plot_grid(canopy_plot, water_plot, year_plot, sub_veg_plot, lag_rain_plot, water_temp_plot, 
                                    nrow = 2, align = "hv")

cowplot::plot_grid(effects_plots, interaction_plot,
                   nrow = 2, rel_heights = c(2, 1))




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
