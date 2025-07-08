library(brms)
library(lme4)
library(marginaleffects)
library(cowplot)
library(bayesplot)
library(sjPlot)
library(cowplot)
library(priorsense)
library(tidyverse)
library(ggridges)
library(rstan)
library(tidybayes)
library(ggeffects)

scaled_between_year <- read_csv(here::here("data", "scaled_between_year.csv")) %>% 
  mutate(water_flow = as.factor(water_flow),
         water_regime = as.factor(water_regime),
         water_flow = fct_infreq(water_flow),
         LocationID = as.factor(LocationID),
         County = as.factor(County))

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
  prior(normal(0.0, 0.5), coef = mean_percent_water_scaled), # add squared term
  # prior(normal(0.5, 1), coef = water_flowlentic), 
  prior(normal(0, 1), coef = water_flowboth),
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
    # BRDYEAR_uncentered +
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
    # (1 | Watershed/LocationID) +
    (BRDYEAR_scaled || Watershed/LocationID) +
    (BRDYEAR_scaled || County),
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

ranef(mod.zi.no.salinity.linear)

get_variables(mod.zi.no.salinity.linear)
#### prior and posterior plots ####
prior_dist <- prior_draws(mod.zi.no.salinity.linear,
                          variable = c("b_BRDYEAR_scaled", 
                                       "b_mean_percent_water_scaled",
                                       "b_interpolated_canopy_scaled",
                                       "b_WaterTemp_scaled",
                                       "b_mean_percent_sub_scaled",
                                       "b_yearly_rain_scaled",
                                       "b_yearly_rain_lag_scaled",
                                       "b_water_regimeseasonal",
                                       # "b_water_flowboth",
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
    # "b_water_flowboth" ~ "Intermediate Water Flow",
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
                    # "b_water_flowboth",
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
    # "b_water_flowboth" ~ "Intermediate Water Flow",
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
  # filter(parameter == "Lentic Water Flow" | parameter == "Lotic Water Flow")

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


##### sjPlot effects plots: mod.zi.no.salinity.linear #####

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
    # BRDYEAR_unscaled = BRDYEAR_uncentered + 2003,
    mean_percent_water_unscaled = (mean_percent_water_scaled * col_sd$mean_percent_water) + col_means$mean_percent_water,
    lagged_rain_unscaled = (yearly_rain_lag_scaled * col_sd$yearly_rain_lag) + col_means$yearly_rain_lag,
    mean_percent_sub_unscaled = (mean_percent_sub_scaled * col_sd$mean_percent_sub) + col_means$mean_percent_sub,
    rain_unscaled = (yearly_rain_scaled * col_sd$yearly_rain) + col_means$yearly_rain,
    water_temp_unscaled = (WaterTemp_scaled * col_sd$WaterTemp) + col_means$WaterTemp
  )

# color palette because i want the plots to look pretty

main_color <- "#49741A"
background <- "#7DC82D"
# background2 <- "#91D548"

main_color_2 <- "#CC5803"
background_2 <- "#FF9505"
# background2 <- "#FFB627"

main_color_3 <- "#0070CC"
background_3 <- "#47ACFF"

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
         grey = {
           mc <- "black"
           bg <- "#54494B"
         },
         blue = {
           mc <- main_color_3
           bg <- background_3
         })
  plot_data <- as.data.frame(get_model_data(mod.zi.no.salinity.linear, type = "pred", terms = paste0(term, "_scaled [" , min_val, ":", max_val, ", by = 0.01]"), interval = "confidence", ci.lvl = 0.89)) %>% 
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
year_plot <- sjPlot_effects("BRDYEAR", "Breeding year", "grey", " ")
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

# geom_smooth()# using sjPlot
posterior <- as.matrix(mod.zi.no.salinity.linear)
palette_green <- c(
  "#A5DD69",
  "#87D237",
  "#68A626",
  "#49741A",
  "#2A430F",
   "black")

color_scheme_set(palette_green)

#### forest plot for appendix ####
mcmc_intervals(posterior, point_est = "mean", prob = 0.89, prob_outer = 0.89,
               inner_size = 1, 
               point_size = 2,
               pars = c(
                 # "b_water_flowlentic",
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
    # "b_BRDYEAR_uncentered" = "Breeding year",
    "b_mean_percent_water_scaled" = "Mean percent water",
    "b_interpolated_canopy_scaled" = "Interpolated canopy cover",
    "b_WaterTemp_scaled" = "Water temperature",
    "b_mean_percent_sub_scaled" = "Mean percent submergent veg.",
    "b_yearly_rain_scaled" = "Yearly rain",
    "b_yearly_rain_lag_scaled" = "Lagged yearly rain",
    "b_water_regimeseasonal" = "Seasonal water regime",
    # "b_water_flowlentic" = "Lentic flow",
    "b_water_flowlotic" = "Lotic flow",
    "b_yearly_rain_scaled:water_regimeseasonal" = "Yearly rain × Seasonal water",
    "b_yearly_rain_lag_scaled:water_regimeseasonal" = "Lagged rain × Seasonal water"
  )) +
  labs(x = "Estimate") +
  theme_minimal()

# plot_model(mod.zi.no.salinity.linear, type = "pred", terms = c("yearly_rain_lag_scaled"))



#### plotting random effects -- watersheds ####
# intercepts <- as.data.frame(ranef(mod.zi.no.salinity.linear)$Watershed) %>% 
#   mutate(Watershed = rownames(.))
# colnames(intercepts) <- str_remove_all(colnames(intercepts), ".Intercept")

re <- as.matrix(mod.zi.no.salinity.linear) %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  # filter(str_detect(parameter, "(r_Watershed)\\[")) %>%  # filter for Watershed
  filter(str_detect(parameter, "(r_Watershed)\\["), 
         !str_detect(parameter, ",BRDYEAR_scaled\\]")) %>% 
          # filter for Watershed INTERCEPTS (NOT Year random slopes)
  mutate(parameter = str_remove_all(parameter, "r_Watershed\\[|,Intercept\\]"),
         parameter = str_replace(parameter, "[.]", " ")) %>% 
  rename(Watershed = parameter) %>% 
  group_by(Watershed) %>% 
  mutate(mean = mean(value),
         lower = quantile(value, 0.055),
         upper = quantile(value, 0.945)) %>% # 89% CIs
  ungroup() %>% 
  arrange(desc(mean)) %>% 
  mutate(Watershed = fct_inorder(Watershed, ordered = T))

levels(re$Watershed)

ggplot(re, aes(x = value, y = Watershed)) +
  geom_density_ridges(alpha = 0.7, rel_min_height = 0.01, scale = 0.8) +
  geom_point(aes(x = mean)) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.06))) +
  labs(x = "Distribution") +
  theme_ridges(center_axis_labels = TRUE)

#### plotting random effects -- sites within watersheds ####
# intercepts <- as.data.frame(ranef(mod.zi.no.salinity.linear)$Watershed) %>% 
#   mutate(Watershed = rownames(.))
# colnames(intercepts) <- str_remove_all(colnames(intercepts), ".Intercept")

re <- as.matrix(mod.zi.no.salinity.linear) %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>% 
  filter(str_detect(parameter, "r_Watershed\\:"), 
         !str_detect(parameter, ",BRDYEAR_scaled\\]")) %>%
        # filter for Watershed:LocationID random INTERCEPTS, exclude random slopes for BRDYEAR
  mutate(parameter = str_remove_all(parameter, "r_Watershed:LocationID\\[|,Intercept\\]"),
         parameter = str_replace(parameter, "[.]", " ")) %>% 
  # rename(Watershed = parameter) %>% 
  separate_wider_delim(parameter,  delim = "_", names = c("Watershed", "Site")) %>% 
  group_by(Watershed, Site) %>% 
  mutate(mean = mean(value),
         lower = quantile(value, 0.055),
         upper = quantile(value, 0.945)) %>% # 89% CIs
  ungroup() %>% 
  arrange(desc(Watershed), desc(mean)) %>% 
  mutate(Site = fct_inorder(Site, ordered = T),
         Watershed = fct(Watershed))
  

ggplot(re, aes(x = value, y = Site)) +
  geom_density_ridges(alpha = 0.7, rel_min_height = 0.01, scale = 0.8, aes(fill = Watershed)) +
  geom_point(aes(x = mean)) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  geom_vline(aes(xintercept = 0), color = "black", linetype = 2) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.06))) +
  labs(x = "Random Intercept (log scale)") +
  theme_ridges(center_axis_labels = TRUE)


#### (added July 3) plotting random slopes -- BRDYEAR for LocationID ####
re_locID <- as.matrix(mod.zi.no.salinity.linear) %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  filter(str_detect(parameter, "r_Watershed:LocationID\\[.*?,BRDYEAR_scaled\\]")) %>% 
  #filter only random slopes for Watershed:LocationID
  mutate(parameter = str_remove_all(parameter, "r_Watershed:LocationID\\[|,BRDYEAR_scaled\\]"),
         parameter = str_replace_all(parameter, "[.]", " ")) %>%
  rename(Location = parameter) %>%
  group_by(Location) %>%
  mutate(
    mean = mean(value),
    lower = quantile(value, 0.055),
    upper = quantile(value, 0.945)
  ) %>%
  ungroup() %>%
  arrange(desc(mean)) %>%
  mutate(Location = fct_inorder(Location, ordered = TRUE))

## following the plots Robin made
re_locID <- re_locID %>%
  separate(Location, into = c("Watershed", "Site"), sep = "_", remove = FALSE)

ggplot(re_locID, aes(x = value, y = fct_rev(Site), fill=Watershed)) +  # reverse to show top at top
  geom_density_ridges(
    alpha = 0.7,
    rel_min_height = 0.01,
    scale = 0.8,
    aes(fill = Watershed)
  ) +
  geom_point(aes(x = mean), color = "black", size = 1) +
  geom_linerange(aes(xmin = lower, xmax = upper), color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  scale_x_continuous(limits = c(-3, 3)) +  # adjust based on your slope scale
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.06))) +
  labs(
    x = "Random Slope Deviation (Year)",
    y = "Site"
    # title = "Posterior Distributions of Random Slopes by Location",
    # subtitle = "89% Credible Intervals and Density Ridges"
  ) +
  theme_ridges(center_axis_labels = TRUE) +
  theme(legend.position = "right")

#### NEW JULY 3: plotting random slopes by County for BRDYEAR ####
re_county <- as.matrix(mod.zi.no.salinity.linear) %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  filter(str_detect(parameter, "r_County\\[.*?,BRDYEAR_scaled\\]")) %>% 
  mutate(parameter = str_remove_all(parameter, "r_County\\[|,BRDYEAR_scaled\\]"),
         parameter = str_replace_all(parameter, "[.]", " ")) %>%
  rename(County = parameter) %>%
  group_by(County) %>%
  mutate(
    mean = mean(value),
    lower = quantile(value, 0.055),
    upper = quantile(value, 0.945)
  ) %>%
  ungroup() %>%
  arrange(desc(mean)) %>%
  mutate(County = fct_inorder(County, ordered = TRUE))

# palette
marin_color = c("yellow4")
sanmateo_color = c('orchid4')

ggplot(re_county, aes(x = value, y = fct_rev(County))) +  # reverse to show top at top
  geom_density_ridges(
    alpha = 0.7,
    rel_min_height = 0.01,
    scale = 0.8,
    aes(fill=County)) +
  geom_point(aes(x = mean), color = "black", size = 1) +
  geom_linerange(aes(xmin = lower, xmax = upper), color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  scale_x_continuous(limits = c(-3, 3)) +  # adjust based on your slope scale
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.06))) +
  scale_fill_manual(values = c("Marin" = marin_color, "San Mateo" = sanmateo_color))+
  labs(
    x = "Random Slope Deviation (Year)",
    y = "County",
    # title = "Posterior Distributions of Random Slopes by County",
    # subtitle = "89% Credible Intervals and Density Ridges"
  ) +
  theme_ridges(center_axis_labels = TRUE) +
  theme(legend.position = "right")


#### marginaleffects by hand plots -- old ####
# canopy -- significant
canopy_plot <- ggplot(pred_unscaled, aes(x = interpolated_canopy_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(data = scaled_between_year, aes(x = interpolated_canopy, y = num_egg_masses), color = background, alpha = 0.3) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color, linewidth = 1.5) +
  labs(x = "Percent canopy cover", y = "Number of egg masses")

# percent open water -- significant
water_plot <- ggplot(pred_unscaled, aes(x = mean_percent_water_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(data = scaled_between_year, aes(x = mean_percent_water, y = num_egg_masses), color = background, alpha = 0.3) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color, linewidth = 1.5) +
  labs(x = "Percent open water cover", y = " ")

# BRDYEAR -- almost significant, nice to see trends over time
year_plot <- ggplot(pred_unscaled, aes(x = BRDYEAR_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(data = scaled_between_year, aes(x = BRDYEAR, y = num_egg_masses), color = background_2, alpha = 0.3) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color_2, linewidth = 1.5) +
  labs(x = "Water year", y = " ")

# lagged yearly rain -- almost significant
lag_rain_plot <- ggplot(pred_unscaled, aes(x = lagged_rain_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(data = scaled_between_year, aes(x = yearly_rain_lag, y = num_egg_masses), color = background_3, alpha = 0.3) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color_3, linewidth = 1.5) +
  labs(x = "Lagged yearly rainfall (cm)", y = " ")

# percent submergent vegetation -- almost significant
sub_veg_plot <- ggplot(pred_unscaled, aes(x = mean_percent_sub_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(data = scaled_between_year, aes(x = mean_percent_sub, y = num_egg_masses), color = background_3, alpha = 0.3) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", color = main_color_3, linewidth = 1.5) +
  labs(x = "Percent submergent vegetation", y = "Number of egg masses")

# water temp
water_temp_plot <- ggplot(pred_unscaled, aes(x = water_temp_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(data = scaled_between_year, aes(x = WaterTemp, y = num_egg_masses), color = background_3, alpha = 0.3) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(stat = "smooth", linewidth = 1.5, color = main_color_3) +
  labs(x = "Water temperature (°C)", y = " ")

cowplot::plot_grid(canopy_plot, water_plot, year_plot, sub_veg_plot, lag_rain_plot, water_temp_plot, nrow = 2, align = "hv")

# water temp
rain_water_regime_plot <- ggplot(pred_unscaled, aes(x = rain_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(-1, 175)) +
  theme_bw() +
  geom_point(aes(y = num_egg_masses), color = background, alpha = 0.035) +
  geom_line(aes(y = conf.low), color = 'black', stat = "smooth", alpha = 0.5) +
  geom_line(aes(y = conf.high), color = 'black', stat = "smooth", alpha = 0.5) +
  geom_line(stat = "smooth", linewidth = 1.5, color=  main_color) +
  labs(x = "Yearly rainfall (cm)", y = "Number of egg masses") +
  facet_wrap(~water_regime, labeller = labeller(water_regime = c(perennial = "Perennial Water Regime", seasonal = "Seasonal Water Regime")))
rain_water_regime_plot
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



