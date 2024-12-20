library(brms)
library(tidyverse)
library(marginaleffects)

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

## unscaling
# for unscaling later
col_means_lag <- lag_between_year %>% 
  filter(complete_case == TRUE) %>%
  select(-Watershed, -LocationID, -dry_year, -CoastalSite, -water_flow, -water_regime, -complete_case) %>% 
  colMeans() %>% 
  t() %>% 
  as.data.frame()

write_csv(col_means_lag, here::here("data", "between_year_col_means_lag.csv"))

col_sd_lag <- lag_between_year %>% 
  filter(complete_case == TRUE) %>%
  select(-Watershed, -LocationID, -dry_year, -CoastalSite, -water_flow, -water_regime, -complete_case) %>% 
  apply(., 2, sd) %>% 
  t() %>% 
  as.data.frame()

write_csv(col_sd_lag, here::here("data", "between_year_col_sd_lag.csv"))






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

summary(lag.zi.linear, prob = 0.89) # year and canopy sig! Nov 20, 2024

#### plots ####
conditional_effects(lag.zi.linear, surface = FALSE, prob = 0.89)

# Try using the package sjPlot and function plot_model(MODELNAME, type = “eff”, terms = c(“TERM1”, “TERM2…)
# Also try type = “int”
# library(sjPlot)
# plot_model(lag.zi.linear, type = "pred", term="BRDYEAR_scaled")
# plot_model(lag.zi.linear, type = "pred", term="interpolated_canopy_scaled")
# plot_model(lag.zi.linear, type = "pred", term="yearly_rain_lag_scaled")
# plot_model(lag.zi.linear, type = "pred", term="water_regime")
# 
# plot_model(lag.zi.linear, type = "int", mdrt.values = "meansd")    


##### marginal effects #####
# using marginaleffects
plot_predictions(lag.zi.linear, by = "BRDYEAR_scaled", conf_level = 0.89)

pred_lag <- predictions(lag.zi.linear, conf_level = 0.89, type = "prediction", ndraws = 10, re_formula = NA)
pred_lag <- get_draws(pred_lag)

write_csv(pred_lag, here::here("data", "pred_lag.csv"))

# trying to unscale response variables for plotting
col_means_lag <- read_csv(here::here("data", "between_year_col_means_lag.csv"))
col_sd_lag <- read_csv(here::here("data", "between_year_col_sd_lag.csv"))

## TODO: DL to go through and put the relevant columns and names
# sig variables = year and canopy
pred_unscaled_lag <- pred_lag %>% 
  mutate(
    interpolated_canopy_unscaled = (interpolated_canopy_scaled * col_sd_lag$interpolated_canopy) + col_means_lag$interpolated_canopy,
    BRDYEAR_unscaled = (BRDYEAR_scaled * col_sd_lag$BRDYEAR) + col_means_lag$BRDYEAR
  )

# color palette because i want the plots to look pretty
main_color <- "#CC5803"
background <- "#FF9505"
background2 <- "#FFB627"

# OR
main_color <- "#0B5563"
background <- "#5299D3"
background2 <- "#BEB8EB"


# canopy -- significant
canopy_plot_lag <- ggplot(pred_unscaled_lag, aes(x = interpolated_canopy_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(0, 200)) +
  theme_bw() +
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_point(aes(y = num_egg_masses), color = background, alpha = 0.035) +
  geom_line(stat = "smooth", color = main_color, linewidth = 1.5) +
  labs(x = "Percent canopy cover", y = "Number of egg masses")


plot_grid(canopy_plot, water_plot, nrow = 1)

# BRDYEAR -- almost significant, nice to see trends over time
year_plot_lag <- ggplot(pred_unscaled_lag, aes(x = BRDYEAR_unscaled, y = estimate)) +
  scale_y_continuous(limits = c(0, 200)) +
  theme_bw() +
  geom_line(aes(y = conf.low), stat = "smooth", color = "black", alpha = 0.5) +
  geom_line(aes(y = conf.high), stat = "smooth", color = "black", alpha = 0.5) +
  geom_point(aes(y = num_egg_masses), color = background, alpha = 0.035) +
  geom_line(stat = "smooth", color = main_color, linewidth = 1.5) +
  labs(x = "Water year", y = "Number of egg masses")
            