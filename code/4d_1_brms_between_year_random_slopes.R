# random slopes for year by county & location within watershed
# uses ALL data (i.e. since year is only covariate, so complete cases retains more data)
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

# data ####
# all site-year combinations, columns: num_egg_masses, BRDYEAR, Watershed, LocationID, County
scaled_btw_year_data_random_slopes <- read_csv(here::here("data", "scaled_btw_year_random_slopes.csv"))

# to help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# priors ####
bprior.zi.random.slopes.year <- c(
  prior(normal(0, 0.5), coef = BRDYEAR_scaled)
)

# model ####
mod.zi.random.slopes.year <- brm(
  num_egg_masses ~ 
    BRDYEAR_scaled +
    (BRDYEAR_scaled || Watershed/LocationID) +
    (BRDYEAR_scaled || County),
  data = scaled_btw_year_data_random_slopes, 
  family = zero_inflated_negbinomial(),
  prior = bprior.zi.random.slopes.year,
  chains = 3, 
  cores = 3,
  iter = 13000,
  warmup = 12000,
  sample_prior = TRUE,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod.zi.random.slopes.year, prob = 0.89)

# plots ####
ranef(mod.zi.random.slopes.year)
get_variables(mod.zi.random.slopes.year)

colors <- c(
  "Audubon Canyon"     = "#6ee7c3",  
  "Easkoot Creek"      = "#ff7f00",   
  "Garden Club Canyon" = "#c153c1",   
  "Kanoff Creek"       = "darkolivegreen",   
  "Laguna Salada"      = "#e7298a",   
  "Milagra Creek"      = "#cde15b",   
  "Oakwood Valley"     = "#d93102",  
  "Olema Creek"        = "#fb9a99",  
  "Redwood Creek"      = "#1f78b4",  
  "Rodeo Lagoon"       = "#cab2d6",  
  "San Mateo Creek"    = "#e6ab02",  
  "San Pedro Creek"    = "#501c87",  
  "Tennessee Valley"   = "#7eb42d",  
  "West Union"         = "#00747a",  
  "Wilkins Gulch"      = "#fdbf6f"  
)


## by county ####
rs_county <- as.matrix(mod.zi.random.slopes.year) %>%
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

ggplot(rs_county, aes(x = value, y = fct_rev(County))) +  # reverse to show top at top
  geom_density_ridges(
    alpha = 0.7,
    rel_min_height = 0.01,
    scale = 0.8,
    aes(fill=County)) +
  geom_point(aes(x = mean), color = "black", size = 1) +
  geom_linerange(aes(xmin = lower, xmax = upper), color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  scale_x_continuous(limits = c(-4, 4)) +  # adjust based on your slope scale
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

ggsave("year random slopes by county.png", path = here::here('Output'), width = 6, height = 4, units = "in")

### raw data vs trend for county ####
# plot faceted by watershed; x-axis = year, y-axis = eggs?
# based on sjPlot effects plots from 4d_brms_between_year.R
pred_RS <- predictions(mod.zi.random.slopes.year, conf_level = 0.89, type = "prediction", ndraws = 10, re_formula = NA)
pred_RS <- get_draws(pred_RS)

write_csv(pred_RS, here::here("data", "pred_RS.csv"))

# unscaling response variables for plotting
col_means_RS <- read_csv(here::here("data", "btw_year_random_slopes_col_means.csv"))
col_sd_RS <- read_csv(here::here("data", "btw_year_random_slopes_col_sd.csv"))

pred_unscaled_RS <- pred_RS %>% 
  mutate(
    BRDYEAR_unscaled = (BRDYEAR_scaled * col_sd$BRDYEAR) + col_means$BRDYEAR,
  )


## by location in watershed ####
rs_location <- as.matrix(mod.zi.random.slopes.year) %>%
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
rs_location <- rs_location %>%
  separate(Location, into = c("Watershed", "Site"), sep = "_", remove = FALSE) %>% 
  arrange(Watershed, Site) %>%
  mutate(Site = factor(Site, levels = unique(Site)))

ggplot(rs_location, aes(x = value, y = fct_rev(Site), fill=Watershed)) +  # reverse to show top at top
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
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.02))) +
  scale_fill_manual(values=colors)+
  labs(
    x = "Random Slope Deviation (Year)",
    y = "Site"
    # title = "Posterior Distributions of Random Slopes by Location",
    # subtitle = "89% Credible Intervals and Density Ridges"
  ) +
  theme_ridges(center_axis_labels = TRUE) +
  theme(legend.position = "right")

ggsave("year random slopes by site.png", path = here::here('Output'), width = 7, height = 10, units = "in")

## watershed ####
rs_watershed <- as.matrix(mod.zi.random.slopes.year) %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  filter(str_detect(parameter, "r_Watershed\\[.*?,BRDYEAR_scaled\\]")) %>% 
  mutate(parameter = str_remove_all(parameter, "r_Watershed\\[|,BRDYEAR_scaled\\]"),
         parameter = str_replace_all(parameter, "[.]", " ")) %>%
  rename(Watershed = parameter) %>%
  group_by(Watershed) %>%
  mutate(
    mean = mean(value),
    lower = quantile(value, 0.055),
    upper = quantile(value, 0.945)
  ) %>%
  ungroup() %>%
  arrange(desc(mean)) %>%
  mutate(Watershed = fct_inorder(Watershed, ordered = TRUE))

ggplot(rs_watershed, aes(x = value, y = fct_rev(Watershed))) +  # reverse to show top at top
  geom_density_ridges(
    alpha = 0.7,
    rel_min_height = 0.01,
    scale = 0.8,
    aes(fill=Watershed)) +
  geom_point(aes(x = mean), color = "black", size = 1) +
  geom_linerange(aes(xmin = lower, xmax = upper), color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = 2) +
  scale_x_continuous(limits = c(-4, 4)) +  # adjust based on your slope scale
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.06))) +
  scale_fill_manual(values=colors)+
  labs(
    x = "Random Slope Deviation (Year)",
    y = "Watershed"
  ) +
  theme_ridges(center_axis_labels = TRUE) +
  theme(legend.position = "right")

ggsave("year random slopes by watershed.png", path = here::here('Output'), width = 8.5, height = 6, units = "in")
