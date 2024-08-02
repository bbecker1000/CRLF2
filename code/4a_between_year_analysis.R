library(glmmTMB)
library(tidyverse)
library(here)
#library(bbmle) 
library(lme4)
library(sjPlot)
library(mgcv)
library(gamlss)
library(gratia)
library(gam.hp)
library(patchwork)

between_year_data <- read_csv(here::here("data", "between_year_data.csv"))  %>% 
  mutate(LocationID = as.factor(LocationID),
         Watershed = as.factor(Watershed), 
         LocationInWatershed = interaction(Watershed, LocationID))

#### scaling covariates ####
# creating a "complete case" column
between_year_data$complete_case <- complete.cases(between_year_data)
complete_btw_data <- between_year_data %>% filter(complete_case == TRUE)

# write to CSV
write_csv(complete_btw_data, here::here("data", "complete_btw_data.csv"))

# scaling covariates
scaled_between_year <- complete_btw_data %>% 
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
    WaterTemp_scaled = as.vector(scale(WaterTemp))
    # mean_salinity_scaled = as.vector(scale(mean_salinity)),
    # max_salinity_scaled = as.vector(scale(max_salinity)),
  )

### zero-inflated GAM model ####
# gamlss uses pb() instead of s()

# testing correlations, removed mean_percent_emerg because it was highly correlated with percent water
cor(scaled_between_year[, c("BRDYEAR", "mean_percent_emerg", "mean_percent_sub", 
                            "mean_percent_water", "interpolated_canopy", "yearly_rain", 
                            "max_depth", "WaterTemp", "max_salinity")])


between_year_gamlss <- gamlss(formula = 
                                num_egg_masses ~ pb(BRDYEAR_scaled) + 
                                pb(mean_percent_sub_scaled) +
                                pb(mean_percent_water_scaled) +
                                pb(interpolated_canopy_scaled) +
                                max_depth_scaled +
                                pb(WaterTemp_scaled) +
                                max_salinity_scaled:as.factor(CoastalSite) +
                                pb(yearly_rain_scaled)+
                                water_regime + # added new covariate & interaction
                                yearly_rain_scaled:as.factor(water_regime)+
                                re(random = ~1 | Watershed/LocationID),
                              nu.formula = ~ 
                                max_depth_scaled +
                                pb(yearly_rain_scaled) +
                                re(random = ~1 | Watershed/LocationID),
                              data = scaled_between_year,
                              family = ZINBI,
                           control = gamlss.control(n.cyc = 500))

# model summary and diagnostics
summary(between_year_gamlss)
plot(between_year_gamlss)


# plotting model
predictions <- lpred(between_year_gamlss, what = "mu", type = "response", se.fit = TRUE)
plot_df <- data.frame(scaled_between_year, fv =  predictions$fit, se = predictions$se.fit)

predictions <- predict(between_year_gamlss, what = "mu", type = "response", se.fit = TRUE)
plot_df2 <- data.frame(
  BRDYEAR = plot_df$BRDYEAR, 
  fv = plot_df$fv, 
  predicted_fv = predictions$fit, 
  se = predictions$se.fit
)

brdyear_plot0 <- ggplot(data = plot_df2, aes(x = BRDYEAR)) + 
  geom_point(aes(y = fv), color = "red3", alpha = 0.5) + 
  geom_smooth(aes(y = predicted_fv), method = "loess", se = TRUE, color = "blue", span = 0.5) +
  labs(x = "Breeding Year", y = "Number of Egg Masses") + 
  theme_classic()
brdyear_plot0


# breeding year plot
brdyear_plot <- ggplot(data = plot_df, aes(x = BRDYEAR)) + 
  coord_cartesian(ylim = c(0, 100)) +
  geom_point(aes(y = num_egg_masses), alpha = 0.5) + 
  geom_point(aes(y = fv), color = "red3", alpha = 0.5) +
  geom_smooth(aes(y = fv), method = "loess", se = TRUE) +
  labs(x = "Breeding Year", y = "Number of Egg Masses") +
  theme_classic()
brdyear_plot

# yearly rain plot
yearly_rain_plot <- ggplot(data = plot_df, aes(x = yearly_rain)) + 
  coord_cartesian(ylim = c(0, 100)) +
  geom_point(aes(y = num_egg_masses), alpha = 0.5) + 
  geom_point(aes(y = fv), color = "red3", alpha = 0.5) +
  geom_smooth(aes(y = fv), method = "gam", formula = y ~ s(x, k = 10), se = TRUE) +
  labs(x = "Yearly Rainfall", y = "Number of Egg Masses") +
  theme_classic()
yearly_rain_plot

# percent water plot
percent_water_plot <- ggplot(data = plot_df, aes(x = mean_percent_water)) + 
  coord_cartesian(ylim = c(0, 100)) +
  geom_point(aes(y = num_egg_masses), alpha = 0.5) + 
  geom_point(aes(y = fv), color = "red3", alpha = 0.5) +
  geom_smooth(aes(y = fv), method = "gam", se = TRUE) +
  labs(x = "Percent Open Water", y = "Number of Egg Masses") +
  theme_classic()
percent_water_plot

# percent submergent vegetation plot
percent_submergent_plot <- ggplot(data = plot_df, aes(x = mean_percent_sub)) +
  coord_cartesian(ylim = c(0, 100)) +
  geom_point(aes(y = num_egg_masses), alpha = 0.5) + 
  geom_point(aes(y = fv), color = "red3", alpha = 0.5) +
  geom_smooth(aes(y = fv), method = "gam", se = TRUE) +
  labs(x = "Percent Submergent Vegetation", y = "Number of Egg Masses") +
  theme_classic()
percent_submergent_plot

combined_plot <- brdyear_plot + yearly_rain_plot + percent_water_plot + percent_submergent_plot +
  plot_layout(ncol = 2)
combined_plot

# plotting by watershed
# breeding year plot
brdyear_plot_watershed <- ggplot(data = plot_df, aes(x = BRDYEAR, color = Watershed)) + 
  coord_cartesian(ylim = c(0, 100)) +
  geom_point(aes(y = fv), color = "black", alpha = 0.3) +
  geom_smooth(aes(y = fv), method = "gam", se = TRUE) +
  geom_point(aes(y = num_egg_masses), alpha = 0.5) + 
  labs(x = "Breeding Year", y = "Number of Egg Masses") +
  theme_classic()
brdyear_plot_watershed

# yearly rain plot
yearly_rain_plot_watershed <- ggplot(data = plot_df, aes(x = yearly_rain, color = Watershed)) + 
  coord_cartesian(ylim = c(0, 100)) +
  geom_point(aes(y = fv), color = "black", alpha = 0.3) +
  geom_smooth(aes(y = fv), method = "gam", se = TRUE) +
  geom_point(aes(y = num_egg_masses), alpha = 0.5) + 
  labs(x = "Yearly Rainfall", y = "Number of Egg Masses") +
  theme_classic()
yearly_rain_plot_watershed

# combined plot
combined_plot2 <- brdyear_plot_watershed + yearly_rain_plot_watershed +
  plot_layout(ncol = 2)
combined_plot2

# plotting zero inflation
# plot_nu_df <- data.frame(scaled_between_year, fv = lpred(between_year_gamlss, what = "nu"))

nu_predictions <- lpred(between_year_gamlss, what = "nu", type = "response", se.fit = TRUE)
plot_nu_df <- data.frame(scaled_between_year, fv =  nu_predictions$fit, se = nu_predictions$se.fit) %>% 
  mutate(egg_masses_boolean = if_else(num_egg_masses == 0, 0, 1),
         model_probability = 1 - fv)

# yearly rain plot
yearly_rain_plot <- ggplot(data = plot_nu_df, aes(x = yearly_rain)) + 
  geom_jitter(width = 0.5, height = 0.05, aes(y = egg_masses_boolean), alpha = 0.5) +
  # geom_jitter(width = 0.25, height = 0, aes(y = model_probability), color = "black", alpha = 0.5) +
  geom_smooth(aes(y = model_probability), method = "gam") +
  labs(x = "Yearly Rainfall", y = "Number of Egg Masses") +
  theme_classic()
yearly_rain_plot

yearly_rain_plot <- ggplot(data = plot_nu_df, aes(x = yearly_rain)) + 
  geom_jitter(width = 0.5, height = 0.05, aes(y = egg_masses_boolean), alpha = 0.5) +
  # geom_jitter(width = 0.25, height = 0, aes(y = fv), color = "black", alpha = 0.5) +
  geom_smooth(aes(y = fv), method = "gam") +
  labs(x = "Yearly Rainfall", y = "Number of Egg Masses") +
  theme_classic()
yearly_rain_plot

probability_comparison <- plot_nu_df %>% 
  select(egg_masses_boolean, fv)

#### unused models (saving just in case) ####
#### complete case model ####
complete_case_model <- glmmTMB(num_egg_masses ~ BRDYEAR + 
                                 mean_percent_emerg +
                                 mean_percent_sub +
                                 mean_percent_water +
                                 interpolated_canopy +
                                 yearly_rain + 
                                 mean_max_depth +
                                 max_depth +
                                 WaterTemp +
                                 (1 | LocationID),
                               data = scaled_between_year,
                               ziformula = ~ yearly_rain +
                                 max_depth +
                                 WaterTemp +
                                 mean_salinity:CoastalSite +
                                 max_salinity:CoastalSite +
                                 interpolated_canopy +
                                 (1 | LocationID),
                               family = nbinom2)

complete_case_model_glmer <- glmer(num_egg_masses ~ BRDYEAR + 
                                     I(BRDYEAR^2) +
                                     # mean_percent_emerg + 
                                     # mean_percent_sub +
                                     # mean_percent_water +
                                     # interpolated_canopy +
                                     yearly_rain + 
                                     mean_max_depth +
                                     # max_depth +
                                     # AirTemp +
                                     WaterTemp +
                                     # mean_percent_water +
                                     # mean_salinity:CoastalSite +
                                     # max_salinity:CoastalSite +
                                     interpolated_canopy +
                                     (1 | LocationID),# +
                                   #(1 | LocationID),
                                   data = scaled_between_year,
                                   family = negative.binomial(0.88)) 

#### initial model ####
model1 <- glmmTMB(num_egg_masses ~ BRDYEAR +
                    yearly_rain +
                    (1 | LocationID),
                  data = between_year_data,
                  ziformula = ~0,
                  family = nbinom2) 
summary(model1)

#### GAM model ####
model1.gam <- gam(num_egg_masses ~ s(BRDYEAR) + 
                    s(mean_percent_emerg, k = 3) + # small k so that it doesn't get too wigggly for cover data
                    s(mean_percent_sub, k = 3) +
                    s(mean_percent_water, k = 3) +
                    s(interpolated_canopy, k = 3) +
                    s(yearly_rain) + 
                    max_depth +
                    s(WaterTemp) +
                    max_salinity:as.factor(CoastalSite) +
                    s(Watershed, LocationID, bs = 're'),
                  data = scaled_between_year,
                  family = negbin(0.88))
summary(model1.gam)

#### plotting GAM model ####
# check assumptions
appraise(model1.gam)

# smooth terms
draw(model1.gam)

# all terms
plot(model1.gam, pages = 1, all.terms = TRUE, rug = TRUE)

# just over year
draw(model1.gam, select = 1, rug = FALSE)

# just for rainfall
draw(model1.gam, select = 6)

plot_model(model1.gam, terms = c("CoastalSite","max_salinity"), type = "int")

gam.hp(mod = model1.gam, type = "dev")
permu.gamhp(model1.gam,permutations=100)
plot(gam.hp(mod=model1.gam,type="dev"))

#### lme4 model ####
model1 <- glmer(num_egg_masses ~ I(BRDYEAR - 2009) + 
                    # mean_percent_emerg + 
                    # mean_percent_sub +
                    # mean_percent_water +
                    # interpolated_canopy +
                    yearly_rain + #total annual rainfall
                    #mean_max_depth +
                    # max_depth +
                    #AirTemp +
                    # WaterTemp +
                    # mean_salinity:CoastalSite +
                    # max_salinity:CoastalSite +
                    (1 | Watershed/LocationID),
                  data = between_year_data,
                  family = negative.binomial(0.88) )
summary(model1)


#some gams to deal with non-linearity
#need groups to be factors

between_year_data$Watershed <- factor(between_year_data$Watershed)
between_year_data$LocationID <- factor(between_year_data$LocationID)
between_year_data$CoastalSite <- factor(between_year_data$CoastalSite)



#### plotting complete case model ####

# forest plot (excluding salinity)
plot_model(complete_case_model, terms = c("BRDYEAR", "mean_percent_emerg", "mean_percent_sub", "mean_percent_water", 
                                          "interpolated_canopy", "yearly_rain", "mean_max_depth", "max_depth", "AirTemp", 
                                          "WaterTemp"), vline.color = "slategrey", show.p = TRUE)

# forest plot including salinity (not sure how to just take out the CoastalSiteFALSE variable)
plot_model(complete_case_model, vline.color = "slategrey", show.p = TRUE)

## BB look into whether we can fix the coefficient for coastal site = FALSE
 ## PREFERRED data so no variation and slope = 0?  replace zeros with NAs
 ## Robin will do this.
## BB Need to figure out how to remove the coeffient for coastal site from plot.

# plotting random effects
plot_model(complete_case_model, type = "re", vline.color = "slategrey")




# to test model assumptions
plot_model(complete_case_model,  type = "diag", vline.color = "slategrey")

#### active breeding (actBRD) site model ####

# how many sites active per year?
yearly_active_breeding <- between_year_data %>%
  mutate(breeding = if_else(num_egg_masses > 0, TRUE, FALSE)) %>% 
  filter(breeding == TRUE) %>% 
  select(BRDYEAR, Watershed, LocationID, yearly_rain) %>% 
  group_by(BRDYEAR)%>% 
  summarise(
    actBRD_sites = n(), # active breeding sites each year
    yearly_rain = sum(yearly_rain)
  )
yearly_active_breeding <- as.data.frame(active_breeding)

## plot
plot_active_sites <- yearly_active_breeding %>% ggplot(aes(x=BRDYEAR, y=actBRD_sites))+
  geom_line()+
  geom_point()+
  ylab("number of active breeding sites")
plot_active_sites

## TODO: plotting active breeding with rainfall
  
# create 0/1 binomial for breeding (where num_egg_masses == 0)
scaled_between_year <- scaled_between_year %>% 
  mutate(breeding = if_else(num_egg_masses > 0, TRUE, FALSE))

# logistic regression using glm
BRD_model1 <- glmer(breeding ~ BRDYEAR + 
                    mean_percent_emerg +
                    mean_percent_sub +
                    mean_percent_water +
                    interpolated_canopy +
                    yearly_rain + #total annual rainfall
                    mean_max_depth +
                    max_depth +
                    AirTemp +
                    WaterTemp +
                    mean_salinity:CoastalSite +
                    max_salinity:CoastalSite +
                    (1 | Watershed) +
                    (1 | LocationID),
                  data = scaled_between_year,
                  family = binomial)

summary(BRD_model1)
