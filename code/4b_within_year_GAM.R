library(tidyverse)
library(mgcv)
library(here)
library(nlme) 
library(gratia)
library(ggplot2)
library(cowplot)

#### prepping data for analysis ####
setwd(here::here("code"))
#rename file
onset_of_breeding_surv <- read_csv(here::here("data", "onset_of_breeding_gam.csv"))

# scaling covariates
scaled_within_year <- onset_of_breeding_surv %>% 
  mutate(
    BRDYEAR_scaled = as.vector(scale(BRDYEAR)),
    yearly_rain_scaled = as.vector(scale(yearly_rain)),
    rain_to_date_scaled = as.vector(scale(rain_to_date)),
    # max_depth_scaled = as.vector(scale(MaxD_proportion)),
    # AirTemp_scaled = as.vector(scale(AirTemp)),
    # WaterTemp_scaled = as.vector(scale(WaterTemp)), 
    water_flow = as.factor(water_flow),
    water_regime = as.factor(water_regime), 
    Watershed = as.factor(Watershed),
    LocationID = as.factor(LocationID),
    cum_sun_hours_scaled = as.vector(scale(cum_sun_hours)),
    dir_dur_scaled = as.vector(scale(dir_dur))) %>% 
  select(-NumberofEggMasses)

# creating a "complete case" column
scaled_within_year$complete_case <- complete.cases(scaled_within_year)
complete_onset <- scaled_within_year %>% filter(complete_case == TRUE) %>% select(-complete_case)

#### *** GAM MODEL *** ####
#Generative additive model: first look at onset of breeding with fixed variables
#respectively, and plot to see is the line looks linear or curve.
within_year_gam <- gam(dayOfWY ~ 
                         s(rain_to_date_scaled, by = water_regime) +
                         # s(AirTemp_scaled) +
                         # s(WaterTemp_scaled) +
                         s(BRDYEAR_scaled) + 
                         # s(max_depth_scaled) +
                         water_flow +
                         # water_regime +
                         cum_sun_hours_scaled +
                         # dir_dur_scaled +
                         s(LocationID, Watershed, bs = "re"),
                       data = complete_onset)
summary(within_year_gam)
plot(within_year_gam)
AIC(within_year_gam)

library(gam.hp)
gam.hp(mod=within_year_gam,type="dev")
plot(gam.hp(mod=within_year_gam,type="dev"))


##### plotting GAM model ####
# check assumptions
appraise(within_year_gam)
gam.check(within_year_gam)

# smooth terms
draw(within_year_gam)

# all terms
plot(within_year_gam, pages = 1, all.terms = TRUE, rug = TRUE)

#### plotting using newdata ####
#Mark update: I tried to make all plots more smooth, do they look better/ correct?

#Create a data frame with all predictors
# newdata_AirTemp <- scaled_within_year %>%
#   mutate(
#   max_depth_scaled = mean(scaled_within_year$max_depth_scaled, na.rm = TRUE),
#   AirTemp_scaled = seq(min(scaled_within_year$AirTemp_scaled, na.rm = TRUE), max(scaled_within_year$AirTemp_scaled, na.rm = TRUE), length.out = 1000),
#   WaterTemp_scaled = mean(scaled_within_year$WaterTemp_scaled, na.rm = TRUE),
#   BRDYEAR_scaled = mean(scaled_within_year$BRDYEAR_scaled, na.rm = TRUE),
#   rain_to_date_scaled = mean(scaled_within_year$rain_to_date_scaled, na.rm = TRUE),
#   water_flow = factor(levels(scaled_within_year$water_flow)[1], levels = levels(scaled_within_year$water_flow)),
#   water_regime = factor(levels(scaled_within_year$water_regime)[1], levels = levels(scaled_within_year$water_regime)),
#   Watershed = factor(levels(scaled_within_year$Watershed)[1], levels = levels(scaled_within_year$Watershed)),
#   LocationID = factor(levels(scaled_within_year$LocationID)[1], levels = levels(scaled_within_year$LocationID))
# )

# Generate predictions
predictions <- predict(within_year_gam, newdata = newdata, type = "response", se.fit = TRUE)
plot_df <- data.frame(scaled_within_year, 
                      fv =  predictions$fit, 
                      se = predictions$se.fit,
                      lower = predictions$fit - (1.96 * predictions$se.fit),
                      upper = predictions$fit + (1.96 * predictions$se.fit))


# Create a sequence for AirTemp_scaled
newdata_AirTemp <- with(scaled_within_year, 
                        data.frame(
                          max_depth_scaled = mean(max_depth_scaled, na.rm = TRUE),
                          AirTemp_scaled = seq(min(AirTemp_scaled, na.rm = TRUE), max(AirTemp_scaled, na.rm = TRUE), length.out = 1000),
                          WaterTemp_scaled = mean(WaterTemp_scaled, na.rm = TRUE),
                          BRDYEAR_scaled = mean(BRDYEAR_scaled, na.rm = TRUE),
                          rain_to_date_scaled = mean(rain_to_date_scaled, na.rm = TRUE),
                          water_flow = factor(levels(water_flow)[1], levels = levels(water_flow)),
                          water_regime = factor(levels(water_regime)[1], levels = levels(water_regime)),
                          Watershed = factor(levels(Watershed)[1], levels = levels(Watershed)),
                          LocationID = factor(levels(LocationID)[1], levels = levels(LocationID))
                        )
)

# Generate predictions
predictions_AirTemp <- predict(within_year_gam, newdata = newdata_AirTemp, type = "response", se.fit = TRUE)

# Create a new dataframe for plotting
plot_df_AirTemp <- data.frame(
  AirTemp_scaled = newdata_AirTemp$AirTemp_scaled,
  fv = predictions_AirTemp$fit,
  se = predictions_AirTemp$se.fit,
  lower = predictions_AirTemp$fit - (1.96 * predictions_AirTemp$se.fit),
  upper = predictions_AirTemp$fit + (1.96 * predictions_AirTemp$se.fit)
)

# Plot
Air_temp_plot <- ggplot(data = plot_df_AirTemp, aes(x = AirTemp_scaled, y = fv)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.2) +
  geom_line(color = "blue", size = 1) +
  geom_point(data = scaled_within_year, aes(x = AirTemp_scaled, y = first_breeding), color = "darkblue", alpha = 0.5) +
  labs(x = "Air Temperature (scaled)", y = "Predicted First Breeding") +
  theme_classic() +
  theme(text = element_text(size = 12))
Air_temp_plot

# Create a sequence for rain_to_date_scaled
# I just re-run our gam model, and from the summary rain_to_date is no longer significant?
# Keeping these codes for reference
newdata_rain_to_date <- with(scaled_within_year, 
                             data.frame(
                               max_depth_scaled = mean(max_depth_scaled, na.rm = TRUE),
                               AirTemp_scaled = mean(AirTemp_scaled, na.rm = TRUE),
                               WaterTemp_scaled = mean(WaterTemp_scaled, na.rm = TRUE),
                               BRDYEAR_scaled = mean(BRDYEAR_scaled, na.rm = TRUE),
                               rain_to_date_scaled = seq(min(rain_to_date_scaled, na.rm = TRUE), max(rain_to_date_scaled, na.rm = TRUE), length.out = 1000),
                               water_flow = factor(levels(water_flow)[1], levels = levels(water_flow)),
                               water_regime = factor(levels(water_regime)[1], levels = levels(water_regime)),
                               Watershed = factor(levels(Watershed)[1], levels = levels(Watershed)),
                               LocationID = factor(levels(LocationID)[1], levels = levels(LocationID))
                             )
)

# Generate predictions
predictions_rain_to_date <- predict(within_year_gam, newdata = newdata_rain_to_date, type = "response", se.fit = TRUE)

# Create a new dataframe for plotting
plot_df_rain_to_date <- data.frame(
  rain_to_date_scaled = newdata_rain_to_date$rain_to_date_scaled,
  fv = predictions_rain_to_date$fit,
  se = predictions_rain_to_date$se.fit,
  lower = predictions_rain_to_date$fit - (1.96 * predictions_rain_to_date$se.fit),
  upper = predictions_rain_to_date$fit + (1.96 * predictions_rain_to_date$se.fit)
)

# Plot
rain_plot <- ggplot(data = plot_df_rain_to_date, aes(x = rain_to_date_scaled, y = fv)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.2) +
  geom_line(color = "blue", size = 1) +
  geom_point(data = scaled_within_year, aes(x = rain_to_date_scaled, y = first_breeding), color = "darkblue", alpha = 0.5) +
  labs(x = "Rain to Date (scaled)", y = "Predicted First Breeding") +
  theme_classic() +
  theme(text = element_text(size = 12))
rain_plot


# Create a sequence for WaterTemp_scaled
newdata_WaterTemp <- with(scaled_within_year, 
                          data.frame(
                            max_depth_scaled = mean(max_depth_scaled, na.rm = TRUE),
                            AirTemp_scaled = mean(AirTemp_scaled, na.rm = TRUE),
                            WaterTemp_scaled = seq(min(WaterTemp_scaled, na.rm = TRUE), max(WaterTemp_scaled, na.rm = TRUE), length.out = 1000),
                            BRDYEAR_scaled = mean(BRDYEAR_scaled, na.rm = TRUE),
                            rain_to_date_scaled = mean(rain_to_date_scaled, na.rm = TRUE),
                            water_flow = factor(levels(water_flow)[1], levels = levels(water_flow)),
                            water_regime = factor(levels(water_regime)[1], levels = levels(water_regime)),
                            Watershed = factor(levels(Watershed)[1], levels = levels(Watershed)),
                            LocationID = factor(levels(LocationID)[1], levels = levels(LocationID))
                          )
)

# Generate predictions
predictions_WaterTemp <- predict(within_year_gam, newdata = newdata_WaterTemp, type = "response", se.fit = TRUE)

# Create a new dataframe for plotting
plot_df_WaterTemp <- data.frame(
  WaterTemp_scaled = newdata_WaterTemp$WaterTemp_scaled,
  fv = predictions_WaterTemp$fit,
  se = predictions_WaterTemp$se.fit,
  lower = predictions_WaterTemp$fit - (1.96 * predictions_WaterTemp$se.fit),
  upper = predictions_WaterTemp$fit + (1.96 * predictions_WaterTemp$se.fit)
)

# Plot
WaterTemp_plot <- ggplot(data = plot_df_WaterTemp, aes(x = WaterTemp_scaled, y = fv)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.2) +
  geom_line(color = "blue", size = 1) +
  geom_point(data = scaled_within_year, aes(x = WaterTemp_scaled, y = first_breeding), color = "darkblue", alpha = 0.5) +
  labs(x = "Water Temperature (scaled)", y = "Predicted First Breeding") +
  theme_classic() +
  theme(text = element_text(size = 12))
WaterTemp_plot



# Create a sequence for BRDYEAR_scaled
newdata_BRDYEAR <- with(scaled_within_year, 
                        data.frame(
                          max_depth_scaled = mean(max_depth_scaled, na.rm = TRUE),
                          AirTemp_scaled = mean(AirTemp_scaled, na.rm = TRUE),
                          WaterTemp_scaled = mean(WaterTemp_scaled, na.rm = TRUE),
                          BRDYEAR_scaled = seq(min(BRDYEAR_scaled, na.rm = TRUE), max(BRDYEAR_scaled, na.rm = TRUE), length.out = 1000),
                          rain_to_date_scaled = mean(rain_to_date_scaled, na.rm = TRUE),
                          water_flow = factor(levels(water_flow)[1], levels = levels(water_flow)),
                          water_regime = factor(levels(water_regime)[1], levels = levels(water_regime)),
                          Watershed = factor(levels(Watershed)[1], levels = levels(Watershed)),
                          LocationID = factor(levels(LocationID)[1], levels = levels(LocationID))
                        )
)

# Generate predictions
predictions_BRDYEAR <- predict(within_year_gam, newdata = newdata_BRDYEAR, type = "response", se.fit = TRUE)

# Create a new dataframe for plotting
plot_df_BRDYEAR <- data.frame(
  BRDYEAR_scaled = newdata_BRDYEAR$BRDYEAR_scaled,
  fv = predictions_BRDYEAR$fit,
  se = predictions_BRDYEAR$se.fit,
  lower = predictions_BRDYEAR$fit - (1.96 * predictions_BRDYEAR$se.fit),
  upper = predictions_BRDYEAR$fit + (1.96 * predictions_BRDYEAR$se.fit)
)

# Plot
BRD_plot <- ggplot(data = plot_df_BRDYEAR, aes(x = BRDYEAR_scaled, y = fv)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.2) +
  geom_line(color = "blue", size = 1) +
  geom_point(data = scaled_within_year, aes(x = BRDYEAR_scaled, y = first_breeding), color = "darkblue", alpha = 0.5) +
  labs(x = "Breeding Year (scaled)", y = "Predicted First Breeding") +
  theme_classic() +
  theme(text = element_text(size = 12))
BRD_plot


# Create a sequence for max_depth_scaled, focusing on both water regimes
newdata_max_depth_perennial <- with(scaled_within_year, 
                                    data.frame(
                                      max_depth_scaled = seq(min(max_depth_scaled, na.rm = TRUE), max(max_depth_scaled, na.rm = TRUE), length.out = 1000),
                                      AirTemp_scaled = mean(AirTemp_scaled, na.rm = TRUE),
                                      WaterTemp_scaled = mean(WaterTemp_scaled, na.rm = TRUE),
                                      BRDYEAR_scaled = mean(BRDYEAR_scaled, na.rm = TRUE),
                                      rain_to_date_scaled = mean(rain_to_date_scaled, na.rm = TRUE),
                                      water_flow = factor(levels(water_flow)[1], levels = levels(water_flow)),
                                      water_regime = "perennial", # Set the water_regime to perennial
                                      Watershed = factor(levels(Watershed)[1], levels = levels(Watershed)),
                                      LocationID = factor(levels(LocationID)[1], levels = levels(LocationID))
                                    )
)

newdata_max_depth_seasonal <- with(scaled_within_year, 
                                   data.frame(
                                     max_depth_scaled = seq(min(max_depth_scaled, na.rm = TRUE), max(max_depth_scaled, na.rm = TRUE), length.out = 1000),
                                     AirTemp_scaled = mean(AirTemp_scaled, na.rm = TRUE),
                                     WaterTemp_scaled = mean(WaterTemp_scaled, na.rm = TRUE),
                                     BRDYEAR_scaled = mean(BRDYEAR_scaled, na.rm = TRUE),
                                     rain_to_date_scaled = mean(rain_to_date_scaled, na.rm = TRUE),
                                     water_flow = factor(levels(water_flow)[1], levels = levels(water_flow)),
                                     water_regime = "seasonal", # Set the water_regime to seasonal
                                     Watershed = factor(levels(Watershed)[1], levels = levels(Watershed)),
                                     LocationID = factor(levels(LocationID)[1], levels = levels(LocationID))
                                   )
)

# Generate predictions for both regimes
predictions_max_depth_perennial <- predict(within_year_gam, newdata = newdata_max_depth_perennial, type = "response", se.fit = TRUE)
predictions_max_depth_seasonal <- predict(within_year_gam, newdata = newdata_max_depth_seasonal, type = "response", se.fit = TRUE)

# Create new dataframes for plotting
plot_df_max_depth_perennial <- data.frame(
  max_depth_scaled = newdata_max_depth_perennial$max_depth_scaled,
  fv = predictions_max_depth_perennial$fit,
  se = predictions_max_depth_perennial$se.fit,
  lower = predictions_max_depth_perennial$fit - (1.96 * predictions_max_depth_perennial$se.fit),
  upper = predictions_max_depth_perennial$fit + (1.96 * predictions_max_depth_perennial$se.fit),
  water_regime = "perennial"
)

plot_df_max_depth_seasonal <- data.frame(
  max_depth_scaled = newdata_max_depth_seasonal$max_depth_scaled,
  fv = predictions_max_depth_seasonal$fit,
  se = predictions_max_depth_seasonal$se.fit,
  lower = predictions_max_depth_seasonal$fit - (1.96 * predictions_max_depth_seasonal$se.fit),
  upper = predictions_max_depth_seasonal$fit + (1.96 * predictions_max_depth_seasonal$se.fit),
  water_regime = "seasonal"
)

# Combine the dataframes
plot_df_max_depth <- rbind(plot_df_max_depth_perennial, plot_df_max_depth_seasonal)

# Plot
max_depth_plot <- ggplot(data = plot_df_max_depth, aes(x = max_depth_scaled, y = fv, color = water_regime)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = water_regime), alpha = 0.2) +
  geom_line(size = 1) +
  geom_point(data = scaled_within_year, aes(x = max_depth_scaled, y = first_breeding, color = water_regime), alpha = 0.5) +
  labs(x = "Max Depth (scaled)", y = "Predicted First Breeding") +
  theme_classic() +
  theme(text = element_text(size = 12))
max_depth_plot




# Create a sequence for Watershed
newdata_Watershed <- with(scaled_within_year, 
                          data.frame(
                            max_depth_scaled = mean(max_depth_scaled, na.rm = TRUE),
                            AirTemp_scaled = mean(AirTemp_scaled, na.rm = TRUE),
                            WaterTemp_scaled = mean(WaterTemp_scaled, na.rm = TRUE),
                            BRDYEAR_scaled = mean(BRDYEAR_scaled, na.rm = TRUE),
                            rain_to_date_scaled = mean(rain_to_date_scaled, na.rm = TRUE),
                            water_flow = factor(levels(water_flow)[1], levels = levels(water_flow)),
                            water_regime = factor(levels(water_regime)[1], levels = levels(water_regime)),
                            Watershed = factor(levels(Watershed), levels = levels(Watershed)),
                            LocationID = factor(levels(LocationID)[1], levels = levels(LocationID))
                          )
)

# Generate predictions
predictions_Watershed <- predict(within_year_gam, newdata = newdata_Watershed, type = "response", se.fit = TRUE)

# Create a new dataframe for plotting
plot_df_Watershed <- data.frame(
  Watershed = newdata_Watershed$Watershed,
  fv = predictions_Watershed$fit,
  se = predictions_Watershed$se.fit,
  lower = predictions_Watershed$fit - (1.96 * predictions_Watershed$se.fit),
  upper = predictions_Watershed$fit + (1.96 * predictions_Watershed$se.fit)
)

# Plot
ggplot(data = plot_df_Watershed, aes(x = Watershed, y = fv)) + 
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(x = "Watershed", y = "Predicted First Breeding") +
  theme_classic() +
  theme(text = element_text(size = 12))


plot_grid(Air_temp_plot,WaterTemp_plot, BRD_plot,max_depth_plot,nrow =2)


#Below codes just for reference.
#install this package called "gam.hp", it tells the R2 for each fixed variable.
#somehow this results shows water temperature is not significant? very low R2
library(gam.hp)

#the 2D and 3D plots for GAM (cumulative rain $ water temp)
fit_interaction <- gam(first_breeding ~ te(rain_to_date, WaterTemp, k = c(6, 6)), data = onset_of_breeding_surv)
vis.gam(fit_interaction, view = c("rain_to_date", "WaterTemp"), theta = 30, phi = 30, color = "topo")
vis.gam(fit_interaction, color = 'cm', plot.type = 'contour')
points(onset_of_breeding_surv$rain_to_date, onset_of_breeding_surv$WaterTemp, pch = 16)


#gam model with 3 variables, and MaxD turns out to be insignifant, so we probably
#don't need to include that in later analysis.
fit2_test <- gam(first_breeding ~ s(rain_to_date, k = 10) + s(WaterTemp, k = 10) + s(MaxD, k = 10), data = onset_of_breeding_surv)
summary(fit2_test)
plot(fit2_test, select = 1, pch = 20, se = TRUE, rug = TRUE, residuals = TRUE)
plot(fit2_test, select = 2, pch = 20, se = TRUE, rug = TRUE, residuals = TRUE)
plot(fit2_test, select = 3, pch = 20, se = TRUE, rug = TRUE, residuals = TRUE)
vis.gam(fit2_test, view = c("rain_to_date", "WaterTemp"), theta = 30, phi = 30, color = "heat")