# plots of covariate raw data between years
library(ggplot2)
library(tidyverse)
library(cowplot)

setwd(here::here("code"))
between_year_data <- read_csv(here::here("data", "between_year_data.csv"))

#### summary data table ####
summary_table_data <- read_csv(here::here("data", "filtered_raw_data.csv"))


#### mean_percent_emerg ####
# emergent vegetation vs. egg masses
emergent_plot <- between_year_data %>% 
  ggplot(aes(x=mean_percent_emerg, y=num_egg_masses, color=Watershed)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  xlab("Emergent Vegetation (%)") +
  ylab(NULL) + 
  theme(legend.position="none")
# emergent vs. time
emergent_BRD_plot <- between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=mean_percent_emerg, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### mean_percent_sub ####
# submerg vs. egg masses
submergent_plot <- between_year_data %>% 
  ggplot(aes(x=mean_percent_sub, y=num_egg_masses, color=Watershed)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  xlab("Submergent Vegetation (%)")  +
  ylab(NULL) + 
  theme(legend.position="none") +
  ylim(c(0,300))

## submerg vs. time
subemergent_BRD_plot <- between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=mean_percent_sub, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### mean_percent_water ####
# open water vs. egg masses
open_water_plot <- between_year_data %>% 
  ggplot(aes(x=mean_percent_water, y=num_egg_masses, color=Watershed)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  xlab("Open Water (%)") +
  ylab(NULL) + 
  theme(legend.position="none")

## open water vs. time
open_water_BRD_plot <- between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=mean_percent_water, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### interpolated_canopy ####
# canopy vs. egg masses
canopy_plot <- between_year_data %>% 
  ggplot(aes(x=interpolated_canopy, y=num_egg_masses, color=Watershed)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  xlab("Canopy (%)") +
  ylab(NULL) + 
  theme(legend.position="none")
## canopy vs. time
canopy_BRD_plot <- between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=interpolated_canopy, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### yearly_rain ####
# yearly rain vs. egg masses
rain_plot <- between_year_data %>% 
  ggplot(aes(x=yearly_rain, y=num_egg_masses, color=Watershed)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  xlab("Rainfall") +
  ylab(NULL) + 
  theme(legend.position="none")
## yearly rain vs. time
rain_BRD_plot <- between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=yearly_rain, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### mean_max_depth ####
# mean_max_depth vs. egg masses
mean_max_depth_plot <- between_year_data %>% 
  ggplot(aes(x=mean_max_depth, y=num_egg_masses, color=Watershed)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  xlab("Mean Max Depth") +
  ylab(NULL)

## mean_max_depth vs. time
mean_max_depth_BRD_plot <- between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=mean_max_depth, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### max_depth ####
# max_depth vs. egg masses
max_depth_plot <- between_year_data %>% 
  ggplot(aes(x=max_depth, y=num_egg_masses, color=Watershed)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  xlab("Max Depth") +
  ylab(NULL)

## max_depth vs. time
max_depth_BRD_plot <- between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=max_depth, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### AirTemp ####
# AirTemp vs. egg masses
AirTemp_plot <- between_year_data %>% 
  ggplot(aes(x=AirTemp, y=num_egg_masses, color=Watershed)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  xlab("Air Temperature") +
  ylab(NULL)

## AirTemp vs. time
AirTemp_BRD_plot <- between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=AirTemp, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### WaterTemp ####
# WaterTemp vs. egg masses
WaterTemp_plot <- between_year_data %>% 
  ggplot(aes(x=WaterTemp, y=num_egg_masses, color=Watershed)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  xlab("Water Temperature") +
  ylab(NULL) + 
  theme(legend.position="none")

## WaterTemp vs. time
WaterTemp_BRD_plot <- between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=WaterTemp, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### mean_salinity:CoastalSite ####
#TODO: do these look right?
coastal_between_year_data <- between_year_data %>% 
  filter(CoastalSite == TRUE)

# mean_salinity vs. egg masses
mean_salinity_plot <- coastal_between_year_data %>% 
  ggplot(aes(x=mean_salinity, y=num_egg_masses)) +
  geom_point(aes(color=Watershed)) +
  geom_smooth(se=FALSE) +
  xlab("Mean Salinity (Coastal)") +
  ylab(NULL) +
  ylim(0,150)

## mean_salinity vs. time
mean_salinity_BRD_plot <- coastal_between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=mean_salinity, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### max_salinity:CoastalSite ####
# max_salinity vs. egg masses
max_salinity_plot <- coastal_between_year_data %>% 
  ggplot(aes(x=max_salinity, y=num_egg_masses)) +
  geom_point(aes(color=Watershed)) +
  geom_smooth(se=FALSE) +
  xlab("Max Salinity (Coastal)") +
  ylab(NULL) +
  ylim(0,150)

## max_salinity vs. time
max_salinity_BRD_plot <- coastal_between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=max_salinity, color=Watershed)) +
  geom_point() +
  geom_smooth()

#### year ####
# year vs. egg masses
egg_BRD_plot <- between_year_data %>% 
  ggplot(aes(x=BRDYEAR, y=num_egg_masses, color=Watershed)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  xlab("Year") +
  ylab(NULL)  + 
  theme(legend.position="none")


#### combined output with egg masses####
# # saved in Ouput folder
# plot_grid(emergent_plot, submergent_plot, open_water_plot, canopy_plot, rain_plot, mean_max_depth_plot, max_depth_plot, AirTemp_plot, WaterTemp_plot, mean_salinity_plot, max_salinity_plot, egg_BRD_plot, nrow=4)


# without salinity or depth
plot_grid(emergent_plot, submergent_plot, open_water_plot, canopy_plot, rain_plot, AirTemp_plot, WaterTemp_plot, egg_BRD_plot, nrow=2)
