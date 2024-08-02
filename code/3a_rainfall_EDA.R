# rainfall calculations that are helpful for EDA but not part of data prep or analysis
# main goal is to chow high correlations between corte madera data and muir woods/half moon bay data

library(tidyverse)
library(dplyr)
library(here)
library(lubridate)
library(readxl)

setwd(here::here("CRLF", "code"))

### ~~~ *** DAILY RAINFALL EDA GRAPHS *** ~~~ ###
rainfall_daily_long <- read_csv(here::here("data", "cm_daily_rain.csv")) %>% 
  pivot_longer(cols = starts_with("day_"), names_to = "day_of_year", values_to =  "rainfall") %>% 
  mutate(day_of_year = as.numeric(gsub("day_", "", day_of_year))) %>% 
  filter(Water_Year > 2009, Water_Year < 2024)

rainfall_cum_long <- rainfall_daily_long %>% 
  group_by(Water_Year) %>% 
  mutate(cum_rain = cumsum(rainfall))

# plot: day of water year vs. daily rainfall
ggplot(data = rainfall_daily_long, aes(x = day_of_year, y = rainfall, color = factor(Water_Year))) + 
  geom_point(alpha = 0.2) + 
  stat_summary(fun = "mean", geom = "point",color = "black", alpha = 0.5) +
  geom_smooth(method = "loess", color = "red4", se = FALSE, size = 2) +
  labs(title = "Daily Rainfall", x = "Day of water year", y = "Rainfall (inches)")

merged_df_new <- left_join(eggTiming_new, rainfall_cum_long, by = "Water_Year")
merged_df_new <- merged_df_new 

ggplot(data = merged_df_new, aes(x = day_of_year, y = cum_rain)) + 
  geom_line()+facet_wrap(~Water_Year)+
  geom_vline(aes(xintercept = merged_df_new$firstEgg),color = "darkolivegreen3")+
  geom_vline(aes(xintercept = merged_df_new$lastEgg), color = "cornflowerblue")


# plot: day of water year vs. cumulative rainfall
ggplot(data = rainfall_cum_long, aes(x = day_of_year, y = cum_rain, color = factor(Water_Year))) + 
  geom_line() + 
  stat_summary(fun = "mean",geom = "line",color = "black", size = 2) +
  labs(title = "Cumulative Rainfall", x = "Day of water year", y = "Cumulative Rainfall (inches)")

merged_df <- left_join(eggTiming_new, rainfall_daily_long, by = "Water_Year")
merged_df

# plot: facet wrap by year vs. daily rainfall
ggplot(data = merged_df, aes(x = day_of_year, y = rainfall)) + 
  geom_point(alpha = 0.1) + geom_smooth(method = "loess", color = "red4", se = FALSE, size = 1.5) +
  scale_y_continuous(trans = "log10")+
  geom_vline(aes(xintercept = merged_df$firstEgg),color = "darkolivegreen3")+
  geom_vline(aes(xintercept = merged_df$lastEgg), color = "cornflowerblue")+
  facet_wrap(~Water_Year)

### ~~~ *** COMPARING YEARLY RAIN ACROSS LOCATIONS *** ~~~ ###

# reading in corte madera data
cm_rain<- read_csv(here::here("data", "cm_yearly_rain.csv"))

# reading in half moon bay data
hmb_rain <- read_csv(here::here("data", "hmb_rain.csv")) %>% 
  select("DATE", "PRCP") %>% 
  mutate(DATE = trimws(DATE), DATE = ym(DATE)) %>% 
  rename(monthly_rain = PRCP, date = DATE) %>% 
  mutate(
    beginningWY = case_when(
      month(date) > 9 ~ floor_date(date, unit = "year") + months(9),
      TRUE ~ floor_date(date, unit = "year") - months(3) # gets the beginning of the water year for each date 
    ),
    Water_Year = as.numeric(year(beginningWY)) + 1) %>% 
  select(-beginningWY) %>% 
  mutate(month = format(date, "%b")) %>% 
  group_by(Water_Year) %>% 
  summarise(yearly_rain = sum(monthly_rain))

# reading in muir woods data
muwo_rain <- read_excel(here::here("data", "muwo_rain.xlsx")) %>% 
  mutate(Oct = as.double(Oct), Nov = as.double(Nov), Dec = as.double(Dec), Jan = as.double(Jan), Feb = as.double(Feb), 
         Mar = as.double(Mar), Apr = as.double(Apr), May = as.double(May), Jun = as.double(Jun), Jul = as.double(Jul), 
         Aug = as.double(Aug), Sep = as.double(Sep), yearly_rain = TOTALS) %>% 
  select(Water_Year, yearly_rain)

# for yearly rain -- all locations

rain_to_compare_wide <- cm_rain %>% 
  left_join(muwo_rain, by = c("Water_Year"), suffix = c("_cm", "_muwo")) %>% 
  left_join(hmb_rain, by = c("Water_Year"), suffix = c("", "_hmb")) %>% 
  rename(half_moon_bay = yearly_rain, corte_madera = yearly_rain_cm, muir_woods = yearly_rain_muwo)

# rain_to_compare_wide <- merge(merge(cm_rain, muwo_rain, all = TRUE), hmb_rain, all = TRUE) %>% 
#   select(Water_Year, yearly_rain, TOTALS, yearly_rain) %>% 
#   rename(corte_madera = yearly_rain_cortemadera, muir_woods = TOTALS, half_moon_bay = yearly_rain_halfmoonbay) %>% 
#   filter(Water_Year > 1998, Water_Year < 2022)

rain_to_compare <- rain_to_compare_wide %>% pivot_longer(cols = 2:4, names_to = "location", values_to = "rainfall")

rain_to_compare_after_missing_data <- rain_to_compare_wide %>% 
  filter(Water_Year > 2003)

# plot rainfall by location
ggplot(data = rain_to_compare, aes(x = Water_Year, y = rainfall, color = location)) + geom_line()

# plot corte madera on x axis, muir woods on y axis, geom_smooth
ggplot(data = rain_to_compare_wide, aes(x = corte_madera, y= muir_woods)) + geom_smooth(method = "lm") + geom_point()

# plot corte madera on x axis, half moon bay on Y axis, geom_smooth
ggplot(data = rain_to_compare_wide, aes(x = corte_madera, y= half_moon_bay)) + geom_smooth(method = "lm") + geom_point()

# calculating correlation coefficient
model_rain = lm (rainfall ~ location, data = rain_to_compare)
summary(model_rain)

ggscatter(rain_to_compare_wide, x = "corte_madera", y = "muir_woods",
          add = "reg.line",cor.coef=TRUE,coor.method=" ",color = "orange")

summary(lm(corte_madera ~ muir_woods, data = rain_to_compare_after_missing_data))

summary(lm(half_moon_bay ~ corte_madera, data = rain_to_compare_after_missing_data))