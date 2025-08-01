library(tidyverse)
library(dplyr)
library(here)
library(lubridate)
library(reshape2)
library(readxl)

# READ IN DATA ####
# reading in data from spreadsheets
setwd(here::here("code"))
# for V6 data:
raw_data <- read_csv(here::here("data", "CRLF_EGG_RAWDATA_V6.csv"))
rainfall_daily <- read_csv(here::here("data", "cm_daily_rain.csv")) %>% 
  mutate(across(-Water_Year, ~ . * 2.54))
rainfall_yearly <- read_csv(here::here("data", "cm_yearly_rain.csv")) %>% 
  mutate(across(-Water_Year, ~ . * 2.54))
land_cover <- read_csv(here::here("data", "cover_estimates.csv"))
location_type <- read_csv(here::here("data", "CRLF_tblLocations.csv")) %>% 
  select("LocationID", 'Lotic_Lentic', 'WaterRegime') %>% 
  rename(
    water_flow = Lotic_Lentic,
    water_regime = WaterRegime
  ) %>% 
  mutate(water_regime = as.factor(water_regime), water_flow = as.factor(water_flow))
# sun hours data
sun_hours <- read_csv(here::here("data", "Sun_Hours_WY2024.csv")) %>% 
  mutate(str_time = substr(str_time, 1, 10),
         str_time = ymd(str_time),
         BRDYEAR = if_else(month(str_time) > 9, year(str_time) + 1, year(str_time)),
         beginning_WY = ymd(paste0(BRDYEAR - 1, "1001")),
         day_number = as.numeric(str_time - beginning_WY)) %>% 
  select(-OBJECTID, -global_ave, -direct_ave, -diff_ave, -BRDYEAR) %>% 
  group_by(LocationID) %>% 
  mutate(cum_sun_hours = cumsum(dir_dur)) %>% 
  ungroup()

# CLEAN DATA ####

# removing unnecessary columns, making new column for total vegetation (to make sure it adds to 100), making data types more accurate/easier to use
# the DATA variable that this pipe generates has all validated rows and has not been filtered

# filtered data is denoted below this, and uses unfiltered_data as a starting point
unfiltered_data <- raw_data %>% select(-ParkCode, -ProjectCode, -BTime, -TTime, -USGS_ID, -SEASON, -SvyLength, -SvyWidth, -tblEvents.Comments, 
                            -DateEntered, -EventID, -SpeciesID, -WaterDepth, -EggDepth, -Distance, -EggMassStageID, -AS_UTMSOURCE, -AS_UTMZONE, 
                            -GPS_ID, -tblEggCount_CRLF.Comments, -AttachType) %>% 
  filter(OldMass == "FALSE") %>%
  mutate(Date = strptime(SurveyDate, format = "%m/%d/%Y")) %>%
  mutate(Survey_MONTH = as.integer(format(Date, "%m"))) %>%
  mutate(LocationID = as.factor(LocationID), Watershed = as.factor(Watershed), Date = as.Date(Date), Obsv1 = as.factor(Obsv1), Obsv2 = as.factor(Obsv2), Obsv3 = as.factor(Obsv3), 
         Weather = as.factor(Weather), Wind = as.integer(Wind), HabType= as.factor(HabType), SurveyMethodID = as.integer(SurveyMethodID), SalinityMethodID = as.integer(SalinityMethodID), 
         WaterFlowID = as.integer(WaterFlowID), MassID = as.integer(MassID), NumberofEggMasses = as.integer(NumberofEggMasses), BRDYEAR = as.integer(BRDYEAR), DrySurvey = as.logical(DrySurvey)) %>%
  mutate(
    beginningWY = case_when(
      month(Date) > 9 ~ floor_date(Date, unit = "year") + months(9),
      TRUE ~ floor_date(Date, unit = "year") - months(3) # gets the beginning of the water year for each date 
    ) # might be relevant to add that day of water year is zero-indexed, so October 1st is the 0th day of the water year. 
  ) %>%
  mutate(dayOfWY = as.numeric(Date - beginningWY)) %>% # adds column for number of days after the beginning of the water year
  mutate(CoastalSite = if_else((LocationID == "LS01" | LocationID == "LS02" | LocationID == "LS03" | LocationID == "LS04" | LocationID == "LS08" | LocationID == "LS11" | 
                                  LocationID == "RC14" | LocationID == "RC20" | LocationID == "RC21" | LocationID == "RL04" | LocationID == "RL05" | LocationID == "TV06" | LocationID == "WG01"), TRUE, FALSE)) %>% 
  # mutate(WaterSalinity = if_else(!CoastalSite & is.na(WaterSalinity), 0, WaterSalinity)) %>% 
  group_by(EggCountGUID) %>% 
  mutate(obsv_total = sum(!is.na(Obsv1), !is.na(Obsv2), !is.na(Obsv3))) %>% 
  ungroup() %>% 
  select(-Obsv1, -Obsv2, -Obsv3, -OldMass, -SurveyDate) %>% 
  left_join(., rainfall_yearly, join_by(BRDYEAR == Water_Year)) %>% 
  left_join(., land_cover, join_by(LocationID, BRDYEAR == year_numeric)) %>% 
  left_join(., location_type, join_by(LocationID)) %>% 
  rename(
    ground_sub = PercentSubVeg,
    ground_emerg = PercentEmergVeg,
    ground_open_water = PercentOpenWater
  ) %>% 
  mutate(ground_percent_cover_validation = if_else(ground_sub + ground_emerg + ground_open_water == 100 & !is.na(ground_sub) & !is.na(ground_emerg) & !is.na(ground_open_water), TRUE, FALSE),
         interpolated_percent_cover_validation = !is.na(interpolated_sub)) %>% 
  rowwise() %>% 
  mutate(
    mean_percent_sub = if_else(ground_percent_cover_validation == TRUE, if_else(interpolated_percent_cover_validation, mean(c_across(all_of(c("ground_sub", "interpolated_sub"))), na.rm = TRUE), ground_sub), interpolated_sub),
    mean_percent_emerg = if_else(ground_percent_cover_validation == TRUE, if_else(interpolated_percent_cover_validation, mean(c_across(all_of(c("ground_emerg", "interpolated_emerg"))), na.rm = TRUE), ground_emerg), interpolated_emerg),
    mean_percent_water = if_else(ground_percent_cover_validation == TRUE, if_else(interpolated_percent_cover_validation, mean(c_across(all_of(c("ground_open_water", "interpolated_openwater"))), na.rm = TRUE), ground_open_water), interpolated_openwater),
    LocationID = as.factor(LocationID),
    WaterVis_continuous = WaterVis,
    WaterVis = as.integer(if_else(is.na(WaterVis), NA, if_else(WaterVis < 0.3, 0, 1)))
  ) %>% 
  mutate(County = case_when(
    Watershed %in% c("Redwood Creek", "Tennessee Valley", "Rodeo Lagoon", "Oakwood Valley", "Wilkins Gulch", "Easkoot Creek", "Audubon Canyon", "Garden Club Canyon", "Olema Creek") ~ "Marin",
    Watershed %in% c("Milagra Creek", "Laguna Salada", "Kanoff Creek", "West Union", "San Mateo Creek", "San Pedro Creek") ~ "San Mateo",
    TRUE ~ NA_character_
  ))


# ~~~ *** DATA FILTERING *** ~~~ #####

# filter to only include [commented OUT: the 7 watersheds that Darren said had the most data] and only sites that had at least 2 surveys in a given year
# [commented OUT: and only after 2009]
data <- unfiltered_data %>% 
  group_by(LocationID, BRDYEAR) %>% 
  summarize(survey_count_site_yr = n_distinct(EventGUID), .groups = "drop") %>% 
  ungroup() %>% 
  full_join(unfiltered_data, by = c("LocationID" = "LocationID", "BRDYEAR" = "BRDYEAR")) %>% 
  filter(survey_count_site_yr > 1 | (DrySurvey == TRUE)) %>% 
  # filter(Watershed == "Kanoff Creek" | Watershed == "Laguna Salada" | Watershed =="Milagra Creek"|
  #          Watershed == "Redwood Creek" | Watershed == "Rodeo Lagoon" | Watershed=="Tennessee Valley" |
  #          Watershed == "Wilkins Gulch") %>%
  # filter(BRDYEAR > 2009) %>%
  group_by(EventGUID) %>% 
  mutate(NumberofEggMasses = sum(NumberofEggMasses)) %>% 
  ungroup() %>% 
  select(-EggCountGUID, -MassID) %>% 
  distinct(EventGUID, .keep_all = TRUE)  %>% 
  mutate(
    Watershed = droplevels(Watershed),
    LocationID = droplevels(LocationID)
  ) %>% 
  mutate(LocationID = case_when( # fix lowercase LocationIDs (ex: rc01 to RC01)
    LocationID == "rc01" ~ "RC01",
    LocationID == "rc02" ~ "RC02",
    LocationID == "rc03" ~ "RC03",
    TRUE ~ LocationID #keeps other values unchanged
  )) %>%
  mutate(Watershed = case_when(
    Watershed == "Redwood creek" ~ "Redwood Creek",
    TRUE ~ Watershed
  )) %>% 
  mutate( # codes LS02 and LS03 as LS11
    LocationID = case_when(
      LocationID == "LS02" | LocationID == "LS03" ~ "LS11",
      TRUE ~ LocationID
    )
  ) %>%
  filter(BRDYEAR != 2025)

write_csv(data, here::here("data", "filtered_raw_data.csv"))

# ~~~ *** BETWEEN YEAR DATA *** ~~~ ####

between_year_data <- data %>% 
  select(LocationID, BRDYEAR, Watershed, NumberofEggMasses, AirTemp, WaterTemp, MaxD, WaterSalinity, CoastalSite, yearly_rain, yearly_rain_lag,
         ground_sub, ground_emerg, ground_open_water, interpolated_canopy, water_flow, water_regime, DrySurvey, WaterVis, County) %>% 
  group_by(LocationID, BRDYEAR) %>% 
  summarize(
         mean_max_depth = ifelse(all(is.na(MaxD)), NA, mean(MaxD, na.rm = TRUE)),
         max_depth = ifelse(all(is.na(MaxD)), NA, max(MaxD, na.rm = TRUE)),
         # mean_salinity = ifelse(all(is.na(WaterSalinity)), NA, mean(WaterSalinity, na.rm = TRUE)),
         # max_salinity = ifelse(all(is.na(WaterSalinity)), NA, max(WaterSalinity, na.rm = TRUE)),
         AirTemp = ifelse(all(is.na(AirTemp)), NA, mean(AirTemp, na.rm = TRUE)),
         WaterTemp = ifelse(all(is.na(WaterTemp)), NA, mean(WaterTemp, na.rm = TRUE)),
         num_egg_masses = sum(NumberofEggMasses, na.rm = TRUE), 
         mean_percent_sub = ifelse(all(is.na(ground_sub)), NA, mean(ground_sub, na.rm = TRUE)),
         mean_percent_emerg = ifelse(all(is.na(ground_emerg)), NA, mean(ground_emerg, na.rm = TRUE)),
         mean_percent_water = ifelse(all(is.na(ground_open_water)), NA, mean(ground_open_water, na.rm = TRUE)),
         dry_year = if_else(all(DrySurvey), TRUE, FALSE),
         proportion_high_water_vis = sum(WaterVis, na.rm = TRUE) / (n() - sum(is.na(WaterVis))),
         proportion_na_water_vis = sum(is.na(WaterVis)) / n(),
         proportion_high_water_vis = if_else(is.nan(proportion_high_water_vis), 0, proportion_high_water_vis),
         across(everything(), ~first(.))) %>% 
  select(-MaxD, -WaterSalinity, -NumberofEggMasses, -ground_sub, -ground_emerg, -ground_open_water, 
         -DrySurvey, -WaterVis, -AirTemp, -mean_max_depth, -mean_percent_emerg) %>% 
  mutate(
    # mean_salinity = if_else(CoastalSite, mean_salinity, 0),
         # max_salinity = if_else(CoastalSite, max_salinity, 0),
         LocationID = as.factor(LocationID),
         Watershed = as.factor(Watershed)) %>% 
  # select(-max_salinity, -mean_salinity) %>% # deleting salinity from between year
  ungroup()

between_year_data$complete_case <- complete.cases(between_year_data)

## (un)scaling ####
scaled_between_year <- between_year_data %>% 
  filter(complete_case == TRUE) %>% 
  mutate( BRDYEAR_scaled = as.vector(scale(BRDYEAR)),
          mean_percent_sub_scaled = as.vector(scale(mean_percent_sub)),
          # mean_percent_emerg_scaled = as.vector(scale(mean_percent_emerg)),
          mean_percent_water_scaled = as.vector(scale(mean_percent_water)),
          interpolated_canopy_scaled = as.vector(scale(interpolated_canopy)),
          yearly_rain_scaled = as.vector(scale(yearly_rain)),
          # mean_max_depth_scaled = as.vector(scale(max_depth)),
          max_depth_scaled = as.vector(scale(max_depth)),
          # AirTemp_scaled = as.vector(scale(AirTemp)),
          WaterTemp_scaled = as.vector(scale(WaterTemp)),
          yearly_rain_lag_scaled = as.vector(scale(yearly_rain_lag)))

nrow(scaled_between_year)

# for unscaling later
col_means <- between_year_data %>% 
  filter(complete_case == TRUE) %>%
  select(-Watershed, -LocationID, -dry_year, -CoastalSite, -water_flow, -water_regime, -complete_case, -County) %>% 
  colMeans() %>% 
  t() %>% 
  as.data.frame()

write_csv(col_means, here::here("data", "between_year_col_means.csv"))

col_sd <- between_year_data %>% 
  filter(complete_case == TRUE) %>%
  select(-Watershed, -LocationID, -dry_year, -CoastalSite, -water_flow, -water_regime, -complete_case, -County) %>% 
  apply(., 2, sd) %>% 
  t() %>% 
  as.data.frame()

write_csv(col_sd, here::here("data", "between_year_col_sd.csv"))

## lag effect ####
# if we want to include lagged egg masses, this is the code to do so
# holding off for now because it produces so many NA's
## see 4e_brms_lag_betwee_year.R for other filters and use
between_year_data_lagged <- scaled_between_year %>%
  group_by(LocationID) %>%
  arrange(BRDYEAR_scaled) %>%
  mutate(num_egg_masses_lag = lag(num_egg_masses, n = 3)) %>%
  ungroup()

## random slopes data ####
# no covariates, only year 
btw_year_data_random_slopes <- data %>% 
  select(LocationID, BRDYEAR, Watershed, County, NumberofEggMasses) %>% 
  group_by(LocationID, BRDYEAR) %>% 
  summarize(num_egg_masses = sum(NumberofEggMasses, na.rm = TRUE),
            across(everything(), ~first(.))) %>% 
  mutate(
    Watershed = as.factor(Watershed),
    LocationID = as.factor(LocationID),
    County = as.factor(County)
  ) %>% 
  ungroup()

btw_year_data_random_slopes$complete_case <- complete.cases(btw_year_data_random_slopes)

nrow(btw_year_data_random_slopes)

# scale data
scaled_btw_year_data_random_slopes <- btw_year_data_random_slopes %>% 
  filter(complete_case == TRUE) %>% 
  mutate(BRDYEAR_scaled = as.vector(scale(BRDYEAR)))

nrow(scaled_btw_year_data_random_slopes)


## write to CSV ####
write_csv(between_year_data, here::here("data", "between_year_data.csv"))
write_csv(scaled_between_year, here::here("data", "scaled_between_year.csv"))
write_csv(between_year_data_lagged, here::here("data", "lag_between_year_data.csv"))
write_csv(scaled_btw_year_data_random_slopes, here::here("data", "scaled_btw_year_random_slopes.csv"))

## (old) cover comparison ####
# between_year_data_for_cover_comparison <- data %>% 
#   select(LocationID, BRDYEAR, Watershed, NumberofEggMasses, AirTemp, WaterTemp, MaxD, WaterSalinity, CoastalSite, yearly_rain, mean_percent_sub, 
#          mean_percent_emerg, mean_percent_water, ground_sub, ground_emerg, ground_open_water, interpolated_sub, interpolated_emerg, interpolated_openwater) %>% 
#   group_by(LocationID, BRDYEAR) %>% 
#   summarize(
#     mean_max_depth = ifelse(all(is.na(MaxD)), NA, mean(MaxD, na.rm = TRUE)),
#     max_depth = ifelse(all(is.na(MaxD)), NA, max(MaxD, na.rm = TRUE)),
#     mean_salinity = ifelse(all(is.na(WaterSalinity)), NA, mean(WaterSalinity, na.rm = TRUE)),
#     max_salinity = ifelse(all(is.na(WaterSalinity)), NA, max(WaterSalinity, na.rm = TRUE)),
#     AirTemp = ifelse(all(is.na(AirTemp)), NA, mean(AirTemp, na.rm = TRUE)),
#     WaterTemp = ifelse(all(is.na(WaterTemp)), NA, mean(WaterTemp, na.rm = TRUE)),
#     num_egg_masses = sum(NumberofEggMasses, na.rm = TRUE), 
#     mean_percent_sub = ifelse(all(is.na(mean_percent_sub)), NA, mean(mean_percent_sub, na.rm = TRUE)),
#     mean_percent_emerg = ifelse(all(is.na(mean_percent_emerg)), NA, mean(mean_percent_emerg, na.rm = TRUE)),
#     mean_percent_water = ifelse(all(is.na(mean_percent_water)), NA, mean(mean_percent_water, na.rm = TRUE)),
#     mean_ground_sub = ifelse(all(is.na(ground_sub)), NA, mean(ground_sub, na.rm = TRUE)),
#     mean_ground_emerg = ifelse(all(is.na(ground_emerg)), NA, mean(ground_emerg, na.rm = TRUE)),
#     mean_ground_open_water = ifelse(all(is.na(ground_open_water)), NA, mean(ground_open_water, na.rm = TRUE)),
#     mean_interpolated_sub = ifelse(all(is.na(interpolated_sub)), NA, mean(interpolated_sub, na.rm = TRUE)),
#     mean_interpolated_emerg = ifelse(all(is.na(interpolated_emerg)), NA, mean(interpolated_emerg, na.rm = TRUE)),
#     mean_interpolated_open_water = ifelse(all(is.na(interpolated_openwater)), NA, mean(interpolated_openwater, na.rm = TRUE)),
#     across(everything(), ~first(.))) %>% 
#   select(-MaxD, -WaterSalinity, -NumberofEggMasses, -ground_sub, -ground_emerg, -ground_open_water, -interpolated_sub, -interpolated_emerg, -interpolated_openwater) %>% 
#   ungroup()



# ~~~ *** WITHIN YEAR DATA *** ~~~ ####

## survival model ####

# getting rain to date for each day of the water year for every year
rainfall_daily_transposed <- as.data.frame(t(rainfall_daily))
colnames(rainfall_daily_transposed) <- as.character(rainfall_daily_transposed[1,])
rainfall_daily_transposed <- rainfall_daily_transposed[-1,]

for(i in 1:ncol(rainfall_daily_transposed)) {
  x <- rainfall_daily_transposed[,i]
  rainfall_daily_transposed[,i] <- cumsum(x)
}

threshold <- 0.7 # threshold value, we can choose later

first_rainfall <- rownames_to_column(rainfall_daily_transposed, var = "first_rainfall") %>% 
  mutate(first_rainfall = as.integer(substr(first_rainfall, 5, 7))) %>% 
  pivot_longer(cols = -first_rainfall, names_to = "Year", values_to = "rainfall") %>%
  mutate(Year = as.integer(Year)) %>% 
  group_by(Year) %>% 
  filter(rainfall > threshold) %>% 
  arrange(rainfall) %>% 
  slice_head() %>% 
  ungroup()

# getting first breeding entries for each site and year
onset_of_breeding <- data %>% 
  select(LocationID, Watershed, BRDYEAR, dayOfWY, NumberofEggMasses, yearly_rain, yearly_rain_lag, water_flow, water_regime, interpolated_canopy, WaterTemp) %>% 
  group_by(BRDYEAR, LocationID) %>%
  arrange(BRDYEAR, LocationID, dayOfWY) %>% 
  mutate(next_survey = lead(dayOfWY)) %>% 
  filter(NumberofEggMasses > 0) %>% 
  arrange(BRDYEAR, LocationID, dayOfWY) %>% 
  slice(1) %>% 
  mutate(first_breeding = dayOfWY) %>% 
  select(BRDYEAR, LocationID, first_breeding, next_survey) %>% 
  left_join(data, by = c("BRDYEAR", "LocationID")) %>% 
  select(LocationID, Watershed, BRDYEAR, dayOfWY, next_survey, first_breeding, NumberofEggMasses, yearly_rain, yearly_rain_lag, water_flow, water_regime, interpolated_canopy, WaterTemp) %>%
  filter(dayOfWY <= first_breeding) %>% 
  arrange(BRDYEAR, LocationID, dayOfWY) %>% 
  mutate(breeding_status = if_else(first_breeding == dayOfWY, 1, 0),
         rain_to_date = NA,
         next_survey_2 = lead(dayOfWY),
         next_survey = if_else(is.na(next_survey_2), if_else(is.na(next_survey), max(dayOfWY) + 1, next_survey), next_survey_2)) %>% 
  select(-next_survey_2) %>% 
  left_join(first_rainfall, by = c("BRDYEAR" = "Year"))

# this can easily be optimized using the first part of the first_rainfall pipe to put rainfall
# data into tidy format
for (i in 1:nrow(onset_of_breeding)) {
  row <- onset_of_breeding[i,]
  day_of_wy <- as.numeric(row$dayOfWY)
  year <- as.character(row$BRDYEAR)
  y <- c(paste0("day_", day_of_wy))
  z <- rainfall_daily_transposed[rownames(rainfall_daily_transposed) %in% y, ][year]
  onset_of_breeding$rain_to_date[i] <- z
}

onset_of_breeding <- onset_of_breeding %>% 
  mutate(rain_to_date = as.numeric(rain_to_date)) %>% 
  mutate(days_since_first_rain = pmax(0, dayOfWY - first_rainfall)) %>% 
  select(-rainfall, -first_rainfall, -yearly_rain, -yearly_rain_lag, -first_breeding, -NumberofEggMasses)

# adding sun hours
within_year_with_sun_hours <- left_join(onset_of_breeding, sun_hours, by = c("LocationID" = "LocationID", "dayOfWY" = "day_number")) %>% 
  select(-beginning_WY, -str_time)

# write to CSV
write_csv(within_year_with_sun_hours, here::here("data", "onset_of_breeding.csv"))

# looking into earliest and latest onsets for results section
marin_sites <- within_year_with_sun_hours %>% 
  ungroup() %>% 
  filter(breeding_status == 1,
         Watershed %in% c("Redwood Creek", "Rodeo Lagoon", "Tennessee Valley", "Wilkins Gulch")) %>% 
  arrange(dayOfWY)
  summarize(mean_onset = mean(dayOfWY),
            median_onset = median(dayOfWY))

sm_sites <- within_year_with_sun_hours %>% 
  ungroup() %>% 
  filter(breeding_status == 1,
         Watershed %in% c("Milagra Creek", "Kanoff Creek", "Laguna Salada")) %>% 
  summarize(mean_onset = mean(dayOfWY),
            median_onset = median(dayOfWY))

# this is not updated to include days since first rainfall. I'll update it if we end up using timing
# instead of onset
# within year data -- all breeding, not just onset
breeding_timing <- data %>% 
  select(LocationID, Watershed, BRDYEAR, dayOfWY, NumberofEggMasses, yearly_rain, yearly_rain_lag, water_flow, water_regime) %>% 
  group_by(BRDYEAR, LocationID) %>% 
  filter(NumberofEggMasses > 0) %>% 
  select(-NumberofEggMasses) %>% 
  arrange(BRDYEAR, LocationID, dayOfWY) %>% 
  mutate(breeding_date = dayOfWY) %>% 
  mutate(rain_to_date = NA) %>% 
  group_by(BRDYEAR, LocationID) %>% 
  uncount(dayOfWY + 1, .id = "dayOfWY") %>% 
  mutate(dayOfWY = dayOfWY - 1, # need to zero index day of water year to match
         breeding_status = if_else(breeding_date == dayOfWY, 1, 0)) %>% 
  select(-breeding_date) %>% 
  distinct() %>% 
  group_by(LocationID, BRDYEAR, dayOfWY) %>% 
  slice_max(n = 1, breeding_status, with_ties = FALSE) %>% 
  ungroup()

for (i in 1:nrow(breeding_timing)) {
  row <- breeding_timing[i,]
  day_of_wy <- as.numeric(row$dayOfWY)
  year <- as.character(row$BRDYEAR)
  y <- c(paste0("day_", day_of_wy))
  z <- rainfall_daily_transposed[rownames(rainfall_daily_transposed) %in% y, ][year]
  breeding_timing$rain_to_date[i] <- z
}

breeding_timing <- breeding_timing %>% 
  mutate(rain_to_date = as.numeric(rain_to_date))

# adding sun hours
breeding_timing_with_sun_hours <- left_join(breeding_timing, sun_hours, by = c("LocationID" = "LocationID", "dayOfWY" = "day_number")) %>% 
  select(-beginning_WY, -str_time)

# write to CSV
write_csv(breeding_timing_with_sun_hours, here::here("data", "breeding_timing.csv"))

## GAM ####

# getting rain to date for each day of the water year for every year
rainfall_daily_transposed <- as.data.frame(t(rainfall_daily))
colnames(rainfall_daily_transposed) <- as.character(rainfall_daily_transposed[1,])
rainfall_daily_transposed <- rainfall_daily_transposed[-1,]

for(i in 1:ncol(rainfall_daily_transposed)) {
  x <- rainfall_daily_transposed[,i]
  rainfall_daily_transposed[,i] <- cumsum(x)
}

# getting first breeding entries for each site and year
onset_of_breeding <- data %>% 
  select(Watershed, LocationID, BRDYEAR, dayOfWY, NumberofEggMasses, yearly_rain, yearly_rain_lag, water_flow, water_regime) %>% 
  group_by(BRDYEAR, LocationID) %>% 
  filter(NumberofEggMasses > 0) %>% 
  arrange(BRDYEAR, LocationID, dayOfWY) %>% 
  slice(1) %>% 
  mutate(rain_to_date = NA) %>% 
  group_by(BRDYEAR, LocationID) 

for (i in 1:nrow(onset_of_breeding)) {
  row <- onset_of_breeding[i,]
  day_of_wy <- as.numeric(row$dayOfWY)
  year <- as.character(row$BRDYEAR)
  y <- c(paste0("day_", day_of_wy))
  z <- rainfall_daily_transposed[rownames(rainfall_daily_transposed) %in% y, ][year]
  onset_of_breeding$rain_to_date[i] <- z
}

onset_of_breeding <- onset_of_breeding %>% 
  mutate(rain_to_date = as.numeric(rain_to_date))

# adding sun hours
within_year_with_sun_hours <- left_join(onset_of_breeding, sun_hours, by = c("LocationID" = "LocationID", "dayOfWY" = "day_number")) %>% 
  select(-beginning_WY, -str_time)

# write to CSV
write_csv(within_year_with_sun_hours, here::here("data", "onset_of_breeding_gam.csv"))


