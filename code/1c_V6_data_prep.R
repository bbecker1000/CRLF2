# follows same process as 1_data_prep, but with CRLF V6 RAW DATA
library(tidyverse)
library(dplyr)

setwd(here::here("code"))

# V4 = current working data
current_data <- read_csv(here::here("data", "filtered_raw_data.csv"))

#### prepping V6 NEW DATA - follows same process as 1_data_prep ####
# reading in V6 NEW data from spreadsheets
V6_raw_data <- read_csv(here::here("data", "CRLF_EGG_RAWDATA_V6.csv"))
rainfall_daily <- read_csv(here::here("data", "cm_daily_rain.csv"))
rainfall_yearly <- read_csv(here::here("data", "cm_yearly_rain.csv"))
land_cover <- read_csv(here::here("data", "cover_estimates.csv"))
location_type <- read_csv(here::here("data", "CRLF_tblLocations.csv")) %>% 
  select("LocationID", 'Lotic_Lentic', 'WaterRegime') %>% 
  rename(
    water_flow = Lotic_Lentic,
    water_regime = WaterRegime
  ) %>% 
  mutate(water_regime = as.factor(water_regime), water_flow = as.factor(water_flow))

# filtered data is denoted below this, and uses unfiltered_data as a starting point
V6_unfiltered_data <- V6_raw_data %>% select(-ParkCode, -ProjectCode, -BTime, -TTime, -USGS_ID, -SEASON, -SvyLength, -SvyWidth, 
                                       -DateEntered, -EventID, -SpeciesID, -WaterDepth, -EggDepth, -Distance, -EggMassStageID, -AS_UTMSOURCE, -AS_UTMZONE, 
                                       -GPS_ID, -AttachType) %>% 
  filter(OldMass == "FALSE") %>%
  mutate(Date = strptime(SurveyDate, format = "%m/%d/%Y")) %>%
  mutate(Survey_MONTH = as.integer(format(Date, "%m"))) %>%
  mutate(LocationID = as.factor(LocationID), Watershed = as.factor(Watershed), Date = as.Date(Date), Obsv1 = as.factor(Obsv1), Obsv2 = as.factor(Obsv2), Obsv3 = as.factor(Obsv3), 
         Weather = as.factor(Weather), Wind = as.integer(Wind), HabType= as.factor(HabType), SurveyMethodID = as.integer(SurveyMethodID), SalinityMethodID = as.integer(SalinityMethodID), 
         WaterFlowID = as.integer(WaterFlowID), MassID = as.integer(MassID), NumberofEggMasses = as.integer(NumberofEggMasses)) %>%
  mutate(
    beginningWY = case_when(
      month(Date) > 9 ~ floor_date(Date, unit = "year") + months(9),
      TRUE ~ floor_date(Date, unit = "year") - months(3) # gets the beginning of the water year for each date 
    ) # might be relevant to add that day of water year is zero-indexed, so October 1st is the 0th day of the water year. 
  ) %>%
  mutate(dayOfWY = as.numeric(Date - beginningWY)) %>% # adds column for number of days after the beginning of the water year
  mutate(CoastalSite = if_else((LocationID == "LS01" | LocationID == "LS02" | LocationID == "LS03" | LocationID == "LS04" | LocationID == "LS08" | LocationID == "LS11" | 
                                  LocationID == "RC14" | LocationID == "RC20" | LocationID == "RC21" | LocationID == "RL04" | LocationID == "RL05" | LocationID == "TV06" | LocationID == "WG01"), TRUE, FALSE)) %>% 
  mutate(WaterSalinity = if_else(!CoastalSite & is.na(WaterSalinity), 0, WaterSalinity)) %>% 
  group_by(EggCountGUID) %>% 
  mutate(obsv_total = sum(!is.na(Obsv1), !is.na(Obsv2), !is.na(Obsv3))) %>% 
  ungroup() %>% 
  select(-Obsv1, -Obsv2, -Obsv3, -OldMass) %>% 
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
    LocationID = as.factor(LocationID)
  )




# ~~~ *** DATA FILTERING


# filter to only include [commented OUT: the 7 watersheds that Darren said had the most data] and only sites that had at least 2 surveys in a given year
# [commented OUT: and only after 2009]
V6_data <- V6_unfiltered_data %>% 
  group_by(LocationID, BRDYEAR) %>% 
  summarize(survey_count_site_yr = n_distinct(EventGUID), .groups = "drop") %>% 
  ungroup() %>% 
  full_join(unfiltered_data, by = c("LocationID" = "LocationID", "BRDYEAR" = "BRDYEAR")) %>% 
  # filter(survey_count_site_yr > 1) %>% 
  # filter(Watershed == "Kanoff Creek" | Watershed == "Laguna Salada" | Watershed =="Milagra Creek"|
  #          Watershed == "Redwood Creek" | Watershed == "Rodeo Lagoon" | Watershed=="Tennessee Valley" |
  #          Watershed == "Wilkins Gulch") %>%
  # filter(BRDYEAR > 2009) %>%
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
  )

write_csv(data, here::here("data", "filtered_V6_raw_data.csv"))

#### comparison #####
# x = V6_data %>% select(LocationID, BRDYEAR, EventGUID, Date)
# y = current_data %>% select(LocationID, BRDYEAR, EventGUID, Date)

diff <- anti_join(V6_data, current_data)
view(diff)

nrow_diff <- nrow(V6_data) - nrow(current_data)
nrow_diff
