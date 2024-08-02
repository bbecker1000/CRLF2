library(tidyverse)
library(purrr)
library(readxl)
library(lubridate)
library(ggpubr)
library(here)

setwd(here::here("code"))

#reading in corte madera data
marin_rain_folder <- here::here("data", "MMWD_RainfallRecords2")
rain_files <- list.files(path = marin_rain_folder, pattern = "\\.xlsx$", full.names = TRUE)
sheet_name = "CM"

## rain_data_list = list of matrices, each one = 1 year
rain_data_list <- rain_files %>% map( ~ {
  dates_col <- read_xlsx(.x, sheet = sheet_name, range = cell_rows(3:4), col_types = "date") %>% 
    set_names(c("Jul", "Aug", "Sept", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun")) %>% 
    t()
  days_table <- read_xlsx(.x, sheet = sheet_name, range = cell_rows(5:37)) %>% 
    set_names(c("dayOfMonth", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "dayOfMonthAgain")) %>%
    select(-dayOfMonth, -dayOfMonthAgain) %>% 
    t()
  as.data.frame(cbind(dates_col, days_table))
  
})

## CM: combine list of matrices into one single dataframe
cm_rain <- rain_data_list %>% bind_rows()

colnames(cm_rain) <- c("date", c(1:31), "monthly_rain")
cm_rain <- cm_rain %>% mutate(date = as.Date(date), monthly_rain = as.double(monthly_rain)) %>% 
  mutate(
    beginningOfWY = case_when(
      month(date) > 9 ~ floor_date(date, unit = "year") + months(9),
      TRUE ~ floor_date(date, unit = "year") - months(3) # gets the beginning of the water year for each date 
    )
  ) %>% 
  mutate(beginningOfWY = as.Date(beginningOfWY)) %>% 
  mutate(Water_Year = as.numeric(format(beginningOfWY, "%Y")) + 1) %>% 
  select(-beginningOfWY) %>% 
  mutate(month = format(date, "%b"))

## CM: reformat data frame so there is one column per year, one row per day of water year

cm_daily_rain <- cm_rain %>% 
  select(-month, -monthly_rain) %>% 
  pivot_longer(cols = -c(date, Water_Year), names_to = "day", values_to = "daily_rain") %>% 
  mutate(day = as.numeric(day), date = date + (day - 1)) %>% 
  filter(!is.na(daily_rain)) %>% 
  mutate(
    beginningWY = case_when(
      month(date) > 9 ~ floor_date(date, unit = "year") + months(9),
      TRUE ~ floor_date(date, unit = "year") - months(3) # gets the beginning of the water year for each date 
    )
  ) %>%
  mutate(dayOfWY = as.numeric(date - beginningWY)) %>% 
  select(-day, -beginningWY, -date) %>% 
  pivot_wider(names_from = dayOfWY, values_from = daily_rain, values_fn = max) %>% 
  # values_fn is there because sometimes days that do not exist (e.g. Feb 31st) are inputted as 0 instead of NA on the spreadsheets >:(
  # by taking the max, we will ignore those zeroes
  select(Water_Year, order(as.numeric(colnames(.)))) %>% 
  rename_with(~ paste0("day_", .), -1)

# for yearly rain
cm_yearly_rain <- cm_rain %>% 
  group_by(Water_Year) %>% 
  summarise(yearly_rain = sum(monthly_rain))

# write to CSV for daily and yearly rain

write_csv(cm_daily_rain, here::here("data", "cm_daily_rain.csv"))
write_csv(cm_yearly_rain, here::here("data", "cm_yearly_rain.csv"))
