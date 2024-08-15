library(tidyverse)
library(dplyr)
library(here)
library(lubridate)
library(readxl)
library(here)

setwd(here::here("code"))

# read in cover data file, data filtering (not done yet, need to set data types)
cover_data <- read_xlsx(here::here("data", "canopy_cover.xlsx"), sheet = "Data") %>% 
  mutate(GoogleEarthPhotoDate = as.Date(GoogleEarthPhotoDate))

# making columns for whether it's the first or last date (for later filtering)
cover_data <- cover_data %>%
  group_by(LocationID) %>% 
  mutate(first_date = if_else(GoogleEarthPhotoDate == min(GoogleEarthPhotoDate), TRUE, FALSE),
         last_date = if_else(GoogleEarthPhotoDate == max(GoogleEarthPhotoDate), TRUE, FALSE)) %>% 
  ungroup()

# create a filtered cover data table with only first and last dates (to try to see a linear model for all sites)
cover_data_first_and_last <- cover_data %>% 
  filter(first_date | last_date)

##### visualization purpose #### 
# open water and veg. proportion of the 5 site with continuous data over time.: LS01, RC07, RC10, RL02,TV02
## LS01 
cover_data_LS01 <- cover_data |>
  filter(LocationID == "LS01")
cover_data_LS01$GoogleEarthPhotoDate <- as.Date(cover_data_LS01$GoogleEarthPhotoDate)
cover_data_LS01$year_numeric <- as.numeric(format(cover_data_LS01$GoogleEarthPhotoDate, "%Y"))
cover_data_LS01 <- cover_data_LS01 |>
  select(year_numeric,OpenWater_percent,SubmergentVegetation_percent,EmergentVegetation_percent)
cover_data_LS01_long <- gather(cover_data_LS01, key = "component", value = "proportion", -year_numeric)
ggplot(cover_data_LS01_long,aes(x=year_numeric,y=proportion,fill=component)) + geom_area() + 
  scale_fill_manual(values = c("darkolivegreen3", "steelblue", "darkgreen"))

## RC07
cover_data_RC07 <- cover_data |>
  filter(LocationID == "RC07")
cover_data_RC07$GoogleEarthPhotoDate <- as.Date(cover_data_RC07$GoogleEarthPhotoDate)
cover_data_RC07$year_numeric <- as.numeric(format(cover_data_RC07$GoogleEarthPhotoDate, "%Y"))
cover_data_RC07 <- cover_data_RC07 |>
  select(year_numeric,OpenWater_percent,SubmergentVegetation_percent,EmergentVegetation_percent,TreeCover_percent)
cover_data_RC07_long <- gather(cover_data_RC07, key = "component", value = "proportion", -year_numeric)
ggplot(cover_data_RC07_long,aes(x=year_numeric,y=proportion,fill=component)) + geom_area() + 
  scale_fill_manual(values = c("darkolivegreen3", "steelblue", "grey","darkgreen"))

## RC10
cover_data_RC10 <- cover_data |>
  filter(LocationID == "RC10")
cover_data_RC10$GoogleEarthPhotoDate <- as.Date(cover_data_RC10$GoogleEarthPhotoDate)
cover_data_RC10$year_numeric <- as.numeric(format(cover_data_RC10$GoogleEarthPhotoDate, "%Y"))
cover_data_RC10 <- cover_data_RC10 |>
  select(year_numeric,OpenWater_percent,SubmergentVegetation_percent,EmergentVegetation_percent,TreeCover_percent)
cover_data_RC10_long <- gather(cover_data_RC10, key = "component", value = "proportion", -year_numeric)
ggplot(cover_data_RC10_long,aes(x=year_numeric,y=proportion,fill=component)) + geom_area() + 
  scale_fill_manual(values = c("darkolivegreen3", "steelblue", "grey","darkgreen"))


## RL02
cover_data_RL02 <- cover_data |>
  filter(LocationID == "RL02")
cover_data_RL02$GoogleEarthPhotoDate <- as.Date(cover_data_RL02$GoogleEarthPhotoDate)
cover_data_RL02$year_numeric <- as.numeric(format(cover_data_RL02$GoogleEarthPhotoDate, "%Y"))
cover_data_RL02 <- cover_data_RL02 |>
  select(year_numeric,OpenWater_percent,SubmergentVegetation_percent,EmergentVegetation_percent,TreeCover_percent)
cover_data_RL02_long <- gather(cover_data_RL02, key = "component", value = "proportion", -year_numeric)
ggplot(cover_data_RL02_long,aes(x=year_numeric,y=proportion,fill=component)) + geom_area() + 
  scale_fill_manual(values = c("darkolivegreen3", "steelblue", "grey","darkgreen"))

## TV02
cover_data_TV02 <- cover_data |>
  filter(LocationID == "TV02")
cover_data_TV02$GoogleEarthPhotoDate <- as.Date(cover_data_TV02$GoogleEarthPhotoDate)
cover_data_TV02$year_numeric <- as.numeric(format(cover_data_TV02$GoogleEarthPhotoDate, "%Y"))
cover_data_TV02 <- cover_data_TV02 |>
  select(year_numeric,OpenWater_percent,SubmergentVegetation_percent,EmergentVegetation_percent,TreeCover_percent)
cover_data_TV02_long <- gather(cover_data_TV02, key = "component", value = "proportion", -year_numeric)
ggplot(cover_data_TV02_long,aes(x=year_numeric,y=proportion,fill=component)) + geom_area() + 
  scale_fill_manual(values = c("darkolivegreen3", "steelblue", "grey","darkgreen"))


cover_data$GoogleEarthPhotoDate <- as.Date(cover_data$GoogleEarthPhotoDate)
cover_data$year_numeric <- as.numeric(format(cover_data$GoogleEarthPhotoDate, "%Y"))

# #  linear models for 5 sites
# summary(lm(OpenWater_percent~year_numeric, data = cover_data_LS01))
# plot(lm(OpenWater_percent~year_numeric, data = cover_data_LS01))
# 
# summary(lm(OpenWater_percent~year_numeric, data = cover_data_RC07))
# plot(lm(OpenWater_percent~year_numeric, data = cover_data_RC07))
# 
# summary(lm(OpenWater_percent~year_numeric, data = cover_data_RC10))
# plot(lm(OpenWater_percent~year_numeric, data = cover_data_RC10))  
# 
# summary(lm(OpenWater_percent~year_numeric, data = cover_data_RL02))
# plot(lm(OpenWater_percent~year_numeric, data = cover_data_RL02))

#### ~~~ *** LINEAR INTERPOLATION (all sites) *** ~~~ ####
# 1.INITIAL MANUAL CODE (see below for automated process)
KC01x <- c(0, 12)     
KC01y <- c(50, 35)   

## Apply approx function 
data_approxKC01 <- approx(KC01x, KC01y)        
  
## Draw output of approx function 
plot(data_approxKC01$x,                  
     data_approxKC01$y) 

KC02x <- c(1, 11) 
KC02y <- c(20, 0)   

data_approxKC02 <- approx(KC02x, KC02y)
plot(data_approxKC02$x,                  
     data_approxKC02$y) 


KC03x <- c(1, 12) 
KC03y <- c(20, 0)   

data_approxKC03 <- approx(KC03x, KC03y)
plot(data_approxKC03$x,                  
     data_approxKC03$y) 


LS05x <- c(1, 19) 
LS05y <- c(66.67, 0)   

data_approxLS05 <- approx(LS05x, LS05y)
plot(data_approxLS05$x,                  
     data_approxLS05$y) 

RC10x <- c(2010,2023)
RC10y <- c(100,8)
data_approxRC10 <- approx(RC10x,RC10y, n=14)
plot(data_approxRC10$x, data_approxRC10$y)

Df_new <-data.frame(data_approxRC10)
Df_new
## create function for correct number of year length


# 2. AUTOMATING APPROXIMATION & ADDING TO DATA TABLE
# figuring out year range
start_year <- cover_data %>% 
  filter(first_date == T) %>% 
  select(LocationID, year_numeric) %>% 
  rename(first_year_numeric = year_numeric)

last_year <- cover_data %>% 
  filter(last_date == T) %>% 
  select(LocationID, year_numeric) %>% 
  rename(last_year_numeric = year_numeric)

## creating column "range" = difference between first and last years
range_table <- merge(start_year, last_year) %>% 
  mutate(year_range = last_year_numeric - first_year_numeric)

# start and end year cover values
start_cover <- cover_data %>% 
  filter(first_date == T) %>% 
  select(LocationID, OpenWater_percent, year_numeric)

## checking to see if years align
start_year_check <- merge(start_cover, start_year) %>% 
  mutate(same_start_year = year_numeric == first_year_numeric) # all TRUE!

# start and last year cover values
## start
start_cover <- cover_data %>% 
  filter(first_date == T) %>% 
  select(LocationID, OpenWater_percent, year_numeric) %>% 
  rename(start_open_water = OpenWater_percent)

### checking to see if start years align
check_start <- merge(start_cover, start_year) %>% 
  mutate(same_start_year = year_numeric == first_year_numeric) # all TRUE!

## last
last_cover <- cover_data %>% 
  filter(last_date == T) %>% 
  select(LocationID, OpenWater_percent, year_numeric) %>% 
  rename(last_open_water = OpenWater_percent)

### checking to see if last years align
check_last <- merge(last_cover, last_year) %>% 
  mutate(same_last_year = year_numeric == last_year_numeric) # all TRUE!

## merging start & final cover
#TODO: combine start_cover, last_cover, and range_table$year_range
approx_table <- start_cover %>% 
  full_join(last_cover) ## not doing what i want it to! 
 ## goal: table with site, year, percent open water


#### ~ WORKING METHOD FOR INTERPOLATIONS ~ ####
## BB try using expand.grid
# create all combinations of year and locationID
d1 <- expand.grid(
  year_numeric = c((min(cover_data$year_numeric)):(max(cover_data$year_numeric))),
  LocationID = c(unique(cover_data$LocationID))
  )

## OPEN WATER
#select columns wanted from cover_data
d2OW <- cover_data %>%
  select(LocationID, year_numeric, OpenWater_percent)

#join using unique combinations of locationID and year
d3OW <- as_tibble(left_join(d1, d2OW, by = c('LocationID'='LocationID', 'year_numeric'='year_numeric')))

library(zoo) #for na.approx function
d4OW <- d3OW %>%
  group_by(LocationID) %>%
  mutate(OpenWater_percent = na.approx(OpenWater_percent, na.rm=FALSE)) %>%  #interpolate NA by LocationID
  group_by(LocationID) %>%
  fill(OpenWater_percent, .direction = "downup") %>%                         #fill in leading and trailing NAs
  ungroup() %>%
  # filter(year_numeric > 2009) %>%                                                 #remove pre 2010 data
  rename(interpolated_openwater=OpenWater_percent)

#check data
ggplot(d4OW, aes(year_numeric, interpolated_openwater)) +
  geom_point() +
  geom_line() +
  facet_wrap(.~LocationID)


## EMERGENT VEGETATION
#select columns wanted from cover_data
d2EV <- cover_data %>% 
  select(LocationID, year_numeric, EmergentVegetation_percent)

#join using unique combinations of locationID and year
d3EV <- as_tibble(left_join(d1, d2EV, by=c('LocationID'='LocationID', 'year_numeric'='year_numeric')))

# following BB process from above: interpolating & filling in NAs
d4EV <- d3EV %>% 
  group_by(LocationID) %>% 
  mutate(EmergentVegetation_percent = na.approx(EmergentVegetation_percent, na.rm=FALSE)) %>% 
  group_by(LocationID) %>% 
  fill(EmergentVegetation_percent, .direction="downup") %>% 
  ungroup() %>% 
  # filter(year_numeric > 2009) %>% 
  rename(interpolated_emerg = EmergentVegetation_percent)

# check data
ggplot(d4EV, aes(year_numeric, interpolated_emerg)) +
  geom_point() +
  geom_line() +
  facet_wrap(.~LocationID)


# SUBMERGENT VEGETATION
# select columns from cover_data
d2SV <- cover_data %>% 
  select(LocationID, year_numeric, SubmergentVegetation_percent)

# join using unique combinations of LocationID and year
d3SV <- as_tibble(left_join(d1, d2SV, by=c('LocationID'='LocationID', 'year_numeric'='year_numeric')))

## BB process from above: interpolating & filling in NAs
d4SV <- d3SV %>% 
  group_by(LocationID) %>% 
  mutate(SubmergentVegetation_percent = na.approx(SubmergentVegetation_percent, na.rm = FALSE)) %>% 
  group_by(LocationID) %>% 
  fill(SubmergentVegetation_percent, .direction = 'downup') %>% 
  ungroup() %>% 
  # filter(year_numeric >2009) %>% 
  rename(interpolated_sub = SubmergentVegetation_percent)

# check data
ggplot(d4SV, aes(year_numeric, interpolated_sub)) +
  geom_point() +
  geom_line() +
  facet_wrap(.~LocationID)


# CANOPY (TREE) COVER
## canopy cover boolean
canopy_boolean <- cover_data %>% 
  mutate(canopy = if_else(TreeCover_percent > 0, TRUE, FALSE)) %>% 
  group_by(LocationID, canopy) %>% 
  summarize() %>% 
  arrange(LocationID, desc(canopy)) %>% 
  distinct(LocationID, .keep_all = TRUE)

## canopy cover linear interpolation
# select columns from cover_data
d2CC <- cover_data %>% 
  select(LocationID, year_numeric, TreeCover_percent)

# join using unique combinations of LocationID and year
d3CC <- as_tibble(left_join(d1, d2CC, by=c('LocationID'='LocationID', 'year_numeric'='year_numeric')))

## BB process from above: interpolating & filling in NAs
d4CC <- d3CC %>% 
  group_by(LocationID) %>% 
  mutate(TreeCover_percent = na.approx(TreeCover_percent, na.rm = FALSE)) %>% 
  group_by(LocationID) %>% 
  fill(TreeCover_percent, .direction='downup') %>% 
  ungroup() %>% 
  # filter(year_numeric > 2009) %>% 
  rename(interpolated_canopy = TreeCover_percent)

# check data
ggplot(d4CC, aes(year_numeric, interpolated_canopy)) +
  geom_point() +
  geom_line() +
  facet_wrap(.~LocationID)



## COMBINING ALL  COVER DATA
cover_estimates <- d4EV %>% 
  merge(d4OW) %>% # open water
  merge(d4SV) %>%  # submergent veg
  merge(d4CC) # canopy cover (data)

# converting into csv
write_csv(cover_estimates, here::here("data", "cover_estimates.csv"))

# combined plot
cover_colors <- c('interpolated_openwater'='navy', 'interpolated_emerg'="#00a2ad", 'interpolated_submerg'= '#65a300', 'interpolated_canopy'='#b86500')

cover_estimate_plot <- ggplot(cover_estimates, aes(x=year_numeric)) +
  geom_line(aes(y=interpolated_openwater), color="navy")+
  geom_line(aes(y=interpolated_emerg), color="#00a2ad") +
  geom_line(aes(y=interpolated_sub), color='#65a300')+
  geom_line(aes(y=interpolated_canopy), color='#b86500')+
  facet_wrap(.~LocationID)+
  labs(y="Percent Cover",
     x= "Year",
     color = "Legend") 
  # + scale_color_manual(labels = c("Open Water", "Emergent", "Submergent", 'Canopy'), 
  #                    values = c("navy", "#00a2ad", "#65a300", "#b86500"))

#TODO: add legend to this plot

cover_estimate_plot


