library(ggplot2)
setwd(here::here("code"))
source("1_data_prep.R")

# plotting number of eggs seen by water visibility

data_wv <- data %>% 
  mutate(WaterVis = as.factor(if_else(WaterVis >= 0.3, "high", "low", missing = "missing"))) %>% 
  group_by(WaterVis) %>% 
  summarize(avg_num_egg_masses = mean(NumberofEggMasses))

vis_plot <- ggplot(data = data_wv, aes(x = WaterVis, y = avg_num_egg_masses)) + geom_bar(position = "dodge", stat = "identity")

# visualizing the variability in water visibility within sites

wv_var <- data %>% 
  select(Watershed, LocationID, BRDYEAR, WaterVis, WaterVis_continuous, NumberofEggMasses) %>% 
  group_by(LocationID) %>% 
  filter(!is.na(WaterVis_continuous))

wv_var_plot <- ggplot(wv_var, aes(x = LocationID, y = WaterVis_continuous)) + 
  geom_jitter(alpha = 0.2, aes(color = Watershed)) +
  geom_boxplot(alpha = 0.8) +
  theme_bw()+
  labs(y = "Water Visibility (feet)", x = "Site")

# might be interesting to look at whether this is related to number of egg masses found?

# generates table of number of NA values for each covariate
data_na_count <- data %>% 
  summarize(NA_open_water = sum(is.na(ground_open_water)),
            NA_emergent_cover = sum(is.na(ground_emerg)),
            NA_submergent_cover = sum(is.na(ground_sub)),
            NA_water_temp = sum(is.na(WaterTemp)),
            NA_water_depth_max = sum(is.na(MaxD)),
            NA_water_depth_avg = sum(is.na(AvgD)),
            NA_water_vis = sum(is.na(WaterVis)),
            NA_salinity = sum(is.na(WaterSalinity)))

### ~~~ *** DATA MANIPULATION (tables that are helpful for creating many graphs) *** ~~~ ###

#summary information of new egg masses count
new_egg <- data|>
  filter(OldMass == "FALSE")
new_egg
typeof(new_egg)
new_total_egg_masses <- new_egg$NumberofEggMasses
#total number of egg masses is 4000
sum(new_total_egg_masses, na.rm = TRUE)
#maximum number of egg masses for a given survey is 25, while the mean is 0.84
summary(new_total_egg_masses,na.rm = TRUE)
new_egg$BRDYEAR <- factor(new_egg$BRDYEAR)

#stats for 1.number of survey, 2. mean number of new egg masses per survey, 3.
#total number of egg masses per watershed per site per year, 4. first and last egg mass detected per year, 
# 5. length of breeding season
statistics <-new_egg|>
  group_by(BRDYEAR,Watershed,LocationID)|>
  summarize(survey_count = n_distinct(EventGUID),
            mean_egg_num = (sum(NumberofEggMasses, na.rm = TRUE)/survey_count),
            total_egg_num =sum(NumberofEggMasses, na.rm = TRUE),
            firstEgg = min(dayOfWY), 
            lastEgg = max(dayOfWY), 
            breedingLength = max(dayOfWY) - min(dayOfWY)
            )

### ~~~ *** PLOTS *** ~~~ ###

# plot total number of eggs per year by site (shortened x-axis for view)
## is this the true sum? double counts?
number_of_eggs_per_year_by_site <-
  ggplot(data = totalEggsPerYear, aes(x = BRDYEAR, y = totalEggs)) +
  geom_point() + facet_wrap(totalEggsPerYear$LocationID) + 
  scale_x_continuous(breaks = seq(min(totalEggsPerYear$BRDYEAR), max(totalEggsPerYear$BRDYEAR), by = 6))
number_of_eggs_per_year_by_site 

# plot timing (days since october 1) of egg laying by year and watershed
eggTimingNoZero <- eggTiming %>% filter(breedingLength > 0) %>% group_by(Watershed, BRDYEAR) %>% 
  summarize(meanFirstEgg = mean(firstEgg), meanLastEgg = mean(lastEgg), meanLength = mean(lastEgg) - mean(firstEgg))

## TODO: change meanLength by SITE 

ggplot(data = eggTimingNoZero, aes(x = BRDYEAR, y = meanLength, color = Watershed)) + geom_line()

## timing graph improvements
ggplot(data = eggTimingNoZero, aes(x = BRDYEAR, y = meanLength, color = Watershed)) + geom_point() + geom_smooth(method="lm", se=F)

## TODO: correlate breeding season length with rainfall
### start of season with rainfall


#table of how many surveys were done each year (not relevant)
survey_count_by_year <- data|>
  group_by(BRDYEAR)|>
  summarize(survey_count=n_distinct(EventGUID))
survey_count_by_year

#plots survey count by year across all watersheds/sites
ggplot(data = survey_count_by_year,aes(x=BRDYEAR,y=survey_count))+geom_point()+geom_smooth()

#table of number of surveys at each watershed per year
abundance_counts <- data %>%
  group_by(Watershed, BRDYEAR) %>%
  summarize(survey_count = n_distinct(EventGUID))
abundance_counts

#Plot showing number of surveys per year
rich_graph <- ggplot(abundance_counts,aes(x=BRDYEAR,y=survey_count,colour=Watershed,group=Watershed))+geom_point()
rich_graph <- rich_graph+ facet_wrap(~Watershed)
rich_graph

#heatmap of survey count per year by watershed
## edit x-axis labels
ggplot(abundance_counts, aes(x = BRDYEAR, y = Watershed, fill = survey_count)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgreen", high = "darkblue") +  # Adjust color gradient
  labs(x = "Breeding Year", y = "Watershed", title = "Abundance Heatmap") +  # Add axis labels and title
  theme_minimal()


#plots of mean number of NEW egg masses per survey across all watersheds of all years
ggplot(data = statistics)+ 
  stat_summary(mapping = aes(x = BRDYEAR, y = mean_egg_num),
               fun = "mean",geom = "point",color = "blue", size= 2)+
  labs(x = "Year", y = "mean number of new egg masses")

## TODO: alternate plot of NEW egg masses = all observations with a TRENDLINE (=mean)
### ALL watersheds

### by watershed



#above graph only looking at 7 target watersheds, and filter by having at least 2 new egg masses per year
statistics_rich <- statistics |>
  filter(Watershed == "Kanoff Creek" | Watershed == "Laguna Salada" | Watershed =="Milagra Creek"|
                  Watershed == "Redwood Creek" | Watershed == "Rodeo Lagoon" | Watershed=="Tennessee Valley" |
                  Watershed == "Wilkins Gulch") |>
  filter(survey_count > 1)
statistics_rich$count <- as.factor(statistics_rich$survey_count)

#plot of total egg masses over time per site at the 7 "rich"watersheds
ggplot(data = statistics_rich, mapping = aes(x = BRDYEAR, y = total_egg_num)) +
  geom_point(aes(size = survey_count), alpha = 1/2) + facet_wrap(~LocationID) +
  labs(x = "Year", y = "total new egg masses")
  
#plot of mean egg masses over time per site at the 7 "rich"watersheds
ggplot(data = statistics_rich, mapping = aes(x = BRDYEAR, y = mean_egg_num)) +
  geom_point(alpha = 1/3) + facet_wrap(~LocationID)+
  labs(x = "Year", y = "mean new egg masses")

#example: number of survey for all sites in Redwood Creek watershed in 2016
data|>
  filter(Watershed=="Redwood Creek")|>
  filter(BRDYEAR == "2016")|>
  group_by(Watershed,LocationID)|>
  summarize(survey_count=n_distinct(EventGUID))

#number of survey count for each number of new egg masses of all watersheds by year
ggplot(data = statistics, mapping = aes(x = total_egg_num)) + 
  geom_freqpoly(mapping = aes(colour = BRDYEAR))+
  scale_x_continuous(limits = c(0, 170))

#table showing number of different sites for all watersheds
number_of_sites_within_watershed <- data|>
  group_by(Watershed)|>
  summarize(distinct_count = n_distinct(LocationID))
number_of_sites_within_watershed

#plot of the above table.
ggplot(number_of_sites_within_watershed,aes(x=Watershed,y=distinct_count))+
  geom_point()

#counting the number of non-NA values for the environmental variables we are looking at 
data$WaterSalinity
non_na_Salinity <- na.omit(data$WaterSalinity)
length(non_na_Salinity)
sum(non_na_Salinity != 0)
#there are 1200 salinity data, only 92 of them is non-zero. 

data |>
  filter(data$WaterSalinity>0)|>
  group_by(Watershed)|>
  summarise(count=n())
#the distribution of these92 non-zero salinity data by watershed.

data$PercentOpenWater
non_na_open <- na.omit(data$PercentOpenWater)
length(non_na_open)
sum(non_na_open != 0)
#there are 4937 non-zero open water data

data$AvgD
non_na_depth <- na.omit(data$AvgD)
length(non_na_depth)
#there are 3404 non-zero open water data

data$MaxD
non_na_depth2 <- na.omit(data$MaxD)
length(non_na_depth2)
#there are 5091 non-zero open water data

#correlation between average depth and max depth is significantly positive 
regression_depth=lm(AvgD~MaxD,data=new_egg)
summary(regression_depth)


ggscatter(new_egg, x = "MaxD", y = "AvgD",
          add = "reg.line",cor.coef=TRUE,coor.method=" ")



# investigating number of observers (obsv_total)
## histogram of obsv_total
ggplot(data = data, aes(x = obsv_total)) + geom_histogram()


