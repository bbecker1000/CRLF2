# Cover Data Regressions

library(sjPlot)
library(lme4) #glmer won't run today...having issue with matrix as of 2024-04-09
library(glmmTMB) #instead of lme4
library(cowplot)

setwd(here::here("code"))
source("1a_cover_data_prep.R")
# from 1a_cover_data_prep.R

# sites with >2 time points
cover_multi_date <- c("LS01", "RC07", "RC10", "RL02", "TV02")

# make df
cover_multi_df <- cover_data %>%
  filter(LocationID %in% cover_multi_date) 
    
# check filter
unique(cover_multi_df$LocationID)

# binomial TMB for cover
m.cover.bin <- glmmTMB(cbind(OpenWater_percent, 100) ~ year_numeric + (1|LocationID), 
                 family = binomial, 
                 data = cover_multi_df)
m.cover.gaussian <- lmer(OpenWater_percent ~ year_numeric + (1|LocationID), 
                   data = cover_multi_df)

summary(m.cover.gaussian)

plot_model(m.cover.gaussian, type = "diag")

#plot model
p.model <- plot_model(m.cover.gaussian, type = "pred",
           show.data=TRUE) + #wrong scale from raw input data...won't plot
  theme_gray(base_size = 18) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab(" ") +
  xlab("Year") 
p.model


#plot model and raw data
p.model <- plot_model(m.cover.gaussian, type = "pred") +
  theme_gray(base_size = 18) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab(" ") +
  xlab("Year")
  
p.model + 
  geom_point(data=cover_multi_df, 
             aes(x=year_numeric, 
                 y = OpenWater_percent, 
                 group = LocationID)) + 
  geom_smooth(data=cover_multi_df, 
             aes(x=year_numeric, 
                 y = OpenWater_percent, 
                 group = LocationID),
             method = "lm", alpha = 0.2, level = 0.8) +
  ylim(0,100)

#plot raw data
p.raw <- ggplot(cover_multi_df, aes(year_numeric, OpenWater_percent/100, color = LocationID)) +
  geom_point(size = 4) +
  geom_smooth(method = "glm",
              method.args = list(family = binomial()),
              alpha = 0.2, size = 2, level = 0.8) +
  theme_gray(base_size = 18) +
  theme(legend.position = c(0.85, 0.75)) +
  ylab("Percent open water") +
  xlab("Year")
p.raw 

#plot Model and raw data side by side
plot_grid(p.raw, p.model)

#report glmmTMB results and linear interpolation reasonable.

