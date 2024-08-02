#plot stan models
#doesn't like the missing data function.


library(tidyverse)
library(magrittr)
library(dplyr)
library(purrr)
library(tidyr)
library(cowplot)
library(rethinking)
library(ggrepel)
library(RColorBrewer)
library(brms)
library(rstan)
#library(tidybayes.rethinking)
library(ggthemes)

theme_set(theme_few(base_size = 20))

names(as_tibble(link(fit))) #from 4e, can't extract b/c of missing values code

post <- as.data.frame(link(fit))
nrow(post)
#just get post-warmup values
post <- post[c(3001:6000),]   # if iter = 3000 post[c(2501:5000),]

# Egg predictions
d.mu <- post %>% select(starts_with("mu")) %>%
  pivot_longer(
    cols = starts_with("mu"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "mu")

d.WaterTemp <- post %>% select(starts_with("WaterTemp")) %>%
  pivot_longer(
    cols = starts_with("WaterTemp"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "WaterTemp") %>%
  select(-'case')

d.Depth <- post %>% 
  #select(-contains("_2")) %>%     #remove_2
  select(starts_with("mean_max_depth")) %>%
  pivot_longer(
    cols = starts_with("mean_max_depth"),
    names_to = "case",
    #names_prefix = "wk",
    values_to = "Depth") %>%
  select(-'case')


d.all.post <- bind_cols(d.mu, d.WaterTemp, d.Depth)
## 2024-05-16
## need to add raw data from dat file for unmodeled data
# Year
# Rain
# Temp_2
# BreachDays_2
# Wind
# Micro
# Zone
# Substrate
# Area

d.all.post$Year <- rep(dat_complete$Year, nrow(d.all.post)/nrow(dat_complete))
d.all.post$Year_2 <- rep(dat_complete$Year_2, nrow(d.all.post)/nrow(dat_complete))
d.all.post$Rain <- rep(dat_complete$Rain, nrow(d.all.post)/nrow(dat_complete))
d.all.post$Rain_2 <- rep(dat_complete$Rain_2, nrow(d.all.post)/nrow(dat_complete))
d.all.post$CoastalSite <- rep(dat_complete$CoastalSite, nrow(d.all.post)/nrow(dat_complete))
d.all.post$Location <- rep(dat_complete$Location, nrow(d.all.post)/nrow(dat_complete))
d.all.post$Watershed <- rep(dat_complete$Watershed, nrow(d.all.post)/nrow(dat_complete))

# add raw egg and other covaeriate data
d.all.post$EggMasses <- rep(dat_complete$EggMasses, nrow(d.all.post)/nrow(dat_complete))




#### DO LATER

## now add some original data for plotting on original scale
d.all.post$Year.unscaled <- 
  rep(dat.for.effects.plot$Year.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 
d.all.post$SAV.unscaled <- 
  rep(dat.for.effects.plot$SAV.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 
d.all.post$Temp.unscaled <- 
  rep(dat.for.effects.plot$Temp.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 
d.all.post$DO.unscaled <- 
  rep(dat.for.effects.plot$DO.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 
d.all.post$Rain.unscaled <- 
  rep(dat.for.effects.plot$Rain.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 

d.all.post$Breach.unscaled <- 
  rep(dat.for.effects.plot$Breach.unscaled, nrow(d.all.post)/nrow(dat.for.effects.plot)) 


#and also to dat.plot
dat.plot <- dat

dat.plot$Year.unscaled <- 
  rep(dat.for.effects.plot$Year.unscaled)
dat.plot$SAV.unscaled <- 
  rep(dat.for.effects.plot$SAV.unscaled)
dat.plot$Temp.unscaled <- 
  rep(dat.for.effects.plot$Temp.unscaled) 
dat.plot$DO.unscaled <- 
  rep(dat.for.effects.plot$DO.unscaled) 
dat.plot$Rain.unscaled <- 
  rep(dat.for.effects.plot$Rain.unscaled) 
dat.plot$Breach.unscaled <- 
  rep(dat.for.effects.plot$Breach.unscaled) 




#index for which data sample (cases = 314)
d.all.post$SAMPLE <- rep(1:3000, each = nrow(dat_complete))


#Subset data for the fit lines
# grab random 100 samples near middle of the chain
d <- d.all.post %>% filter(between(SAMPLE, 1900 , 2001) )

 
#Year^2 effects plot
p.Year2 <-  ggplot(data = d, aes(x = Year_2, y = mu)) + 
  geom_point(alpha = 0.05, color = "gray") + #posterior data
  geom_point(data = d, aes(x = Year_2, y = EggMasses), alpha = 0.25, 
             color = "blue") + #raw data
  
  stat_smooth(data = d, method = "lm", 
               formula = y~poly(x,2),
               geom="line", aes(group = SAMPLE), alpha=0.05, linewidth=0.75, color = "red") +
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  xlab("Year")# +
  #guides(x = guide_axis(minor.ticks = TRUE)) #+
  #scale_x_continuous(breaks = c(2000, 2020)) 
p.Year2


hist(d$Year)

#Year effects plot
p.Year <-  ggplot(data = d, aes(x = Year+2009, y = mu)) + 
  geom_point(alpha = 0.05, color = "gray") + #posterior data
  geom_point(data = d, aes(x = Year+2009, y = EggMasses), alpha = 0.25, 
             color = "blue") + #raw data
  
  stat_smooth(data = d, method = "lm", 
              formula = y~poly(x,2),
              geom="line", aes(group = SAMPLE), alpha=0.05, linewidth=0.75, color = "red") +
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  xlab("Year")# +
#guides(x = guide_axis(minor.ticks = TRUE)) #+
#scale_x_continuous(breaks = c(2000, 2020)) 
p.Year




#Rain effects plot
p.Rain <-  ggplot(data = d, aes(x = Rain, y = mu)) + 
  geom_jitter(alpha = 0.05, color = "gray") + #posterior data
  geom_jitter(data = d, aes(x = Rain, y = EggMasses), alpha = 0.25, 
             color = "blue") + #raw data
  stat_smooth(data = d, method = "lm", 
             formula = y~poly(x,2),
              geom="line", aes(group = SAMPLE), alpha=0.05, linewidth=0.75, color = "red") +
  ylim(0,200) + 
  ylab(" ") +
  xlab("Rain")# +
#guides(x = guide_axis(minor.ticks = TRUE)) #+
#scale_x_continuous(breaks = c(2000, 2020)) 
p.Rain


#Rain_2 effects plot
p.Rain2 <-  ggplot(data = d, aes(x = Rain_2, y = mu)) + 
  geom_point(alpha = 0.05, color = "gray") + #posterior data
  geom_point(data = d, aes(x = Rain_2, y = EggMasses), alpha = 0.25, 
              color = "blue") + #raw data
  
  stat_smooth(data = d, method = "lm", 
              formula = y~poly(x,2),
              geom="line", aes(group = SAMPLE), alpha=0.05, linewidth=0.75, color = "red") +
  #geom_smooth(method = "loess", se = FALSE, alpha = 0.25) +
  ylim(0,200) + 
  ylab(" ") +
  xlab("Rain")# +
#guides(x = guide_axis(minor.ticks = TRUE)) #+
#scale_x_continuous(breaks = c(2000, 2020)) 
p.Rain2

hist(d$EggMasses)




#make sure run binary code before logistic model !
p.CoastalSite <- ggplot(data = d, aes(x = as.factor(CoastalSite), y = mu)) + #, group = SAMPLE
  geom_jitter(alpha = 0.05, color = "gray", width = 0.05) + #posterior data
  #stat_smooth (data = d, method = "loess", geom="line", aes(group = SAMPLE), 
  #            alpha=0.05, linewidth=0.75, 
  #           color = "red") +
  # stat_smooth (data = d, method = "lm", 
  #              formula = y~poly(x,2), 
  #              geom="line", aes(group = SAMPLE), 
  #              alpha=0.05, size=0.5) +
  geom_point(data = d, aes(x = CoastalSite, y = EggMasses), 
              alpha = 0.25, color = "blue", width = 0.1) + #raw data
  ylim(0,200) + 
  ylab(" ") +
  xlab("Coastal Site (2=Coastal)") 
p.CoastalSite



p.all.effects <- cowplot::plot_grid(p.Rain2,
                                    p.Year,
                                    p.CoastalSite,
                                    ncol=1, labels="auto", scale = 0.9, 
                                    vjust = 3, hjust = -2.2
)
p.all.effects
ggsave("Output/p.all.effects.jpg", width = 20, height = 35, units = "cm")


#random effects groups
fit %>%
  spread_draws(a_Goby[Zone]) %>%
  median_qi()

# plot RE.
fit %>%
  spread_draws(a_Goby[Zone]) %>%
  median_qi() %>%
  ggplot(aes(y = Zone, x = a_Goby, xmin = .lower, xmax = .upper)) +
  geom_pointinterval()