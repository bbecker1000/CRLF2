#tidybayes
library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(tidybayes)
library(tidybayes.rethinking)
library(ggplot2)
library(cowplot)

library(rethinking)
library(ggrepel)
library(RColorBrewer)
library(rstan)
#library(gganimate)

theme_set(theme_tidybayes() + panel_border())

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


m<-fit
m %>%
  spread_draws(a_Location[Location]) %>%
  head(10)


m %>%
  recover_types(dat) %>%
  spread_draws(a_Location[Location]) %>%
  head(10)


m %<>% recover_types(dat)


m %>%
  spread_draws(beta_Year, beta_Year_2, beta_CoastalSite, beta_WaterTemp, 
               beta_Rain) %>%
  hist(beta_Rain)



dat_complete %>%
  group_by(Location) %>%
  data_grid(Rain = seq_range(Rain, n = 14)) %>%
  add_linpred_draws(m, ndraws = 25) %>%
  ggplot(aes(x = Rain, y = EggMasses, color = as.factor(Location))) +
  geom_line(aes(y = .linpred, group = paste(Location, .draw)), alpha = .1) +
  geom_point(data = dat_complete) #+
  #scale_color_brewer(palette = "Dark2")

plot(m, depth=2)


dat_complete %>%
  data_grid(Rain = seq_range(Rain, n = 20)) %>%
  add_linpred_draws(m@stanfit) %>%
  ggplot(aes(x = Rain, y = .linpred)) +
  stat_lineribbon(color = "red") +
  scale_fill_brewer(palette = "Greys")
