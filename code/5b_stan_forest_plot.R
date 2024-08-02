#causal plots
library(ggthemes)

rstan::extract(fit@stanfit)

post <- extract.samples(fit)
PROBS = c(0.1, 0.5, 0.9) ##80 % CIs
#Effect of Rain --> Breach --> DO --> Goby
Rain_MaxDepth <- as_tibble(quantile( with(post, beta_Rain*beta_mean_max_depth), probs = PROBS)) 
Canopy_WaterTemp <- as_tibble(quantile( with(post, beta_Canopy*beta_WaterTemp), probs = PROBS)) 
Year <- as_tibble(quantile( with(post, beta_Year), probs = PROBS)) 
Year_2 <- as_tibble(quantile( with(post, beta_Year_2), probs = PROBS)) 
Rain <- as_tibble(quantile( with(post, beta_Rain), probs = PROBS)) 
Rain_2 <- as_tibble(quantile( with(post, beta_Rain_2), probs = PROBS))
MaxDepth <- as_tibble(quantile( with(post, beta_mean_max_depth), probs = PROBS)) 
Canopy <- as_tibble(quantile( with(post, beta_Canopy), probs = PROBS)) 
Coastal <- as_tibble(quantile( with(post, beta_CoastalSite), probs = PROBS))
WaterTemp <- as_tibble(quantile( with(post, beta_WaterTemp), probs = PROBS))


names<- c("Rain -> MaxDepth","Canopy -> WaterTemp", "Year", "Year_2", "Rain", "Rain_2", "MaxDepth", 
          "Canopy", "Coastal","WaterTemp")
plot.posteriors<-rbind(Rain_MaxDepth,Canopy_WaterTemp, Year, Year_2, Rain, Rain_2,MaxDepth, 
                       Canopy, Coastal,WaterTemp)


plot.posteriors$names <- rep(names, each=3)
#add probabilities names
plot.posteriors$probability <- rep(c("lower", "median", "upper"), times = length(names))

plot.posteriors$row <- rep(1:10, each = 3) #make more general later

plot.posteriors.wide <- plot.posteriors %>% 
  group_by(probability) %>%
  #mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = probability, values_from = value) %>%
  select(-row)

#add codes for positive or negative coefficients
plot.posteriors.wide$effect <- ifelse(
  plot.posteriors.wide$lower<0 & 
    plot.posteriors.wide$median <0 & 
    plot.posteriors.wide$upper <0, 
  "negative", ifelse(plot.posteriors.wide$lower>0 
                     & plot.posteriors.wide$median >0 & 
                       plot.posteriors.wide$upper >0, 
                     "positive",
                     "neutral"))




print(plot.posteriors.wide)

COLORS = c("red", "black", "blue")

plot.posteriors.wide %>%
  mutate(
    names = fct_reorder(names, -median)
  ) %>%
  ggplot(aes(x = names, y = median, color = effect)) +
  #geom_point(effect = c("red", "black", "blue"))) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("Causal Path") + 
  ylab("Causal Effect on Annual Egg Count") +
  scale_color_manual(breaks = c("negative", "neutral", "positive"),
                     values=c("red", "darkgray", "green3")) + 
  coord_flip() +
  scale_x_discrete(limits=rev) +
  theme_few(base_size = 16)

ggsave("Output/forest.plot.png", width = 20, height = 20, units = "cm")
str(rethinking::extract.samples(fit))




