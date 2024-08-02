#plot chains

chain.data <- as_tibble(extract.samples(fit))



#remove warmup
chain.data <- chain.data[c(3001:6000),]
#add rownames
chain.data$sample <- 1:3000

chain.data$chain <- as.factor(ifelse(
  chain.data$sample<1001, 1, 
  ifelse(chain.data$sample>2000, 3, 2)))

chain.data$sample <- rep(1:1000, times = 3)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_Year))
chain.data %>%
  ggplot(aes(x=sample, y=beta_Year, color = chain)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_Canopy))
chain.data %>%
  ggplot(aes(x=sample, y=beta_Canopy, color = chain)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)


MEAN <- chain.data %>% summarize(MEAN = mean(beta_Year_2))
chain.data %>%
  ggplot(aes(x=sample, y=beta_Year_2, color = chain)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_Rain))
chain.data %>%
  ggplot(aes(x=sample, y=beta_Rain, color = chain)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)

MEAN <- chain.data %>% summarize(MEAN = mean(beta_WaterTemp))
chain.data %>%
  ggplot(aes(x=sample, y=beta_WaterTemp, color = chain)) +
  geom_line(alpha = 0.2) +
  geom_hline(yintercept = as.numeric(MEAN), linetype = 2)


#pairs plots
cowplot::plot_grid(
  
  chain.data %>%
    ggplot(aes(x=beta_Year, y=beta_WaterTemp)) +
    geom_point(color = "blue",alpha = 0.2)
  ,
  
  chain.data %>%
    ggplot(aes(x=beta_Year_2, y=beta_Rain)) +
    geom_point(color = "blue",alpha = 0.2)
  ,
  
  chain.data %>%
    ggplot(aes(x=beta_Rain, y=beta_WaterTemp)) +
    geom_point(color = "blue",alpha = 0.2)
  
)

