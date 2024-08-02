#mvgam

library(mvgam)
data("portal_data")

portal_data %>%
  
  # mvgam requires a 'time' variable be present in the data to index
  # the temporal observations. This is especially important when tracking 
  # multiple time series. In the Portal data, the 'moon' variable indexes the
  # lunar monthly timestep of the trapping sessions
  dplyr::mutate(time = moon - (min(moon)) + 1) %>%
  
  # We can also provide a more informative name for the outcome variable, which 
  # is counts of the 'PP' species (Chaetodipus penicillatus) across all control
  # plots
  dplyr::mutate(count = PP) %>%
  
  # The other requirement for mvgam is a 'series' variable, which needs to be a
  # factor variable to index which time series each row in the data belongs to.
  # Again, this is more useful when you have multiple time series in the data
  dplyr::mutate(series = as.factor('PP')) %>%
  
  # Select the variables of interest to keep in the model_data
  dplyr::select(series, year, time, count, mintemp, ndvi) -> model_data


plot_mvgam_series(data = model_data, series = 1, y = 'count')


model_data %>%
  
  # Create a 'year_fac' factor version of 'year'
  dplyr::mutate(year_fac = factor(year)) -> model_data


model1 <- mvgam(count ~ s(year_fac, bs = 're') - 1,
                family = poisson(),
                data = model_data)

summary(model1)
code(model1)
plot(model1, type = 're')
mcmc_plot(object = model1,
          variable = 'betas',
          type = 'areas')
pp_check(object = model1)
plot(model1, type = 'forecast')

model_data %>% 
  dplyr::filter(time <= 160) -> data_train 
model_data %>% 
  dplyr::filter(time > 160) -> data_test
model1b <- mvgam(count ~ s(year_fac, bs = 're') - 1,
                 family = poisson(),
                 data = data_train,
                 newdata = data_test)
plot(model1b, type = 'forecast', newdata = data_test)



model2 <- mvgam(count ~ s(year_fac, bs = 're') + 
                  ndvi - 1,
                family = poisson(),
                data = data_train,
                newdata = data_test)

summary(model2)

beta_post <- as.data.frame(model2, variable = 'betas')
hist(beta_post$ndvi,
     xlim = c(-1 * max(abs(beta_post$ndvi)),
              max(abs(beta_post$ndvi))),
     col = 'darkred',
     border = 'white',
     xlab = expression(beta[NDVI]),
     ylab = '',
     yaxt = 'n',
     main = '',
     lwd = 2)
abline(v = 0, lwd = 2.5)

conditional_effects(model2)

model3 <- mvgam(count ~ s(time, bs = 'bs', k = 15) + 
                  ndvi,
                family = poisson(),
                data = data_train,
                newdata = data_test)
summary(model3)
conditional_effects(model3, type = 'link')
plot(model3, type = 'forecast', newdata = data_test)

model4 <- mvgam(count ~ s(ndvi, k = 6),
                family = poisson(),
                data = data_train,
                newdata = data_test,
                trend_model = 'AR1')
plot(model4, type = 'forecast', newdata = data_test)
plot(model4, type = 'trend', newdata = data_test)
loo_compare(model3, model4)
