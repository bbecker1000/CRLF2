#survival analysis of first eggs detected
#install.packages(c("survival", "survminer"))

library(tidyverse)
library(survival)
library(survminer)
library(coxme)
library(here)
library(nlme) 
library(gratia)
library(ggplot2)
library(cowplot)
library(riskRegression)
library(adjustedCurves)


#### prepping data for analysis ####
setwd(here::here("code"))
#rename file = unscaled
onset_of_breeding_surv <- read_csv(here::here("data", "onset_of_breeding.csv"))

# scaling covariates
# scaled_within_year <- onset_of_breeding_surv %>%
#   mutate(
#     BRDYEAR_scaled = as.vector(scale(BRDYEAR)),
#     rain_to_date_scaled = as.vector(scale(rain_to_date)),
#     days_since_first_rain_scaled = as.vector(scale(days_since_first_rain)),
#     water_flow = as.factor(water_flow),
#     water_regime = as.factor(water_regime),
#     LocationID = as.factor(LocationID),
#     Watershed = as.factor(Watershed),
#     cum_sun_hours_scaled = as.vector(scale(cum_sun_hours)),
#     dir_dur_scaled = as.vector(scale(dir_dur)))

# for unscaled data
unscaled_within_year <- onset_of_breeding_surv %>% 
  mutate(
    water_flow = as.factor(water_flow),
    water_regime = as.factor(water_regime), 
    LocationID = as.factor(LocationID),
    Watershed = as.factor(Watershed))

# creating a "complete case" column
# scaled_within_year$complete_case <- complete.cases(scaled_within_year)
# complete_onset <- scaled_within_year %>% filter(complete_case == TRUE) %>% select(-complete_case) %>% 
#   mutate(time2 = dayOfWY + 1)

# run only these lines to prep data for unscaled models
unscaled_within_year$complete_case <- complete.cases(unscaled_within_year)
complete_onset <- unscaled_within_year %>% filter(complete_case == TRUE) %>% select(-complete_case)

# sun residuals -- not using
# sun_lm <- lm(cum_sun_hours ~ dayOfWY, data = complete_onset)
# sun_resid <- resid(sun_lm)
# 
# complete_onset <- complete_onset %>%
#   mutate(sun_resid = sun_resid,
#          BRDYEAR = as.factor(BRDYEAR))
#   mutate(breeding_status = max(breeding_status))
# 
# complete_onset <- distinct(complete_onset, diff = paste(LocationID, dayOfWY), .keep_all = "true")


#### *** SURVIVAL MODELS *** ####

#intercept model of the mean
fit.null <- survfit(Surv(dayOfWY, breeding_status) ~ 1, data = onset_of_breeding_surv)
#survival (breeding probability) curves by watershed

fit.rain <- survfit(Surv(dayOfWY, breeding_status) ~ rain_to_date, data = onset_of_breeding_surv)
fit.sunhours <- survfit(Surv(dayOfWY, breeding_status) ~ cum_sun_hours, data = onset_of_breeding_surv)

#pick one to inspect/plot
fit <- fit.null
fit <- fit.watershed
fit <- fit.rain
fit <- fit.sunhours

print(fit)
# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table

#dataframe for whatever
d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower)
d

# plot
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for 
  risk.table = TRUE,       # Add risk table
  risk.table.col = "strata",
  conf.int.style = "step", #or ribbon
  xlab = "Temp",   # customize X axis label.
  ylab = "p(breeding)",
  break.time.by = 1,      # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  surv.median.line = "hv",  # add the median survival pointer.
  fun = "event"             #flips plot to culumative probability
  
)

#### *** COX MODELS *** ####

## choose one
# no random effects: rain to date + cumulative sun hours
# cox_model_no_random <- coxph(Surv(dayOfWY, time2, breeding_status) ~
#                      rain_to_date,
#                      # cum_sun_hours + 
#                      # sun_resid,
#                    data = complete_onset)
# summary(cox_model_no_random)
# 
# plot(survfit(cox_model_no_random))

#### frailty model + plots ####

# group rain_to_date and days_since_first_rain into groups for survival model
# doesn't work now that days since first rain is only positive
onset_grouped <- complete_onset %>% 
  mutate(
    rain_to_date_groups = cut(rain_to_date, 
                              breaks = 10,
                              include.lowest = TRUE),
    days_since_first_rain_groups = cut(days_since_first_rain, 
                                       breaks = 10,
                                       include.lowest = TRUE),
    cum_sun_hours_groups = cut(cum_sun_hours, 
                               breaks = 10,
                               include.lowest = TRUE)
  )

ggplot(data = onset_grouped, aes(x = days_since_first_rain_groups)) +
  geom_bar(aes(fill = as.factor(breeding_status)))
  
ggplot(data = onset_grouped, aes(x = rain_to_date_groups)) +
  geom_bar(aes(fill = as.factor(breeding_status)))

ggplot(data = onset_grouped, aes(x = cum_sun_hours_groups)) +
  geom_bar(aes(fill = as.factor(breeding_status)))


cox_model_frailty <- coxph(Surv(dayOfWY, next_survey, breeding_status) ~ 
                             rain_to_date +
                             days_since_first_rain +
                             frailty(LocationID), 
                           data = complete_onset, 
                           x = TRUE)
summary(cox_model_frailty)

new_data_days <- expand.grid(
  rain_to_date = mean(complete_onset$rain_to_date, na.rm = TRUE),  # Keep rain constant
  days_since_first_rain = seq(min(complete_onset$days_since_first_rain, na.rm = TRUE),
                              max(complete_onset$days_since_first_rain, na.rm = TRUE),
                              length.out = 162)
)

new_data_rain <- expand.grid(
  rain_to_date = seq(min(complete_onset$rain_to_date, na.rm = TRUE),
                     max(complete_onset$rain_to_date, na.rm = TRUE),
                     length.out = 162),
  days_since_first_rain = mean(complete_onset$days_since_first_rain, na.rm = TRUE)
)

new_data_rain$hazard_ratio <- exp(predict(cox_model_frailty, newdata = new_data_rain, type = "risk"))
new_data_days$hazard_ratio <- exp(predict(cox_model_frailty, newdata = new_data_days, type = "risk"))

main_color <- "#49741A"
background <- "#7DC82D"
background2 <- "#91D548"

rain_hazard_plot <- ggplot(new_data_rain, aes(x = rain_to_date, y = hazard_ratio)) +
  geom_line(color = main_color, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = background) +
  labs(x = "Cumulative Rainfall", y = "Hazard Ratio") +
  theme_bw()

time_hazard_plot <- ggplot(new_data_days, aes(x = days_since_first_rain, y = hazard_ratio)) +
  geom_line(color = main_color, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = background) +
  labs(x = "Days Since First Rainfall", y = "Hazard Ratio") +
  theme_bw()

plot_grid(rain_hazard_plot, time_hazard_plot, nrow = 1)

ggforest(cox_model_frailty, data = complete_onset)


predict_fun <- function(...) {
  1 - predictRisk(...)
}

adjusted_curves <- adjustedsurv(
  data = complete_onset,
  variable = "rain_to_date",
  ev_time = "dayOfWY",
  event = "breeding_status",
  method = "direct",
  outcome_model = cox_model_frailty,
  predict_fun = predict_fun
)

plot(adjusted_curves)

# testing assumptions
# to use the cox model, results of test_assumptions must not be significant
test_cox <- cox.zph(cox_model_frailty, transform = "identity") # put desired model name here
ggcoxzph(test_cox)
print(test_cox)

#### coxme + plots -- no longer using ####
# coxme with random effects: rain to date + cumulative sun hours
coxme_model <- coxme(Surv(dayOfWY, time2, breeding_status) ~ 
                     rain_to_date +
                     # sun_resid +
                     # cum_sun_hours +
                     (1 | Watershed),
                   data = complete_onset)

summary(coxme_model)

baseline_hazard <- coxph(Surv(dayOfWY, time2, breeding_status) ~ 1,
                    data = complete_onset)


# trying to plot survival curves... not working very well right now :(

# extract fixed effects and random effects of coxme model
fixed_effects <- fixef(coxme_model)
random_effects <- t(data.frame(ranef(coxme_model)$Watershed))

# get linear predictors of coxme model
linear_pred <- predict(coxme_model, type = "lp")

# get baseline survival curve using baselline_hazard results
baseline_surv <- survfit(baseline_hazard)

### part I'm having trouble with ###

cumulative_hazard <- -log(baseline_surv$surv)

# Adjust the survival curve using the random effects and fixed effects

adjusted_survival <- data.frame(matrix(nrow = length(cumulative_hazard), ncol = length(random_effects))) # needs to be a df -- one column per random effect (site/year), one row per day of WY
colnames(adjusted_survival) <- colnames(random_effects)

# for each day, calculate adjusted survival curve?
# for (i in 1:length(random_effects)) {
#   # For each random effect, adjust the survival curve
#   adjusted_survival[, i] <- baseline_surv * exp(fixed_effects['rain_to_date'] * random_effects[i])  # Adjust based on fixed effects + random effects
# }

# for each locationID/Year combination, iterate through the baseline survival curve and generate adjusted curves
for (i in seq_along(random_effects)) {
  for (t in seq_along(baseline_surv$surv)) {
    # Adjust the baseline survival using fixed and random effects
    hazard_adjustment <- fixed_effects['rain_to_date'] + random_effects[i]
    adjusted_survival[t, i] <- exp(-exp(hazard_adjustment) * cumulative_hazard[t])
  }
}

plot_data <- adjusted_survival %>% 
  mutate(dayOfWY = row.names(.)) %>% 
  pivot_longer(cols = !dayOfWY, names_to = "Watershed", values_to = "hazard") %>% 
  mutate(dayOfWY = as.integer(dayOfWY),
         LocationID = as.factor(Watershed))

ggplot(data = plot_data, aes(x = dayOfWY, y = hazard)) +
  geom_line(aes(color = Watershed))

# geom_point()# baseline_surv * exp(fixed_effects['rain_to_date'] * random_effects[i])
# S(t)=exp⁡(−H(t)⋅exp⁡(βX+u).

### end troublesome part ###

# generate data for plotting -- for future reference: dayOfWY only goes up to 180 because that's the last day of first breeding
plot_data <- data.frame(
  dayOfWY = rep(surv_summary$time, each = length(random_effects)),
  survival = as.vector(survival_estimates),
  random_effect = rep(names(random_effects), length(surv_summary$time))
) %>% 
  group_by(dayOfWY, random_effect) %>% 
  summarize(survival = mean(survival), .groups = "drop") %>%  #aggregate by dayOfWY
  filter(random_effect == "KC01/2012")


# plot!
ggplot(plot_data, aes(x = dayOfWY, y = survival, color = random_effect)) +
  geom_line() +
  labs(
    x = "Day of Water Year",
    y = "Survival Probability",
    color = "LocationID",
    title = "Survival Curves by locationID (Random Effect)"
  ) + theme_minimal()

#### plot model -- forest plot (this one works!) ####

coef <- fixef(coxme_model)

coefficients <- as.data.frame(summary(coxme_model)$coefficients) %>% 
  rename(
    `hazard_ratio` = `exp(coef)`,
    `se` = `se(coef)`
  ) %>% 
  mutate(
    lower_ci = exp(coef - 1.96*se),
    upper_ci = exp(coef + 1.96*se),
    cov = c("Cumulative Rainfall", "Cumulative Sunlight Hours")
  )

ggplot(coefficients, aes(x = hazard_ratio, y = cov)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0) +
  geom_point() +
  labs(x = "Hazard Ratio", y = "Covariate") +
  theme_minimal()

# eda plots -- trying to figure out patterns in the data

ggplot(complete_onset, aes(x = dayOfWY, y = sun_resid)) + 
  geom_smooth(aes(color = LocationID), alpha = 0.4) +
  geom_point(data = complete_onset %>% filter(breeding_status == 1), alpha = 0.25)

ggplot(data = complete_onset, aes(x = cum_sun_hours, y = breeding_status)) +
  geom_smooth()

ggplot(data = complete_onset, aes(x = rain_to_date, y = breeding_status)) +
  geom_smooth()

ggplot(data = complete_onset %>% filter(breeding_status == 1), aes(x = days_since_first_rain)) +
  geom_histogram(binwidth = 5) +
  facet_wrap(~BRDYEAR)

#### GLM: water temperature ####
water_canopy_glm <- glm(WaterTemp ~ interpolated_canopy +
                          dayOfWY +
                          interpolated_canopy:dayOfWY,
                        data=onset_of_breeding_surv)

summary(water_canopy_glm)
plot(water_canopy_glm)

#TODO: trouble-shooting odd residual plots:
#variables on their own
## canopy
water_glm_canopy <- glm(WaterTemp ~ interpolated_canopy,
                        data=onset_of_breeding_surv)
summary(water_glm_canopy)
plot(water_glm_canopy)

## day of water year
water_glm_day <- glm(WaterTemp ~ dayOfWY,
                     data=onset_of_breeding_surv)
summary(water_glm_day)
plot(water_glm_day) # these look funny ...

## interaction
water_glm_int <-  glm(WaterTemp ~ interpolated_canopy:dayOfWY,
                      data=onset_of_breeding_surv)
summary(water_glm_int)
plot(water_glm_int)

#histogram of variables
ggplot(onset_of_breeding_surv, aes(x=WaterTemp))+
  geom_histogram(fill="cornflowerblue")

ggplot(onset_of_breeding_surv, aes(x=dayOfWY))+
  geom_histogram(fill="darkolivegreen4",binwidth=1)

