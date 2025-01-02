#survival analysis of first eggs detected
#install.packages(c("survival", "survminer"))

library(tidyverse)
library(survival)
library("survminer")
library("mgcv")
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
#rename file
onset_of_breeding_surv <- read_csv(here::here("data", "onset_of_breeding.csv"))

# scaling covariates
# scaled_within_year <- onset_of_breeding_surv %>% 
#   mutate(
#     BRDYEAR_scaled = as.vector(scale(BRDYEAR)),
#     yearly_rain_scaled = as.vector(scale(yearly_rain)),
#     rain_to_date_scaled = as.vector(scale(rain_to_date)),
#     water_flow = as.factor(water_flow),
#     water_regime = as.factor(water_regime), 
#     LocationID = as.factor(LocationID),
#     Watershed = as.factor(Watershed),
#     cum_sun_hours_scaled = as.vector(scale(cum_sun_hours)),
#     dir_dur_scaled = as.vector(scale(dir_dur))) %>% 
#   select(-NumberofEggMasses)

# for unscaled data
unscaled_within_year <- onset_of_breeding_surv %>% 
  mutate(
    water_flow = as.factor(water_flow),
    water_regime = as.factor(water_regime), 
    LocationID = as.factor(LocationID),
    Watershed = as.factor(Watershed)) %>% 
  select(-NumberofEggMasses)

# creating a "complete case" column
# scaled_within_year$complete_case <- complete.cases(scaled_within_year)
# complete_onset <- scaled_within_year %>% filter(complete_case == TRUE) %>% select(-complete_case)

# run only these lines to prep data for unscaled models
unscaled_within_year$complete_case <- complete.cases(unscaled_within_year)
complete_onset <- unscaled_within_year %>% filter(complete_case == TRUE) %>% select(-complete_case)

# testing collinearity between day of year and sun hours
sun_lm <- lm(cum_sun_hours ~ dayOfWY, data = complete_onset)
sun_resid <- resid(sun_lm)

complete_onset <- complete_onset %>% 
  mutate(sun_resid = sun_resid)
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
cox_model_no_random <- coxph(Surv(dayOfWY, breeding_status) ~ 
                     rain_to_date + 
                     # cum_sun_hours + 
                     sun_resid, 
                   data = complete_onset)
summary(cox_model_no_random)

# frailty effect: rain to date + cumulative sun hours
cox_model_frailty <- coxph(Surv(dayOfWY, breeding_status) ~ 
                     rain_to_date + 
                     # cum_sun_hours +
                     sun_resid +
                     frailty(LocationID), 
                   data = complete_onset)
summary(cox_model_frailty)

# coxme with random effects: rain to date + cumulative sun hours
coxme_model <- coxme(Surv(dayOfWY, breeding_status) ~ 
                     rain_to_date +
                     sun_resid +
                     # cum_sun_hours +
                     # (1 | Watershed),
                     (1 | LocationID),
                   data = complete_onset)

summary(coxme_model)

# testing assumptions
# to use the cox model, results of test_assumptions must not be significant
test_cox <- cox.zph(coxme_model) # put desired model name here
test_cox
ggcoxzph(test_cox)
print(test_cox)
plot(test_cox)



# trying to plot survival curves... not working very well atm :(

# extract fixed effects and random effects of coxme model
fixed_effects <- fixef(coxme_model)
random_effects <- ranef(coxme_model)$LocationID

# get linear predictors of coxme model
linear_pred <- predict(coxme_model, type = "lp")

# get baseline survival curve using cox_model_no_random results
baseline_surv <- survfit(cox_model_no_random)

# modify survival curve for each LocationID
surv_summary <- summary(baseline_surv, times = baseline_surv$time)
survival_estimates <- outer(
  surv_summary$surv,
  exp(random_effects),
  function(base, re) base^re
)

# generate data for plotting -- for future reference: dayOfWY only goes up to 180 because that's the last day of first breeding
plot_data <- data.frame(
  dayOfWY = rep(surv_summary$time, each = length(random_effects)),
  survival = as.vector(survival_estimates),
  random_effect = rep(names(random_effects), length(surv_summary$time))
) %>% 
  group_by(dayOfWY, random_effect) %>% 
  summarize(survival = mean(survival), .groups = "drop") #aggregate by dayOfWY


# plot!
ggplot(plot_data, aes(x = dayOfWY, y = survival, color = random_effect)) +
  geom_line() +
  labs(
    x = "Day of Water Year",
    y = "Survival Probability",
    color = "LocationID",
    title = "Survival Curves by LocationID (Random Effect)"
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
