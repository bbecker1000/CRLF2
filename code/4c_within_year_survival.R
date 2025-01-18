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
complete_onset <- unscaled_within_year %>% filter(complete_case == TRUE) %>% select(-complete_case) %>% 
  mutate(time2 = dayOfWY + 1)

# testing collinearity between day of year and sun hours
sun_lm <- lm(cum_sun_hours ~ dayOfWY, data = complete_onset)
sun_resid <- resid(sun_lm)

complete_onset <- complete_onset %>%
  mutate(sun_resid = sun_resid)
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
cox_model_no_random <- coxph(Surv(dayOfWY, time2, breeding_status) ~
                     rain_to_date,
                     # cum_sun_hours + 
                     # sun_resid,
                   data = complete_onset)
summary(cox_model_no_random)

plot(survfit(cox_model_no_random))


# frailty effect: rain to date + cumulative sun hours
cox_model_frailty <- coxph(Surv(dayOfWY, time2, breeding_status) ~ 
                     rain_to_date + 
                     # cum_sun_hours +
                     # sun_resid +
                     frailty(LocationID), 
                   data = complete_onset)
summary(cox_model_frailty)

# coxme with random effects: rain to date + cumulative sun hours
coxme_model <- coxme(Surv(dayOfWY, time2, breeding_status) ~ 
                     rain_to_date +
                     # sun_resid +
                     # cum_sun_hours +
                     (1 | LocationID / BRDYEAR),
                   data = complete_onset)

summary(coxme_model)

# testing assumptions
# to use the cox model, results of test_assumptions must not be significant
test_cox <- cox.zph(coxme_model) # put desired model name here
ggcoxzph(test_cox)
print(test_cox)


# trying to plot survival curves... not working very well atm :(

# extract fixed effects and random effects of coxme model
fixed_effects <- fixef(coxme_model)
random_effects <- ranef(coxme_model)$'LocationID/BRDYEAR'

# get linear predictors of coxme model
linear_pred <- predict(coxme_model, type = "lp")

# get baseline survival curve using cox_model_no_random results
baseline_surv <- survfit(cox_model_no_random)

### part I'm having trouble with ###

cumulative_hazard <- -log(baseline_surv$surv)

# # Adjust the survival curve using the random effects and fixed effects
# survival_estimates <- lapply(names(random_effects), function(cluster) {
#   random_effect <- random_effects[cluster]
#   
#   # Survival curve equation
#   survival_curve <- function(t) {
#     # Find the cumulative hazard for the current time t
#     h_t <- cumulative_hazard[baseline_surv$time <= t]
#     
#     # Survival probability: exp(-sum of hazards) * random effect term
#     exp(-sum(h_t) * exp(fixed_effects %*% complete_onset$rain_to_date + random_effect))
#   }
#   
#   # Apply the survival curve function across all time points
#   survival_values <- sapply(baseline_surv$time, survival_curve)
#   return(survival_values)
# })
# 
# # Combine results into a data frame for plotting
# survival_data <- data.frame(
#   time = rep(baseline_surv$time, length(random_effects)),
#   survival = unlist(survival_estimates),
#   random_effect = rep(names(random_effects), each = length(baseline_surv$time))
# )
# 
# # Generate plot data for a specific random effect location (e.g., "KC01/2012")
# plot_data <- survival_data %>%
#   filter(random_effect == "KC01/2012")


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