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

#### prepping data for analysis ####
setwd(here::here("code"))
#rename file
onset_of_breeding_surv <- read_csv(here::here("data", "onset_of_breeding.csv"))

# scaling covariates
scaled_within_year <- onset_of_breeding_surv %>% 
  mutate(
    BRDYEAR_scaled = as.vector(scale(BRDYEAR)),
    yearly_rain_scaled = as.vector(scale(yearly_rain)),
    rain_to_date_scaled = as.vector(scale(rain_to_date)),
    water_flow = as.factor(water_flow),
    water_regime = as.factor(water_regime), 
    LocationID = as.factor(LocationID),
    cum_sun_hours_scaled = as.vector(scale(cum_sun_hours)),
    dir_dur_scaled = as.vector(scale(dir_dur))) %>% 
  select(-NumberofEggMasses)

# creating a "complete case" column
scaled_within_year$complete_case <- complete.cases(scaled_within_year)
complete_onset <- scaled_within_year %>% filter(complete_case == TRUE) %>% select(-complete_case)

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
# univariate: rain to date
cox_model <- coxph(Surv(dayOfWY, breeding_status) ~ 
                     rain_to_date_scaled, 
                   data = complete_onset)

# univariate: locationID
cox_model <- coxph(Surv(dayOfWY, breeding_status) ~ 
                     LocationID, 
                   data = complete_onset)



# no random effects: rain to date + cumulative sun hours
cox_model_no_random <- coxph(Surv(dayOfWY, breeding_status) ~ 
                     rain_to_date_scaled + 
                     cum_sun_hours_scaled, 
                   data = complete_onset)

summary(cox_model_no_random)
cox.zph(cox_model_no_random)

# frailty effect: rain to date + cumulative sun hours
cox_model_frailty <- coxph(Surv(dayOfWY, breeding_status) ~ 
                     rain_to_date_scaled + 
                     cum_sun_hours_scaled +
                     frailty(LocationID), 
                   data = complete_onset)

summary(cox_model_frailty)
cox.zph(cox_model_frailty)

# coxme with random effects: rain to date + cumulative sun hours
coxme_model <- coxme(Surv(dayOfWY, breeding_status) ~ 
                     rain_to_date_scaled +
                     cum_sun_hours_scaled +
                     (1 | Watershed),
                   data = complete_onset)

summary(coxme_model)
print(coxme_model)
cox.zph(coxme_model)

# summary of model

# testing assumptions
# to use the cox model, results of test_assumptions must not be significant
test_cox <- cox.zph(coxme_model)
test_cox
ggcoxzph(test_cox)
print(test_cox)
plot(test_cox)

ggcoxdiagnostics(cox_model, linear.predictions = TRUE)


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
