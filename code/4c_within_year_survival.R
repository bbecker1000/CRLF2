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

unscaled_within_year <- onset_of_breeding_surv %>% 
  mutate(
    water_flow = as.factor(water_flow),
    water_regime = as.factor(water_regime), 
    LocationID = as.factor(LocationID),
    Watershed = as.factor(Watershed)) %>% 
  select(-WaterTemp) # for complete cases, lots of NA's

unscaled_within_year$complete_case <- complete.cases(unscaled_within_year)
complete_onset <- unscaled_within_year %>% 
  filter(complete_case == TRUE,
         next_survey > dayOfWY) %>% 
  select(-complete_case)

summary_stats_data <- complete_onset %>% 
  filter(breeding_status == 1) %>% 
  mutate(beginningWY = ymd(paste0(BRDYEAR - 1, "-10-01")),
         first_breeding = beginningWY + days(dayOfWY))

ggplot(summary_stats_data, aes(x = dayOfWY)) +
  geom_histogram()

summary(summary_stats_data$dayOfWY)

#### frailty model + plots ####

# group rain_to_date into groups for survival model
onset_grouped <- complete_onset %>% 
  mutate(
    rain_to_date_groups = cut(rain_to_date, 
                              # breaks = quantile(rain_to_date, probs = 0:5 / 5, na.rm = TRUE),
                              breaks = c(0, 25, 50, 75, 100, 125),
                              # breaks = 5,
                              labels = c("0 - 25 cm", "25 - 50 cm", "50 - 75 cm", "75 - 100 cm", "100 - 125 cm"),
                              include.lowest = TRUE),
    canopy_groups = cut(interpolated_canopy,
                           # breaks = quantile(WaterTemp, probs = 0:5 / 5, na.rm = TRUE),
                           breaks = c(0, 1, 25, 50, 75, 100),
                           # breaks = 5,
                           labels = c("0%", "1 - 25%", "26 - 50%", "51 -75%", "75 - 100%"),
                           include.lowest = TRUE))
  
ggplot(data = onset_grouped, aes(x = rain_to_date_groups)) +
  geom_bar(aes(fill = as.factor(breeding_status)))

ggplot(data = onset_grouped, aes(x = interpolated_canopy)) +
  geom_histogram()

ggplot(data = onset_grouped, aes(x = canopy_groups)) +
  geom_bar(aes(fill = as.factor(breeding_status)))

# frailty model with grouped rainfall data
cox_frailty_groups <- coxph(Surv(dayOfWY, next_survey, breeding_status) ~ 
                             rain_to_date_groups +
                             canopy_groups +
                             # water_flow +
                             # water_regime +
                             frailty(LocationID), 
                           data = onset_grouped, 
                           control = coxph.control(iter.max = 50),
                           x = TRUE)


summary(cox_frailty_groups, conf.int = 0.89)
extractAIC(cox_frailty_groups)

test_cox <- cox.zph(cox_frailty_groups, transform = "identity") # put desired model name here
ggcoxzph(test_cox)
print(test_cox)

palette_green <- c(
  "#A5DD69",
  "#87D237",
  "#68A626",
  "#49741A",
  "#2A430F")

palette_blue <- c(
  "#85C8FF",
  "#47ACFF",
  "#0A91FF",
  "#0070CC",
  "#004F8F")

palette_brown <- c(
  "#B6AFA4",
  "#9B9182",
  "#7D7364",
  "#5B5349",
  "#FFB238")


predict_fun <- function(...) {
  1 - predictRisk(...)
}

adjusted_curves <- adjustedsurv(
  data = onset_grouped,
  variable = "canopy_groups",
  ev_time = "dayOfWY",
  event = "breeding_status",
  method = "direct",
  # bootstrap = TRUE,
  conf_level = 0.89,
  outcome_model = cox_frailty_groups,
  predict_fun = predict_fun,
)

plot(adjusted_curves, 
     use_boot = TRUE,
     cif = TRUE,
     xlab = "Day of Water Year",
     ylab = "Cumulative Incidence of Breeding",
     custom_colors = palette_blue) +
  theme_bw() +
  labs(color = "Rainfall Group")

plot(adjusted_curves, 
     use_boot = TRUE,
     cif = TRUE,
     # custom_linetypes = c(3,3,3,3,1),
     xlab = "Day of Water Year",
     ylab = "Cumulative Incidence of Breeding",
     legend.title = NULL,
     legend.position = NULL,
     custom_colors = palette_brown) +
  theme_bw() +
  labs(color = "Canopy Cover")
  # scale_linetype_manual(values = c(3,3,3,3,1), name = NULL, labels = NULL)

# ungrouped rain_to_date model
rain_model <- coxph(Surv(dayOfWY, next_survey, breeding_status) ~ 
                             rain_to_date +
                             # interpolated_canopy +
                             frailty(LocationID), 
                           data = complete_onset, 
                           x = TRUE)

# alt models with rain to date or sun hours
sun_model <- coxph(Surv(dayOfWY, next_survey, breeding_status) ~ 
              cum_sun_hours +
              frailty(LocationID), 
            data = complete_onset, 
            x = TRUE)

days_model <- coxph(Surv(dayOfWY, next_survey, breeding_status) ~ 
              days_since_first_rain +
              frailty(LocationID), 
            data = complete_onset, 
            x = TRUE)

Cand.mods <- list("Cumulative Rainfall Model" = rain_model, "Sun Hours Model" = sun_model, "Days Since First Rain Model" = days_model)

aictab(cand.set = Cand.mods)

test_cox <- cox.zph(cox_model_frailty, transform = "identity") # put desired model name here
ggcoxzph(test_cox)
print(test_cox)

# forest plot
ggforest(cox_model_frailty, data = complete_onset)

new_data_rain <- expand.grid(
  rain_to_date = seq(min(complete_onset$rain_to_date, na.rm = TRUE),
                     max(complete_onset$rain_to_date, na.rm = TRUE),
                     length.out = 162))

new_data_rain$hazard_ratio <- exp(predict(cox_model_frailty, newdata = new_data_rain, type = "risk"))

main_color <- "#49741A"
background <- "#7DC82D"
background2 <- "#91D548"

# hazard ratio for ungrouped rain to date
rain_hazard_plot <- ggplot(new_data_rain, aes(x = rain_to_date, y = hazard_ratio)) +
  geom_line(
  color = main_color, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = background) +
  labs(x = "Cumulative Rainfall (cm)", y = "Hazard Ratio") +
  theme_bw()
rain_hazard_plot

ggforest(cox_model_frailty, data = complete_onset)

# testing assumptions
# to use the cox model, results of test_assumptions must not be significant
test_cox <- cox.zph(cox_frailty_groups, transform = "identity") # put desired model name here
ggcoxzph(test_cox)
print(test_cox)

#### coxme + plots -- no longer using ####
# coxme with random effects: rain to date + cumulative sun hours
coxme_model <- coxme(Surv(dayOfWY, next_survey, breeding_status) ~ 
                     rain_to_date_groups +
                     # sun_resid +
                     # cum_sun_hours +
                     (1 | LocationID),
                   data = onset_grouped)

summary(coxme_model)

new_data <- expand_grid(rain_to_date_groups = unique(onset_grouped$rain_to_date_groups), 
              LocationID = unique(onset_grouped$LocationID),
              dayOfWY = seq(min(onset_grouped$dayOfWY), max(onset_grouped$dayOfWY), length.out = 100))

new_data$hazard_ratio <- predict_coxme(coxme_model, newdata = new_data, type = "risk")

km_fit <- survfit(Surv(dayOfWY, next_survey, breeding_status) ~ LocationID, 
                  data = onset_grouped)


baseline_hazard <- coxph(Surv(dayOfWY, time2, breeding_status) ~ 1,
                    data = complete_onset)


#### maually plotting coxme survival curves -- keeping for reference ####
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

coefficients <- as.data.frame(summary(cox_frailty_groups)$coefficients) %>% 
  rename(
    `se` = `se(coef)`
  ) %>% 
  mutate(
    hazard_ratio = exp(coef),
    lower_ci = exp(coef - 1.96*se),
    upper_ci = exp(coef + 1.96*se)
    # cov = as.factor(c("Low-medium rainfall", "medium rainfall", "medium-high rainfall", "high rainfall", "site effect"))
  ) %>%
  rownames_to_column(var = "cov") %>% 
  filter(cov != "site effect")

ggplot(coefficients, aes(x = hazard_ratio, y = cov)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0) +
  geom_point() +
  labs(x = "Hazard Ratio", y = "Covariate") +
  theme_minimal()

#### *** SURVIVAL MODELS -- no longer using, keep for reference *** ####

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
#### eda plots -- trying to figure out patterns in the data ####

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

ggplot(onset_of_breeding_surv, aes(x = WaterTemp, y = interpolated_canopy)) +
  geom_point() +
  geom_smooth(method = "glm")

cor(complete_onset$WaterTemp, complete_onset$interpolated_canopy)

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

