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
    max_depth_scaled = as.vector(scale(MaxD_proportion)),
    AirTemp_scaled = as.vector(scale(AirTemp)),
    WaterTemp_scaled = as.vector(scale(WaterTemp)), 
    water_flow = as.factor(water_flow),
    water_regime = as.factor(water_regime), 
    Watershed = as.factor(Watershed),
    LocationID = as.factor(LocationID),
    cum_sun_hours_scaled = as.vector(scale(cum_sun_hours)),
    dir_dur_scaled = as.vector(scale(dir_dur))) %>% 
  select(-MaxD, -MaxD_yearly, -MaxD_proportion, -NumberofEggMasses)

# creating a "complete case" column
scaled_within_year$complete_case <- complete.cases(scaled_within_year)
complete_onset <- scaled_within_year %>% filter(complete_case == TRUE) %>% select(-complete_case)

#### *** SURVIVAL MODELS *** ####

#assign "dead" to all known breeders.  no censoring.
complete_onset$status <- 2 

#intercept model of the mean
fit.null <- survfit(Surv(first_breeding, status) ~ 1, data = onset_of_breeding_surv)
#survival (breeding probability) curves by watershed
fit.watershed <- survfit(Surv(rain_to_date, status) ~ Watershed + MaxD_proportion, data = onset_of_breeding_surv)

fit.rain <- survfit(Surv(first_breeding, status) ~ 1, data = onset_of_breeding_surv)
fit.watertemp <- survfit(Surv(WaterTemp, status) ~ 1, data = onset_of_breeding_surv)

#pick one to inspect/plot
fit <- fit.null
fit <- fit.watershed
fit <- fit.rain
fit <- fit.watershed.rw
fit <- fit.watertemp

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

# cox model: univariate (rain to date)
cox_model <- coxph(Surv(first_breeding, status) ~ rain_to_date, data = complete_onset)
summary(cox_model)

# cox model: univariate (rain to date)
cox_model <- coxph(Surv(first_breeding, status) ~ rain_to_date + cum_sun_hours, data = complete_onset)
summary(cox_model)

surv_fit <- survfit(cox_model)
plot(surv_fit)

# cox model: univariate (Site)
Site.cox <- coxph(Surv(rain_to_date, status) ~ LocationID, data = complete_onset)
summary(Site.cox)

# univariate cox models for continuous variables
covariates <- c("MaxD_proportion", "AirTemp", "WaterTemp", "BRDYEAR", "rain_to_date")
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(first_breeding, status) ~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = complete_onset)})

univ_results <- lapply(univ_models, function(x) {
  if (!is.null(x)) {
    x <- summary(x)
    p.value <- signif(x$wald["pvalue"], digits = 2)
    wald.test <- signif(x$wald["test"], digits = 2)
    beta <- signif(x$coef[1], digits = 2)
    HR <- signif(x$coef[2], digits = 2)
    HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
    res <- c(beta, HR, wald.test, p.value)
    names(res) <- c("beta", "HR (95% CI for HR)", "wald.test", "p.value")
    return(res)
  } else {
    return(rep(NA, 4))  # Return NA for missing models
  }
})
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

# multivariate with random effects -- scaled variables
multi.cox <- coxme(Surv(first_breeding, status) ~ 
                     scale(MaxD_proportion) + 
                     scale(AirTemp) + 
                     scale(WaterTemp) + 
                     scale(BRDYEAR) + 
                     scale(rain_to_date) + 
                     water_flow +
                     water_regime +
                     (1 | Watershed/LocationID),
                   data = complete_onset)

# see summary of the model:
summary(multi.cox)

#### plot model -- forest plot (this one works!) ####
coefficients <- as.data.frame(summary(multi.cox)$coefficients) %>% 
  rename(
    `estimate` = `exp(coef)`,
    `se` = `se(coef)`
  )
plot_data <- data.frame(
  covariate = rownames(coefficients),
  hazard_ratio = coefficients$estimate,
  lower_CI = coefficients$estimate - (1.96 * coefficients$se),
  upper_CI = coefficients$estimate + (1.96 * coefficients$se)
)

ggplot(plot_data, aes(x = hazard_ratio, y = covariate)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0) +
  geom_point() +
  labs(x = "Hazard Ratio", y = "Covariate", title = "Coefficient Estimates") +
  theme_minimal()

#### testing assumptions for survival models ####
# to use the cox model, results of test_assumptions must not be significant
test_assumptions <- cox.zph(multi.cox)
test_assumptions
print(test_assumptions)
plot(test_assumptions)
