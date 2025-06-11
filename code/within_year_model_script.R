library(tidyverse)
library(survival)
library(survminer)
library(here)
library(ggplot2)
library(cowplot)
library(adjustedCurves)
library(riskRegression)
library(AICcmodavg)

setwd(here::here("code"))
complete_onset <- read_csv(here::here("data", "complete_onset_of_breeding.csv"))

#### frailty model + plots ####

# group rain_to_date and canopy cover into groups for survival model
onset_grouped <- complete_onset %>% 
  mutate(
    rain_to_date_groups = cut(rain_to_date, 
                              breaks = c(0, 25, 50, 75, 100, 125),
                              labels = c("0 - 25 cm", "25 - 50 cm", "50 - 75 cm", "75 - 100 cm", "100 - 125 cm"),
                              include.lowest = TRUE),
    canopy_groups = cut(interpolated_canopy,
                        breaks = c(0, 1, 25, 50, 75, 100),
                        labels = c("0%", "1 - 25%", "26 - 50%", "51 -75%", "75 - 100%"),
                        include.lowest = TRUE))


# frailty model
cox_frailty_groups <- coxph(Surv(dayOfWY, next_survey, breeding_status) ~ 
                              rain_to_date_groups +
                              canopy_groups +
                              frailty(LocationID), 
                            data = onset_grouped, 
                            control = coxph.control(iter.max = 50),
                            x = TRUE)


summary(cox_frailty_groups, conf.int = 0.89)

test_cox <- cox.zph(cox_frailty_groups, transform = "identity")
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

# rainfall plot
adjusted_curves <- adjustedsurv(
  data = onset_grouped,
  variable = "rain_to_date_groups",
  ev_time = "dayOfWY",
  event = "breeding_status",
  method = "direct",
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


# canopy cover plot
adjusted_curves <- adjustedsurv(
  data = onset_grouped,
  variable = "canopy_groups",
  ev_time = "dayOfWY",
  event = "breeding_status",
  method = "direct",
  conf_level = 0.89,
  outcome_model = cox_frailty_groups,
  predict_fun = predict_fun,
)

plot(adjusted_curves, 
     use_boot = TRUE,
     cif = TRUE,
     xlab = "Day of Water Year",
     ylab = "Cumulative Incidence of Breeding",
     legend.title = NULL,
     legend.position = NULL,
     custom_colors = palette_brown) +
  theme_bw() +
  labs(color = "Canopy Cover")


# alternate models with rain to date or sun hours. Used AICc to compare
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

Cand.mods <- list("Cumulative Rainfall Model" = cox_frailty_groups, "Sun Hours Model" = sun_model, "Days Since First Rain Model" = days_model)

aictab(cand.set = Cand.mods)
