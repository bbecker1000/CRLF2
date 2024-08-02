#ulam / Stan model for joint probabilities and missing data

library(tidyverse)
library(plyr)
library(rethinking)
library(rstan) #rstan after rethinking to ensure correct namespace for "stanfit" objects
library(here)

#Use full dataset with NAs
between_year_data <- read_csv(here::here("data", "between_year_data.csv"))
hist(between_year_data$BRDYEAR)




#setup variables
EggMasses = between_year_data$num_egg_masses 
Year = (between_year_data$BRDYEAR-2009)
Year_2 = ((between_year_data$BRDYEAR-2009)^2)
  # mean_percent_emerg + 
  # mean_percent_sub +
Water = scale(between_year_data$mean_percent_water)
Canopy = scale(between_year_data$interpolated_canopy)
Rain = scale(between_year_data$yearly_rain)
Rain_2 = scale(between_year_data$yearly_rain^2)
mean_max_depth = scale(between_year_data$mean_max_depth)
max_depth = scale(between_year_data$max_depth)
AirTemp = scale(between_year_data$AirTemp)
AirTemp_2 = scale((between_year_data$AirTemp)^2)
WaterTemp = scale(between_year_data$WaterTemp)
WaterTemp_2 = scale((between_year_data$WaterTemp)^2)
CoastalSite = between_year_data$CoastalSite
Mean_salinity = scale(between_year_data$mean_salinity)
Watershed = between_year_data$Watershed #RE
Location = between_year_data$LocationID #RE

unique(Location)

#make data frame
dat=data.frame(
  EggMasses,
  Year,
  Year_2,
  # Water,
  Canopy, 
  Rain,
  Rain_2,
  mean_max_depth,
  #  max_depth,
  # AirTemp,
  #  AirTemp_2,
  WaterTemp,  #leave NAs
  #  WaterTemp_2, 
  CoastalSite,
  #  Mean_salinity,
  #  Watershed,
  Location
)

#if including canopy in model, filter NAs
nrow(dat)
dat <- filter(dat, !is.na(Canopy))
nrow(dat)
unique(dat$Location)




#all data in model = 43 sites
#stan doesn't take strings for factors, must be >0
# Location <- as.integer(revalue(Location,
# c(
#   "KC01"="1", "KC02"="2", "KC03"="3", "KC04"="4", 
#   "LS01"="5", "LS04"="6", "LS05"="7", "LS06"="8", "LS07"="9", "LS08"="10", "LS09"="11","LS11"="12",
#   "MC01"="13", 
#   "RC01"="14", "RC02"="15", "RC03"="16", "RC05"="17", "RC07"="18", "RC08"="19",
#   "RC09"="20", "RC10"="21", "RC11"="22", "RC12"="23", "RC13"="24", "RC14"="25",
#   "RC15"="26", "RC16"="27", "RC17"="28", "RC18"="29", "RC19"="30", "RC20"="31",
#   "RC23"="32", "RC24"="33", "RC25"="34", "RC26"="35", 
#   "RL02"="36", "RL04"="37", "RL06"="38", "RL07"="39", 
#   "TV02"="40", "TV03"="41", "TV06"="42", "WG01"="43"
#   )
# ))

#for canopy in model = 26 sites
dat$Location <- as.integer(revalue(dat$Location,
                               c(
"KC01"=1, "KC02"=2, "KC03"=3, "LS01"=4, "LS05"=5, "LS06"=6, "LS07"=7, "LS08"=8, "LS09"=9, 
"MC01"=10, 
"RC07"=11,"RC10"=12, "RC11"=13, "RC13"=14, "RC14"=15, "RC15"=16, "RC17"=17, "RC20"=18, "RC24"=19, 
"RC25"=20, "RL02"=21, "RL04"=22,
"TV02"=23, "TV03"=24, "TV06"=25, "WG01"=26
)
))

#same for watershed
unique(dat$Watershed)
Watershed <- as.integer(revalue(
  Watershed,
  c(
    "Kanoff Creek" = "1",
    "Laguna Salada" = "2",
    "Milagra Creek" = "3",
    "Redwood Creek" = "4",
    "Rodeo Lagoon" = "5",
    "Tennessee Valley" = "6",
    "Wilkins Gulch" = "7"
  )
))

#same for CoastalSite
dat$CoastalSite <- as.integer(ifelse(dat$CoastalSite == "TRUE", 1, 2))





#run model with missing data and joint probabilities
t0 <- Sys.time()

m1.ulam.cmdstan.F <-  ulam(
  alist(
    #Egg model
    EggMasses ~ dgampois( mu, exp( phi )),
    log(mu) <- 
      a_Location[Location] + #random slope and random intercept 
      beta_Year*Year +
      beta_Year_2+Year_2 +
      #beta_Water*Water,
      beta_Rain*Rain +
      beta_Rain_2*Rain_2 +
      beta_mean_max_depth*mean_max_depth +
      #max_depth,
    # beta_AirTemp*AirTemp +
      #AirTemp_2,
      beta_WaterTemp*WaterTemp +
      #WaterTemp_2, 
      beta_CoastalSite*CoastalSite,
      #Mean_salinity,
      #Watershed,
      
    #Location RE model
    a_Location[Location] ~ dnorm(mu_Location, tau_Location),
    mu_Location ~ dnorm(0, 0.5),
    tau_Location ~ dexp(1),
    
    #Temp model
    WaterTemp ~ normal(WaterTemp_nu , WaterTemp_tau),
    WaterTemp_nu <-
      a_WaterTemp +
      beta_Canopy*Canopy,
    
    #mean_max_depth model
    mean_max_depth ~ normal(mean_max_depth_nu, mean_max_depth_tau),
    mean_max_depth_nu <-
      a_mean_max_depth +
      beta_Rain*Rain + 
      beta_Rain_2*Rain_2,
    
    #fixed effects priors
    c(a_Location, 
      a_WaterTemp, 
      a_mean_max_depth, 
      beta_Year,
      beta_Year_2,
      beta_mean_max_depth,
      beta_CoastalSite,
      beta_Rain_2
  )  
  ~ normal( 0 , 0.5 ),      # uninformed
  
      beta_Rain             ~ normal( 0.25 , 0.25 ),  # more rain = more eggs
      beta_WaterTemp        ~ normal( 0.25 , 0.25 ),  # warmer = more eggs
      beta_Canopy           ~ normal( -0.25 , 0.25 ),  # less canopy = more eggs
      WaterTemp_tau ~ exponential(1),
      mean_max_depth_tau ~ exponential(1),
      phi ~ dnorm(1,5)
  ), 
  data=dat , chains=3 , cores=parallel::detectCores(), iter=4000,
  cmdstan=FALSE # FALSE to get stanfit object
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime

m1.ulam <- m1.ulam.cmdstan.F
precis(m1.ulam)
plot(m1.ulam,
     pars = c("beta_Year", 
              "beta_Year_2",
              "beta_Rain",
              "beta_Rain_2",
              "beta_WaterTemp",
              "beta_mean_max_depth",
              "beta_Canopy",
              "beta_CoastalSite"),
     xlab = "Beta Coefficient", 
     main = "Stan Model")
plot(m1.ulam, depth=2)

save(m1.ulam, file = "Code/Output/m1.ulam.RData")
#load("C:/bbecker/Projects/CRLF/code/Output/m1.ulam.RData")


rethinking::stancode(m1.ulam)

fit <- m1.ulam

plot(fit)

trankplot(fit@stanfit)
traceplot(fit@stanfit)

pairs(fit@stanfit)
par(mfrow=c(1,1))


###-----------------
## same model no missing data:

#Make another file that removes all NAs
dat_complete <- na.omit(dat)
dat_complete
nrow(dat_complete)
unique(dat_complete$Location)

t0 <- Sys.time()

m1.ulam.complete.case <-  ulam(
  alist(
    #Egg model
    EggMasses ~ dgampois( mu, exp( phi )),
    log(mu) <- 
      a_Location[Location] + #random slope and random intercept 
      beta_Year*Year +
      beta_Year_2+Year_2 +
      #beta_Water*Water,
      beta_Rain*Rain +
      beta_Rain_2*Rain_2 +
      beta_mean_max_depth*mean_max_depth +
      #max_depth,
      # beta_AirTemp*AirTemp +
      #AirTemp_2,
      beta_WaterTemp*WaterTemp +
      #WaterTemp_2, 
      beta_CoastalSite*CoastalSite,
    #Mean_salinity,
    #Watershed,
    
    #Location RE model
    a_Location[Location] ~ dnorm(mu_Location, tau_Location),
    mu_Location ~ dnorm(0, 0.5),
    tau_Location ~ dexp(1),
    
    #Temp model
    WaterTemp ~ normal(WaterTemp_nu , WaterTemp_tau),
    WaterTemp_nu <-
      a_WaterTemp +
      beta_Canopy*Canopy,
    
    #mean_max_depth model
    mean_max_depth ~ normal(mean_max_depth_nu, mean_max_depth_tau),
    mean_max_depth_nu <-
      a_mean_max_depth +
      beta_Rain*Rain,
    
    #fixed effects priors
    c(a_Location, 
      a_WaterTemp, 
      a_mean_max_depth, 
      beta_Year,
      beta_Year_2,
      beta_mean_max_depth,
      beta_CoastalSite,
      beta_Rain_2
    )  
    ~ normal( 0 , 0.5 ),      # uninformed
    
    beta_Rain             ~ normal( 0.25 , 0.25 ),  # more rain = more eggs
    beta_WaterTemp        ~ normal( 0.25 , 0.25 ),  # warmer = more eggs
    beta_Canopy           ~ normal( -0.25 , 0.25 ),  # less canopy = more eggs
    WaterTemp_tau ~ exponential(1),
    mean_max_depth_tau ~ exponential(1),
    phi ~ dnorm(1,5)
  ), 
  data=dat_complete , chains=3 , cores=parallel::detectCores(), iter=4000,
  cmdstan=FALSE # FALSE to get stanfit object
)

beepr::beep(0)
t1 <- Sys.time()
runtime <- t1-t0
runtime

fit <- m1.ulam.complete.case

save(fit, file = "Code/Output/m1.ulam.complete.case.RData")
#load("C:/bbecker/Projects/CRLF/code/Output/m1.ulam.complete.case.RData")
fit


precis(fit)
plot(fit,
     pars = c("beta_Year", 
              "beta_Year_2",
              "beta_Rain",
              "beta_Rain_2",
              "beta_WaterTemp",
              "beta_mean_max_depth",
              "beta_Canopy",
              "beta_CoastalSite"),
     xlab = "Beta Coefficient", 
     main = "Stan Model")
plot(fit, depth=2)
#dev.off()


rethinking::stancode(fit)



plot(fit)

trankplot(fit)
traceplot(fit@stanfit)

pairs(fit@stanfit)
par(mfrow=c(1,1))




