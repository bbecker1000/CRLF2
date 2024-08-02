### if anyone is having the same glmmTMB issues as I was, running this script should help!!
### I think with some modifications, this might be able to fix lme4 too? 
### I don't have the lme4 issue, so let me know if the lme4 fix works :)

# install groundhog package
while (!require("groundhog", quietly=TRUE)) install.packages("groundhog")

## retrieve build date of installed version of TMB
bd1 <- as.character(asDateBuilt(packageDescription("TMB",fields="Built")))

#install packages from that earlier date
groundhog.library("TMB", bd1)
groundhog.library("Matrix", bd1)

# retrieve build date of lme4
bd2 <- as.character(asDateBuilt(packageDescription("lme4", fields = "Built")))

#install lme4 from earlier date
groundhog.library("lme4", bd2)
