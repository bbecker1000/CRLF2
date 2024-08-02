##DAGs for Frog models

library(dagitty)
#### most up to date model ####
DAG_FROG <- dagitty(
  "dag{
  Year -> FROG ;
  Rain -> FROG ;
  Site -> FROG ;
  Rain -> Vegetation -> FROG ;
  Rain -> Depth -> FROG ;
  Rain -> Site -> FROG;
  Salinity -> FROG ;
  Rain -> Salinity ;
  Canopy -> Temp -> FROG ;
  Temp -> FROG ;
  Water_Regime -> Rain ;
  Site -> Water_Flow -> FROG ;


  Year [exposure] ;
  Site [exposure] ;
  Rain [exposure] ;
  Vegetation [exposure] ;
  Salinity [exposure] ;
  Canopy [exposure] ;
  Temp [exposure] ;
  Water_Regime [exposure] ;
  Water_Flow [exposure] ;
  FROG [outcome] ;
}"
)

plot(DAG_FROG)

## add coordinates
coordinates(DAG_FROG) <- list(
  x = c(
    Year = 2,
    Rain = 2.5,
    Site = 3,
    Vegetation = 2,
    Depth = 2.6,
    Water_Flow = 3,
    Salinity = 2.75,
    FROG = 2.5,
    Canopy = 2.25,
    Temp = 2.25,
    Water_Regime = 2.5
  ),
  y = c(
    Year = -1,
    Rain = -4,
    Site = -4,
    Depth = -2,
    Vegetation = -4,
    Water_Flow = -2,
    Salinity = -2,
    FROG = 2,
    Canopy = -3,
    Temp = -2,
    Water_Regime = -5
  )
)
plot(DAG_FROG)

##### frog model using ggdag ####
library(ggdag)

coord_dag <- list(
  FROG = c(d = 0, x = 1, y = 4),
  Year = c(d = 0, x = 1, y = 0)
)

FROG_ggdag <- dagify(
  FROG ~ Year,
  FROG ~ Rain,
  FROG ~ Site,
  Rain ~ Vegetation,
  Rain ~ Site,
  Depth ~ Rain,
  FROG ~ Depth,
  FROG ~ Bullfrog,
  FROG ~ Salinity,
  exposure = "Year",
  outcome = "FROG",
  coords = coord_dag
)
ggdag(FROG_ggdag, 
      text_col = "white",
      stylized = TRUE)
  theme_dag()
ggdag_paths(FROG_ggdag, 
            #adjust_for = c("Breach"),
            text_col = "black")

ggdag_adjustment_set(FROG_ggdag, text_col = "black")




#####----------------------------------
#try with ggdag
#still from another project.  can add frog data.
library(ggdag)



Phoca_dag <- dagify(
  POPULATION ~ Year,
  Prey ~ Oceanography,
  Pup ~ Prey,
  Adult ~ Prey,
  POPULATION ~ Pup,
  POPULATION ~ Adult,
  POPULATION ~ Molt,
  Pup ~ Coyote,
  Adult ~ Coyote,
  Coyote ~ Site,
  
  exposure = "Coyote",#"Oceanography",
  outcome = "POPULATION"
  
)
ggdag(Phoca_dag, text_col = "white", 
      stylized = TRUE) + 
  theme_dag_blank()
ggdag_paths(Phoca_dag, 
            #adjust_for = c("Breach"),
            text_col = "black")

ggdag_adjustment_set(Phoca_dag, text_col = "black")

library(dagitty)

adjustmentSets(goby_dag)

