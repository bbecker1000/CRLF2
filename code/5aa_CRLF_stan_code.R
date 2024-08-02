#Frog Stan

stan_program <- "

functions{
    vector merge_missing( array[] int miss_indexes , vector x_obs , vector x_miss ) {
        int N = dims(x_obs)[1];
        int N_miss = dims(x_miss)[1];
        vector[N] merged;
        merged = x_obs;
        for ( i in 1:N_miss )
            merged[ miss_indexes[i] ] = x_miss[i];
        return merged;
    }
}
data{
    array[236] int EggMasses;
    array[236] int CoastalSite;
     vector[236] Year_2;
     vector[236] Year;
    array[236] int Location;
     vector[236] WaterTemp;
    array[18] int WaterTemp_missidx;
     vector[236] Canopy;
     vector[236] mean_max_depth;
    array[22] int mean_max_depth_missidx;
     vector[236] Rain;
}
parameters{
     real mu_Location;
     real<lower=0> tau_Location;
     real beta_CoastalSite;
     real beta_mean_max_depth;
     real beta_Year_2;
     real beta_Year;
     real a_mean_max_depth;
     real a_WaterTemp;
     vector[26] a_Location;
     real beta_Rain;
     real beta_WaterTemp;
     real beta_Canopy;
     real<lower=0> WaterTemp_tau;
     real<lower=0> mean_max_depth_tau;
     real phi;
     vector[18] WaterTemp_impute;
     vector[22] mean_max_depth_impute;
}
model{
     vector[236] mu;
     vector[236] WaterTemp_merge;
     vector[236] WaterTemp_nu;
     vector[236] mean_max_depth_merge;
     vector[236] mean_max_depth_nu;
    phi ~ normal( 1 , 5 );
    mean_max_depth_tau ~ exponential( 1 );
    WaterTemp_tau ~ exponential( 1 );
    beta_Canopy ~ normal( -0.25 , 0.25 );
    beta_WaterTemp ~ normal( 0.25 , 0.25 );
    beta_Rain ~ normal( 0.25 , 0.25 );
    a_Location ~ normal( 0 , 0.5 );
    a_WaterTemp ~ normal( 0 , 0.5 );
    a_mean_max_depth ~ normal( 0 , 0.5 );
    beta_Year ~ normal( 0 , 0.5 );
    beta_Year_2 ~ normal( 0 , 0.5 );
    beta_mean_max_depth ~ normal( 0 , 0.5 );
    beta_CoastalSite ~ normal( 0 , 0.5 );
    for ( i in 1:236 ) {
        mean_max_depth_nu[i] = a_mean_max_depth + beta_Rain * Rain[i];
    }
    mean_max_depth_merge = merge_missing(mean_max_depth_missidx, to_vector(mean_max_depth), mean_max_depth_impute);
    mean_max_depth_merge ~ normal( mean_max_depth_nu , mean_max_depth_tau );
    for ( i in 1:236 ) {
        WaterTemp_nu[i] = a_WaterTemp + beta_Canopy * Canopy[i];
    }
    WaterTemp_merge = merge_missing(WaterTemp_missidx, to_vector(WaterTemp), WaterTemp_impute);
    WaterTemp_merge ~ normal( WaterTemp_nu , WaterTemp_tau );
    tau_Location ~ exponential( 1 );
    mu_Location ~ normal( 0 , 0.5 );
    a_Location ~ normal( mu_Location , tau_Location );
    for ( i in 1:236 ) {
        mu[i] = a_Location[Location[i]] + beta_Year * Year[i] + beta_Year_2 + Year_2[i] + beta_Rain * Rain[i] + beta_mean_max_depth * mean_max_depth_merge[i] + beta_WaterTemp * WaterTemp_merge[i] + beta_CoastalSite * CoastalSite[i];
        mu[i] = exp(mu[i]);
    }
    EggMasses ~ neg_binomial_2( mu , exp(phi) );
}
"


library(rstan)
fit1 <- stan(
  model_code = stan_program,  # Stan program
  data = dat.list,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 250,          # number of warmup iterations per chain
  iter = 500,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  refresh = 1            # no progress shown
)
