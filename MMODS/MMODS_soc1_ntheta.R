####
## Parameters
####

## NOTE: This version fits to continuous mobility data. The data we are using is private. Please email us or use your own form of contious mobility data.

set.seed(10001)               ## Set seed to recreate fits
fitting           <- TRUE     ## Small change in pomp objects if fitting or simulating
fit.minus         <- 0        ## Use data until X days prior to the present
more.params.uncer <- FALSE    ## Fit with more (FALSE) or fewer (TRUE) point estimates for a number of parameters
fit.E0            <- TRUE     ## Also fit initial number of infected individuals that starts the epidemic?
usable.cores      <- 6        ## Number of cores to use to fit
fit.with          <- "D_C"    ## Fit with D (deaths) or H (hospitalizations) -- Need your own data for H -- or both deaths and cases (D_C).
fit_to_sip        <- TRUE     ## Fit beta0 and shelter in place strength simultaneously?
meas.nb           <- TRUE     ## Negative binomial measurement process?
import_cases      <- FALSE    ## Use importation of cases?
## mif2 fitting parameters. 
n.mif_runs        <- 6        ## number of repeated fits (6 used in manuscript, 2 suggested to debug/check code)
n.mif_length      <- 120      ## number of steps (100 used in manuscript, 20 suggested to debug/check code)
n.mif_particles   <- 300      ## number of particles (3000 used in manuscript, 3000 suggested to debug/check code)
n.mif_rw.sd       <- 0.02     ## particle perturbation (0.02 used in manuscript, 0.02 suggested to debug/check code)
nparams           <- 200      ## number of parameter sobol samples (200 used in manuscript, 5 suggested to debug/check code)
nsim              <- 200      ## number of stochastic epidemic simulations for each fitted beta0 for dynamics (300 used in manuscript, 150 suggested to debug/check code)
sim.length        <- 300

## Extra params to explore for mutli data stream fitting process
detect.logis      <- TRUE
fixed.E0          <- FALSE

print("MMODS")
focal.county      <- "NA" 
focal.state_abbr  <- "NA"

## Required packages to run this code
needed_packages <- c(
    "pomp"
  , "plyr"
  , "dplyr"
  , "ggplot2"
  , "magrittr"
  , "scales"
  , "lubridate"
  , "tidyr"
  , "foreach"
  , "doParallel"
  , "data.table"
  , "doRNG")

## load packages. Install all packages that return "FALSE"
lapply(needed_packages, require, character.only = TRUE)

## Be very careful here, adjust according to your machine's capabilities
# registerDoParallel(cores = usable.cores)
 registerDoParallel(cores = (Sys.getenv("SLURM_NTASKS_PER_NODE")))
  
source("COVID_pomp_gammabeta.R")

county.data      <- read.csv("case_data.csv")
county.data$date <- as.Date(county.data$date)
mobility         <- read.csv("mobility_data.csv"); names(mobility)[1] <- "date"
mobility$date    <- as.Date(mobility$date)
mobility         <- mobility %>% mutate(sip_prop = workplaces * -1)

params <- read.csv("../params.csv", stringsAsFactors = FALSE)
params <- params %>% mutate(Value = sapply(est, function(x) eval(parse(text = x))))

fixed_params        <- params$Value
names(fixed_params) <- params$Parameter
if (!import_cases) {fixed_params["import_rate"] <- 0}
fixed_params        <- c(fixed_params, N = 100000)

source("variable_params_less_cont_testrun.R")

## Run parameters
sim_start  <- variable_params$sim_start
sim_end    <- sim_start + sim.length

param_array <- array(
  data = 0
, dim  = c(nparams, n.mif_runs, 10))
dimnames(param_array)[[3]] <- c(
  "log_lik"
, "log_lik_es"
, "beta0est"
, "E_init"
, {
  if (detect.logis) {
 c("detect_k", "detect_mid") 
  } else {
 c("detect_t0", "detect_t1") 
  }
}
, c("detect_max"
, "theta"
, "theta2"
, "beta_min"
)
  )  

startvals <- array(
  data = 0
, dim  = c(8, n.mif_runs, nparams))
dimnames(startvals)[[1]] <- c(
  "beta0"
, "E_init"
, {
  if (detect.logis) {
 c("detect_k", "detect_mid") 
  } else {
 c("detect_t0", "detect_t1") 
  }
}
, c("detect_max"
, "theta"
, "theta2"
, "beta_min"  
)
)

Reff <- matrix(data = 0, nrow = nrow(variable_params), ncol = sim.length)

for (i in 1:nrow(variable_params)) {
  
  ## Adjust the data for the current start date
county.data <- county.data %>% 
  mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) %>% 
  dplyr::filter(day > 0) %>%
  dplyr::filter(!is.na(deaths), !is.na(cases))

county.data <- rbind(
  data.frame(
    day    = seq(1:(min(county.data$day) - 1))
  , date   = as.Date(seq(1:(min(county.data$day) - 1)), origin = variable_params[i, "sim_start"])
  , deaths = NA
  , cases  = NA
  )
, county.data 
  )

 mobility <- mobility %>% 
  # three day moving median
#  mutate(sip_prop = mutate(., 
#                         sip_prop_lag = lag(sip_prop),
#                         sip_prop_lead = lead(sip_prop)) %>% 
#           select(contains("sip_prop")) %>%
#           purrr::pmap_dbl(~median(c(...)))) %>% 
  mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) %>% 
  dplyr::filter(date >= variable_params[i, ]$sim_start) %>%
  dplyr::filter(!is.na(sip_prop)) %>% 
  dplyr::filter(day > 0) %>% 
  dplyr::select(date, sip_prop, day)

if (min(mobility$day) > 1) {
mobility <- rbind(
  data.frame(
    date     = as.Date(seq(1:(min(mobility$day) - 1)), origin = variable_params[i, ]$sim_start)
  , sip_prop = mean(mobility$sip_prop[1:5])
  , day      = seq(1:(min(mobility$day) - 1))
  )
, mobility
  )  
}
 
## Remove dates after one week after movement data ends
county.data <- county.data %>% dplyr::filter(date < (max(mobility$date) + 7))

mobility <- mobility %>% mutate(sip_prop = plogis(scale(sip_prop)))

if (detect.logis) {
mob.covtab <- covariate_table(
    sip_prop     = mobility$sip_prop
 ,  order        = "constant"
 ,  times        = mobility$day
 ,  intervention  = rep(0, nrow(mobility))
 ,  detect_t0     = min(county.data[!is.na(county.data$cases), ]$day) - 1
    )
} else {
mob.covtab <- covariate_table(
    sip_prop     = mobility$sip_prop
 ,  order        = "constant"
 ,  times        = mobility$day
 ,  intervention  = rep(0, nrow(mobility))
 ,  detect_t0     = min(county.data[!is.na(county.data$cases), ]$day) - 1
    )  
}

if (detect.logis) {
covid_mobility <- pomp(
   data       = county.data
 , times      = "day"
 , t0         = 1
 , covar      = mob.covtab
 , rprocess   = euler(sir_step_mobility, delta.t = 1/6)
 , rmeasure   = rmeas_multi_logis
 , dmeasure   = dmeas_multi_logis
 , rinit      = sir_init
 , partrans   = par_trans
 , accumvars  = accum_names
 , paramnames = param_names
 , statenames = state_names
)  
} else {
covid_mobility <- pomp(
   data       = county.data
 , times      = "day"
 , t0         = 1
 , covar      = mob.covtab
 , rprocess   = euler(sir_step_mobility, delta.t = 1/6)
 , rmeasure   = rmeas_multi_pwl
 , dmeasure   = dmeas_multi_pwl
 , rinit      = sir_init
 , partrans   = par_trans
 , accumvars  = accum_names
 , paramnames = param_names
 , statenames = state_names
)  
}

  
## Extra params for the new model:
 ## 1) beta0_sigma can't really be fit, so assume for now
 ## 2) fit with no intervention and then simulate later with an intervention
fixed_params["beta0_sigma"] <- 0.16
fixed_params["beta_catch"]  <- 1
fixed_params["beta_red"]    <- 1 

if (variable_params[i, ]$beta0est == 0) {

  ## Fit runs in parallel based on number of cores
if (detect.logis) {
checktime <- system.time({
## Run mifs_local one time as a "burn-in" run (using MCMC terms...)
registerDoRNG(610408799)
mifs_local <- foreach(j = 1:n.mif_runs, .combine = c) %dopar%  {
    
library(pomp)
library(dplyr)
  
start_vals <- c(
      beta0       = rlnorm(1, log(0.7), 0.3)
    , E_init      = rpois(1, 2) + 1
    , detect_k    = rlnorm(1, log(0.1), 0.2)
    , detect_mid  = rlnorm(1, log(60), 0.2)
    , detect_max   = runif(1, 0.1, 0.7)
    , theta        = rlnorm(1, log(1), 0.4)
    , theta2       = rlnorm(1, log(1), 0.4)
    , beta_min     = runif(1, 0.3, 0.7) 
)

mifs_temp <- covid_mobility %>%
  mif2(
    t0      = 1
  , params  = c(
    c(fixed_params
    , Ca    = variable_params[i, ]$Ca
    , alpha = variable_params[i, ]$alpha
    , delta = variable_params[i, ]$delta
    , mu    = variable_params[i, ]$mu
      )
  , {
## random start for each run
    c(
      start_vals["beta0"]
    , start_vals["E_init"]
    , start_vals["detect_k"]
    , start_vals["detect_mid"]
    , start_vals["detect_max"]
    , start_vals["theta"]
    , start_vals["theta2"]
    , start_vals["beta_min"]
      )
  }
  )
  , Np     = n.mif_particles
  , Nmif   = n.mif_length
  , cooling.fraction.50 = 0.75
  , rw.sd  = rw.sd(
      beta0       = 0.02
    , E_init      = ivp(0.02)
    , detect_max  = 0.02
    , detect_mid  = 0.02
    , detect_k    = 0.02
    , theta       = 0.02
    , theta2      = 0.02
    , beta_min    = 0.02
  )
        ) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.75) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.50) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.50) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.50) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.30) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.30) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.30) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.10) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.10) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.10)
  
# ll <- replicate(10, mifs_temp %>% pfilter(Np = 50000) %>% logLik())
ll <- replicate(10, mifs_temp %>% pfilter(Np = 5000) %>% logLik())
ll <- logmeanexp(ll, se = TRUE)
return(list(mifs_temp, ll, start_vals))

}
})
} else {
checktime <- system.time({
## Run mifs_local one time as a "burn-in" run (using MCMC terms...)
registerDoRNG(610408799)
mifs_local <- foreach(j = 1:n.mif_runs, .combine = c) %dopar%  {
    
library(pomp)
library(dplyr)
  
start_vals <- c(
      beta0        = rlnorm(1, log(0.7), 0.3)
    , E_init       = rpois(1, 2) + 1
#    , detect_t0    = rlnorm(1, log(25), 0.2)
    , detect_t1    = rlnorm(1, log(60), 0.2)
    , detect_max   = runif(1, 0.1, 0.7)
    , theta        = rlnorm(1, log(1), 0.4)
    , theta2       = rlnorm(1, log(1), 0.4)
    , beta_min     = runif(1, 0.3, 0.7) 
)

mifs_temp <- covid_mobility %>%
  mif2(
    t0      = 1
  , params  = c(
    c(fixed_params
    , Ca    = variable_params[i, ]$Ca
    , alpha = variable_params[i, ]$alpha
    , delta = variable_params[i, ]$delta
    , mu    = variable_params[i, ]$mu
      )
  , {
## random start for each run
    c(
      start_vals["beta0"]
    , start_vals["E_init"]
#    , start_vals["detect_t0"]
    , start_vals["detect_t1"]
    , start_vals["detect_max"]
    , start_vals["theta"]
    , start_vals["theta2"]
    , start_vals["beta_min"]
      )
  }
  )
  , Np     = n.mif_particles
  , Nmif   = n.mif_length
  , cooling.fraction.50 = 0.75
  , rw.sd  = rw.sd(
      beta0       = 0.02
    , E_init      = ivp(0.02)
    , detect_max  = 0.02
#    , detect_t0   = 0.02
    , detect_t1   = 0.02
    , theta       = 0.02
    , theta2      = 0.02
    , beta_min    = 0.02
  )
        ) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.75) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.50) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.50) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.50) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.30) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.30) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.30) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.10) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.10) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.10)
  
# ll <- replicate(10, mifs_temp %>% pfilter(Np = 50000) %>% logLik())
ll <- replicate(10, mifs_temp %>% pfilter(Np = 5000) %>% logLik())
ll <- logmeanexp(ll, se = TRUE)
return(list(mifs_temp, ll, start_vals))

}
})  
}

mifs.sv           <- mifs_local[seq(3, (n.mif_runs * 3), by = 3)]
startvals[, , i]  <- sapply(mifs.sv, c, simplify = "array")

mifs.ll    <- mifs_local[seq(2, (n.mif_runs * 3), by = 3)]

loglik.est <- numeric(n.mif_runs)
loglik.se  <- numeric(n.mif_runs)
for (j in seq_along(loglik.est)) {
loglik.est[j] <- mifs.ll[[j]][1]
loglik.se[j] <- mifs.ll[[j]][2]
}

mifs_local <- mifs_local[seq(1, (n.mif_runs * 3), by = 3)]
best.fit   <- which(loglik.est == max(loglik.est))

variable_params[i, "beta0est"]   <- coef(mifs_local[[best.fit]])["beta0"]
variable_params[i, "E_init"]     <- coef(mifs_local[[best.fit]])["E_init"]
variable_params[i, "detect_max"] <- coef(mifs_local[[best.fit]])["detect_max"]
if (detect.logis) {
variable_params[i, "detect_k"]  <- coef(mifs_local[[best.fit]])["detect_k"]
variable_params[i, "detect_mid"]  <- coef(mifs_local[[best.fit]])["detect_mid"]  
} else {
variable_params[i, "detect_t0"]  <- coef(mifs_local[[best.fit]])["detect_t0"]
variable_params[i, "detect_t1"]  <- coef(mifs_local[[best.fit]])["detect_t1"]  
}
variable_params[i, "theta"]      <- coef(mifs_local[[best.fit]])["theta"]
variable_params[i, "theta2"]     <- coef(mifs_local[[best.fit]])["theta2"]
variable_params[i, "beta_min"]   <- coef(mifs_local[[best.fit]])["beta_min"]

all.params <- sapply(mifs_local, coef, simplify = "array")

param_array[i, , "beta0est"]   <- all.params[which(dimnames(all.params)[[1]] == "beta0"), ]
param_array[i, , "E_init"]     <- all.params[which(dimnames(all.params)[[1]] == "E_init"), ]
if (detect.logis) {
param_array[i, , "detect_k"]  <- all.params[which(dimnames(all.params)[[1]] == "detect_k"), ]
param_array[i, , "detect_mid"]  <- all.params[which(dimnames(all.params)[[1]] == "detect_mid"), ]  
} else {
param_array[i, , "detect_t0"]  <- all.params[which(dimnames(all.params)[[1]] == "detect_t0"), ]
param_array[i, , "detect_t1"]  <- all.params[which(dimnames(all.params)[[1]] == "detect_t1"), ]  
}
param_array[i, , "detect_max"] <- all.params[which(dimnames(all.params)[[1]] == "detect_max"), ]
param_array[i, , "theta"]      <- all.params[which(dimnames(all.params)[[1]] == "theta"), ]
param_array[i, , "theta2"]     <- all.params[which(dimnames(all.params)[[1]] == "theta2"), ]
param_array[i, , "beta_min"]   <- all.params[which(dimnames(all.params)[[1]] == "beta_min"), ]

variable_params[i, "log_lik"]      <- loglik.est[best.fit]
variable_params[i, "log_lik_se"]   <- loglik.se[best.fit]
param_array[i, , "log_lik"]        <- loglik.est

}

## Simulate from fits
if (detect.logis) {

SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid_mobility
    , times        = mob.covtab@times
    , params = c(fixed_params
      , c(
        beta0      = variable_params[i, "beta0est"]
      , E_init     = variable_params[i, "E_init"]
      , detect_max = variable_params[i, "detect_max"]
      , detect_mid = variable_params[i, "detect_mid"]
      , detect_k   = variable_params[i, "detect_k"]
      , theta      = variable_params[i, "theta"]
      , theta2     = variable_params[i, "theta2"]
      , beta_min   = variable_params[i, "beta_min"]
      , Ca         = variable_params[i, "Ca"]
      , alpha      = variable_params[i, "alpha"]
      , delta      = variable_params[i, "delta"]
      , mu         = variable_params[i, "mu"]
      , E_init     = variable_params[i, "E_init"]
        ))
    , nsim         = nsim
    , format       = "d"
    , include.data = F
    , seed         = 1001)) %>% {
      rbind(.,
         group_by(., day) %>%
           dplyr::select(-.id) %>%
           summarise_all(median) %>%
                    mutate(.id = "median"))
    }

} else {
  
SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid_mobility
    , times        = mob.covtab@times
    , params = c(fixed_params
      , c(
        beta0      = variable_params[i, "beta0est"]
      , E_init     = variable_params[i, "E_init"]
      , detect_max = variable_params[i, "detect_max"]
      , detect_t0 = variable_params[i, "detect_t0"]
      , detect_t1   = variable_params[i, "detect_t1"]
      , theta      = variable_params[i, "theta"]
      , theta2     = variable_params[i, "theta2"]
      , beta_min   = variable_params[i, "beta_min"]
      , Ca         = variable_params[i, "Ca"]
      , alpha      = variable_params[i, "alpha"]
      , delta      = variable_params[i, "delta"]
      , mu         = variable_params[i, "mu"]
      , E_init     = variable_params[i, "E_init"]
        ))
    , nsim         = nsim
    , format       = "d"
    , include.data = F
    , seed         = 1001)) %>% {
      rbind(.,
         group_by(., day) %>%
           dplyr::select(-.id) %>%
           summarise_all(median) %>%
                    mutate(.id = "median"))
    }
  
}

SEIR.sim.s <- SEIR.sim %>% 
  dplyr::filter(.id != "median") %>%
  group_by(day) %>%
  summarize(S = mean(S), D = mean(D))

if (detect.logis) {
#betat      <- variable_params[i, "beta0est"] * (1 - variable_params[i, "beta_min"] * mob.covtab@table[2, ])  
betat     <- variable_params[i, "beta0est"] * exp(log(variable_params[i, "beta_min"])*mob.covtab@table[2, ])
} else {
#betat      <- variable_params[i, "beta0est"] * (1 - variable_params[i, "beta_min"] * mob.covtab@table[1, ])  
betat     <- variable_params[i, "beta0est"] * exp(log(variable_params[i, "beta_min"])*mob.covtab@table[2, ])
}

Reff.t       <- with(variable_params[i, ], covid_R0(
   beta0est      = betat
 , fixed_params  = c(fixed_params, unlist(variable_params[i, ]))
 , sd_strength   = 1
 , prop_S        = SEIR.sim.s$S / (fixed_params["N"] - SEIR.sim.s$D)
  )
  )
Reff[i, 1:length(Reff.t)]  <- Reff.t

## Save an Rds periodically
saveRDS(
   list(
    variable_params  = variable_params
  , fixed_params     = fixed_params
  , startvals        = startvals
  , Reff             = Reff
  , param_array      = param_array
   ), paste(
     paste("output/"
       , paste(focal.county, fit_to_sip, more.params.uncer, fit.minus, Sys.Date(), "cont", "temp_MMODS_ntheta", sep = "_")
         , sep = "")
     , "Rds", sep = "."))

}

## Save a final Rds
saveRDS(
   list(
    variable_params  = variable_params
  , fixed_params     = fixed_params
  , startvals        = startvals
  , Reff             = Reff
  , param_array      = param_array
   ), paste(
     paste("output/"
       , paste(focal.county, fit_to_sip, more.params.uncer, fit.minus, Sys.Date(), "cont", "final_MMODS_ntheta", sep = "_")
         , sep = "")
     , "Rds", sep = "."))
  