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
n.mif_runs        <- 2        ## number of repeated fits (6 used in manuscript, 2 suggested to debug/check code)
n.mif_length      <- 20       ## number of steps (100 used in manuscript, 20 suggested to debug/check code)
n.mif_particles   <- 3000     ## number of particles (3000 used in manuscript, 3000 suggested to debug/check code)
n.mif_rw.sd       <- 0.02     ## particle perturbation (0.02 used in manuscript, 0.02 suggested to debug/check code)
nparams           <- 5        ## number of parameter sobol samples (200 used in manuscript, 5 suggested to debug/check code)
nsim              <- 150      ## number of stochastic epidemic simulations for each fitted beta0 for dynamics (300 used in manuscript, 150 suggested to debug/check code)
sim.length        <- 200 

focal.county      <- "Miami-Dade" 
focal.state_abbr  <- "FL"

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
  , "data.table")

## load packages. Install all packages that return "FALSE"
lapply(needed_packages, require, character.only = TRUE)

## Be very careful here, adjust according to your machine's capabilities
  registerDoParallel(cores = usable.cores)
# registerDoParallel(cores = (Sys.getenv("SLURM_NTASKS_PER_NODE")))

source("COVID_pomp_gammabeta.R")

if (fit.with == "D_C" | fit.with == "D") {
## Scrape death data from the NYT github repo to stay up to date or load a previously saved dataset
#deaths <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
deaths  <- read.csv("us-counties.txt")
deaths  <- deaths %>% mutate(date = as.Date(date)) %>% filter(county == focal.county)
deaths  <- deaths %>% dplyr::filter(date < max(date) - fit.minus)
} else if (fit.with == "H") {
## Not supported right now for SCC, ony currently have data for Contra Costa County.
 ## To use H, supply your own data and change the path
hospit     <- read.csv("contra_costa/ccc_data.csv")
hospit     <- hospit %>% 
  mutate(date = as.Date(REPORT_DATE)) %>% 
  filter(CURRENT_HOSPITALIZED != "NULL") %>% 
  mutate(ch = as.numeric(as.character(CURRENT_HOSPITALIZED))) %>% 
  dplyr::select(date, ch)
hospit    <- hospit %>% dplyr::filter(date < max(date) - fit.minus)  
}

## Bring in parameters. Two sets used here depending on the desired level of uncertainty and number of parameters to sample over
if (!more.params.uncer) {
params <- read.csv("params.csv", stringsAsFactors = FALSE)
} else {
params <- read.csv("params_wide.csv", stringsAsFactors = FALSE) 
}
params <- params %>% mutate(Value = sapply(est, function(x) eval(parse(text = x))))

fixed_params        <- params$Value
names(fixed_params) <- params$Parameter
if (!import_cases) {fixed_params["import_rate"] <- 0}

## Location-specific parameters
location_params     <- read.csv("location_params.csv", stringsAsFactors = FALSE)
location_params     <- location_params %>% filter(location == focal.county)

## Combine all parameters
fixed_params        <- c(fixed_params
  , N = location_params[location_params$Parameter == "N", ]$est)

source("variable_params_less_cont.R")

## Run parameters
sim_start  <- variable_params$sim_start
sim_end    <- sim_start + sim.length

param_array <- array(
  data = 0
, dim  = c(nparams, n.mif_runs, 10))
dimnames(param_array)[[3]] <- c(
  "log_lik"
, "R0"
, "Reff"
, "beta0est"
, "E_init"
, "detect_k"
, "detect_mid"
, "theta"
, "theta2"
, "beta_min"
  )  

for (i in 1:nrow(variable_params)) {
  
## Adjust data for parameter choices and the start date for the specific parameter set
if (fit.with == "D" | fit.with == "D_C") {
  
  ## Adjust the data for the current start date
county.data <- deaths %>% 
  mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) %>% 
  filter(day > 0) %>%
  dplyr::select(day, date, deaths, cases) %>% 
  mutate(deaths = deaths - lag(deaths),
         cases = cases - lag(cases)) %>% 
    # three day moving median
  mutate(deaths = mutate(., 
                         deaths_lag = lag(deaths),
                         deaths_lead = lead(deaths)) %>% 
           select(contains("deaths")) %>%
           purrr::pmap_dbl(~median(c(...)))) %>% 
  mutate(cases = mutate(., 
                         cases_lag = lag(cases),
                         cases_lead = lead(cases)) %>% 
           select(contains("cases")) %>%
           purrr::pmap_dbl(~median(c(...)))) %>%
  filter(!is.na(deaths), !is.na(cases))

county.data <- rbind(
  data.frame(
    day    = seq(1:(min(county.data$day) - 1))
  , date   = as.Date(seq(1:(min(county.data$day) - 1)), origin = variable_params[i, "sim_start"])
  , deaths = NA
  , cases  = 0
  )
, county.data 
  )

} else if (fit.with == "H") {
  
## Adjust the data for the current start date
county.data <- hospit %>% mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) 
names(county.data)[2] <- "hosp"  

}
  
mobility <- readRDS("unfolded_daily_clean.rds") %>% 
  filter(county_name == focal.county & state_abbr == focal.state_abbr)  %>%
  select(datestr, sip_prop) %>% 
  filter(!is.na(sip_prop)) %>%
  mutate(day = as.numeric(datestr - as.Date("2019-12-31"))) %>% 
  arrange(day) %>% 
  filter(datestr >= variable_params[i, ]$sim_start) %>%
  mutate(day = day - 14)

mobility <- rbind(
  data.frame(
    datestr  = as.Date(seq(1:(min(mobility$day) - 1)), origin = variable_params[i, ]$sim_start)
  , sip_prop = mean(mobility$sip_prop[1:10])
  , day      = seq(1:(min(mobility$day) - 1))
  )
, mobility
  )

## Remove dates after one week after movement data ends
county.data <- county.data %>% filter(date < (max(mobility$datestr) + 7))

mob.covtab <- covariate_table(
    sip_prop     = c(mobility$sip_prop, rep(mean(mobility$sip_prop[(length(mobility$sip_prop) - 10):length(mobility$sip_prop)]), sim.length - length(mobility$sip_prop))),
    order        = "constant",
    times        = seq(mobility$day[1], sim.length + mobility$day[1] - 1, by = 1),
    intervention = c(rep(0, 115), rep(1, sim.length - 115))
    )

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
  
## Extra params for the new model:
 ## 1) beta0_sigma can't really be fit, so assume for now
 ## 2) fit with no intervention and then simulate later with an intervention
fixed_params["beta0_sigma"] <- 1
fixed_params["beta_catch"]  <- 1
fixed_params["beta_red"]    <- 1 ## No intervention if 1

if (variable_params[i, ]$beta0est == 0) {

## Fit runs in parallel based on number of cores
checktime <- system.time({
  
## Run mifs_local one time as a "burn-in" run (using MCMC terms...)
mifs_local <- foreach(j = 1:n.mif_runs, .combine = c) %dopar%  {
    
library(pomp)
library(dplyr)

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
      beta0       = rlnorm(1, log(0.7), 0.17)
    , E_init      = rpois(1, 2) + 1
    , detect_k    = rlnorm(1, log(0.02), 0.5)
    , detect_mid  = rlnorm(1, log(200), .3)
    , theta       = rlnorm(1, log(5), 0.2)
    , theta2      = rlnorm(1, log(10), 0.3)
    , beta_min    = rlnorm(1, log(0.01), 0.5)
      )
  }
  )
  , Np     = n.mif_particles
  , Nmif   = n.mif_length
  , cooling.fraction.50 = 0.5
  , rw.sd  = rw.sd(
      beta0       = 0.02
    , E_init      = 0.05
    , detect_k    = 0.02
    , detect_mid  = 0.15
    , theta       = 0.005
    , theta2      = 0.01
    , beta_min    = 0.02
  )
        )

testt <- pfilter(mifs_temp)

## Then run mifs_local again as the "sampling" period
covid_mobility %>%
  mif2(
    t0      = 1
  , params  = coef(mifs_temp)
  , Np      = n.mif_particles
  , Nmif    = n.mif_length
  , cooling.fraction.50 = 0.5
  , rw.sd   = rw.sd(
      beta0       = 0.02
    , E_init      = 0.05
    , detect_k    = 0.02
    , detect_mid  = 0.15
    , theta       = 0.005
    , theta2      = 0.01
    , beta_min    = 0.02
  )
        )

}

})

variable_params[i, "beta0est"]       <- mean(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "beta0"), ])
variable_params[i, "E_init"]         <- mean(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "E_init"), ])
variable_params[i, "detect_k"]       <- mean(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "detect_k"), ])
variable_params[i, "detect_mid"]     <- mean(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "detect_mid"), ])
variable_params[i, "theta"]          <- mean(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "theta"), ])
variable_params[i, "theta2"]         <- mean(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "theta2"), ])
variable_params[i, "beta_min"]       <- mean(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "beta_min"), ])
variable_params[i, "beta0est_var"]   <- var(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter  == "beta0"), ])
variable_params[i, "E_init_var"]     <- var(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter  == "E_init"), ])
variable_params[i, "detect_k_var"]   <- var(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter  == "detect_k"), ])
variable_params[i, "detect_mid_var"] <- var(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter  == "detect_mid"), ])
variable_params[i, "theta_var"]      <- var(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter  == "theta"), ])
variable_params[i, "theta2_var"]     <- var(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter  == "theta2"), ])
variable_params[i, "beta_min_var"]   <- var(coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter  == "beta_min"), ])

param_array[i, , "beta0est"]   <- coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "beta0"), ]
param_array[i, , "E_init"]     <- coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "E_init"), ]
param_array[i, , "detect_k"]   <- coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "detect_k"), ]
param_array[i, , "detect_mid"] <- coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "detect_mid"), ]
param_array[i, , "theta"]      <- coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "theta"), ]
param_array[i, , "theta2"]     <- coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "theta2"), ]
param_array[i, , "beta_min"]   <- coef(mifs_local)[which(dimnames(coef(mifs_local))$parameter == "beta_min"), ]

## Check mif2 plots of convergence. Adjust according to parameters
gg.fit <- try(
  mifs_local %>%
  traces() %>%
  melt() %>%
  filter(
    variable == "loglik" | 
    variable == "beta0"  | 
    variable == "E_init" |
    variable == "detect_k" |
    variable == "detect_mid" |
    variable == "theta" |
    variable == "theta2" |
    variable == "beta_min" 
    ) %>% 
  ggplot(aes(x = iteration, y = value, group = L1, colour = factor(L1)))+
  geom_line() +
  guides(color = FALSE) +
  facet_wrap(~variable, scales = "free_y") +
  theme_bw()
, silent = TRUE
)

loglik.out    <- numeric(length(mifs_local))
for (k in seq_along(loglik.out)) {
loglik.out[k] <- pfilter(mifs_local[[k]], 10000)@loglik
}

variable_params[i, "log_lik"]      <- mean(loglik.out)
variable_params[i, "log_lik_var"]  <- var(loglik.out)
param_array[i, , "log_lik"]        <- loglik.out

}

## Simulate from fits
SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid_mobility
    , times        = mob.covtab@times
    , params = c(fixed_params
      , c(
        beta0      = variable_params[i, "beta0est"]
      , E_init     = variable_params[i, "E_init"]
      , detect_k   = variable_params[i, "detect_k"]
      , detect_mid = variable_params[i, "detect_mid"]
      , theta      = variable_params[i, "theta"]
      , theta2     = variable_params[i, "theta2"]
      , beta_min   = variable_params[i, "beta_min"]
      , Ca         = variable_params[i, ]$Ca
      , alpha      = variable_params[i, ]$alpha
      , delta      = variable_params[i, ]$delta
      , mu         = variable_params[i, ]$mu
      , E_init     = variable_params[i, ]$E_init
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


SEIR.sim.s <- SEIR.sim %>% 
  dplyr::filter(.id != "median") %>%
  group_by(day) %>%
  summarize(S = mean(S), D = mean(D))

betat <- variable_params[i, "beta0est"] * exp(log(variable_params[i, "beta_min"]) * mob.covtab@table[1, ])

Reff  <- with(variable_params[i, ], covid_R0(
   beta0est      = betat
 , fixed_params  = c(fixed_params, unlist(variable_params[i, ]))
 , sd_strength   = 1
 , prop_S        = SEIR.sim.s$S / (location_params[location_params$Parameter == "N", ]$est - SEIR.sim.s$D)
  )
  )

## Save an Rds periodically
if (((i / 10) %% 1) == 0) {
 saveRDS(
   list(
    variable_params  = variable_params
  , fixed_params     = fixed_params
  , Reff             = Reff
  , param_array      = param_array
   ), paste(
     paste("output/"
       , paste(focal.county, fit_to_sip, more.params.uncer, fit.minus, Sys.Date(), "cont", "temp", sep = "_")
         , sep = "")
     , "Rds", sep = "."))
}

}

## Save a final Rds
saveRDS(
   list(
    variable_params  = variable_params
  , fixed_params     = fixed_params
  , Reff             = Reff
  , param_array      = param_array
   ), paste(
     paste("output/"
       , paste(focal.county, fit_to_sip, more.params.uncer, fit.minus, Sys.Date(), "cont", "final", sep = "_")
         , sep = "")
     , "Rds", sep = "."))
