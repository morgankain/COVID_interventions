##################################################################
## Use pomp fits to simulate dynamics for various interventions ##
##################################################################
## Using continuous human movement data, fitting to cases and deaths and fitting case detection 

####
## Parameters
####

## NOTE: This version fits to continuous mobility data. The data we are using is private. Please email us or use your own form of contious mobility data.

set.seed(10001)               ## Set seed to recreate fits
# fitting         <- TRUE     ## Small change in pomp objects if fitting or simulating
fit.minus         <- 0        ## Use data until X days prior to the present
usable.cores      <- 3        ## Number of cores to use to fit
import_cases      <- TRUE     ## Use importation of cases?
dt                <- 1/6      ## time step used in fractions of a day
con_theta         <- FALSE    ## Use constrained thetas? 
determ_beta0      <- TRUE     ## Use deterministic beta0?
ci.epidemic_cut   <- 100      ## Criteria of throwing away a stochastic realization as not resulting in an epidemic (total # infected)

## mif2 fitting parameters. 
n.mif_runs        <- 6        ## number of repeated fits (6 used in manuscript, 2 suggested to debug/check code)
n.mif_length      <- 50       ## number of steps (100 used in manuscript, 20 suggested to debug/check code)
n.mif_particles   <- 500      ## number of particles (3000 used in manuscript, 3000 suggested to debug/check code)
n.mif_particles_LL<- 5000     ## number of particles for calculating LL (XXXX used in manuscript, 5000 suggested to debug/check code)
n.mif_rw.sd       <- 0.02     ## particle perturbation (0.02 used in manuscript, 0.02 suggested to debug/check code)
nparams           <- 20       ## number of parameter sobol samples (200 used in manuscript, 5 suggested to debug/check code)
nsim              <- 200      ## number of stochastic epidemic simulations for each fitted beta0 for dynamics (300 used in manuscript, 150 suggested to debug/check code)
# sim.length      <- 200 

sim.plot          <- TRUE
ci.epidemic       <- TRUE

# location of interest
focal.county      <- "Santa Clara" 
focal.state_abbr  <- "CA"
focal.state       <- "California"

mobility.file     <- "unfolded_Jun15.rds"
date_origin       <- as.Date("2019-12-31")
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

# set number of boxes to use for multi compartment classes
nE  <- 3 # must be >= 2
nIa <- 7
nIp <- 2 # must be >= 2
nIm <- nIs <- 5 # must be >= 2
source("COVID_pomp_gammabeta.R")

## Be very careful here, adjust according to your machine's capabilities
  registerDoParallel(cores = usable.cores)
# registerDoParallel(cores = (Sys.getenv("SLURM_NTASKS_PER_NODE")))

## Scrape death data from the NYT github repo to stay up to date or load a previously saved dataset
# deaths <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
# write.csv(deaths, "us-counties.txt")

if (focal.county == "Fulton") {

deaths <- read.csv("us-counties.txt") %>%
  mutate(date = as.Date(date)) %>%
  dplyr::filter((county == focal.county | county == "DeKalb") & state == focal.state) %>%
  group_by(date) %>% summarize(cases = sum(cases), deaths = sum(deaths)) %>%
  mutate(county = focal.county, state = focal.state)

} else {
  
deaths  <- read.csv("us-counties.txt") %>% 
  mutate(date = as.Date(date)) %>% 
  dplyr::filter(county == focal.county & state == focal.state) %>% 
  dplyr::filter(date < max(date) - fit.minus) %>% 
  filter(state == focal.state)
  
}

## Bring in parameters
params <- read.csv("params_multi.csv", stringsAsFactors = FALSE)
params <- params %>% mutate(Value = sapply(est, function(x) eval(parse(text = x))))

fixed_params        <- params$Value
names(fixed_params) <- params$Parameter
if (!import_cases) {fixed_params["import_rate"] <- 0}

## Location-specific parameters
location_params     <- read.csv("location_params.csv", stringsAsFactors = FALSE)
location_params     <- location_params %>% filter(location == focal.county)

## Combine all parameters
if (focal.county == "Fulton") {
fixed_params        <- c(fixed_params, N = location_params[location_params$Parameter == "N", ]$est + 691893)
} else {
fixed_params        <- c(fixed_params, N = location_params[location_params$Parameter == "N", ]$est) 
}
 
if(determ_beta0) {
  fixed_params["beta0_k"] <- 0
}

source("variable_params_less_cont.R")

## Extra needed params
fixed_params["beta_catch"] <- NA
fixed_params["catch_eff"]  <- NA

# bring in mobility data and add days relative to date_origin, with date_origin = 0, units in days
if (focal.county == "Fulton") {
  mobility <- readRDS(mobility.file) %>% 
    dplyr::filter((county_name == focal.county | county_name == "DeKalb") & (state_abbr == focal.state_abbr)) %>%
    dplyr::group_by(datestr) %>%
    dplyr::summarize(sip_prop = mean(sip_prop), .groups = "drop") 
} else {
  mobility <- readRDS(mobility.file) %>% 
    dplyr::filter(county_name == focal.county & state_abbr == focal.state_abbr)
}

# clean up mobility data by taking 3 day moving median and adding days since date_origin
mobility %<>% 
  dplyr::select(datestr, sip_prop) %>% 
  filter(sip_prop != 0) %>% 
  {full_join(., # join with full list of dates so NAs will occur for missing days
             data.frame(datestr = seq(pull(., datestr) %>% min, 
                                      pull(., datestr) %>% max, 
                                      by = "day")))} %>% 
  mutate(day = as.numeric(datestr - date_origin)) %>% 
  arrange(day) %>% 
  # three day moving median
  mutate(sip_prop = mutate(., 
                           sip_prop_lag = lag(sip_prop),
                           sip_prop_lead = lead(sip_prop)) %>% 
           select(contains("sip_prop")) %>%
           purrr::pmap_dbl(~median(c(...), na.rm = T))) %>% 
  dplyr::filter(!is.na(sip_prop)) 

# backfill mobility to date_origin
if (min(mobility$day) > 1) {
  mobility <- rbind(
    data.frame(
      datestr  = seq.Date(date_origin + 1, 
                          min(mobility$datestr) - 1,
                          "day")
      , sip_prop = mean(mobility$sip_prop[1:10])
    ) %>% mutate(day = as.numeric(datestr - date_origin)) 
    , mobility
  )  
}

## Adjust the cases/deaths data to anchor dates from date_origin
county.data <- deaths %>% 
  mutate(day = as.numeric(date - date_origin)) %>% 
  dplyr::select(day, date, deaths, cases) %>% 
  arrange(date) %>%
  # convert to daily deaths and cases
  mutate(deaths = deaths - lag(deaths),
         cases = cases - lag(cases)) %>%
  # convert negatve cases and deaths to zeros
  dplyr::rowwise() %>%
  mutate(cases = max(0, cases),
         deaths = max(0, deaths)) %>%
  ungroup() %>%
  dplyr::filter(!is.na(deaths), !is.na(cases)) %>%
  # fill with NAs from the latest considered sim start date to the data, 
  # NAs are to prevent any oddities with the accumulator variables
  # if the last sim_start date considered occurs when or after the data start, problems will still occur
  {if(max(variable_params[, "sim_start"]) < min(pull(., "date"))) {rbind(
    data.frame(
      date   = seq.Date(variable_params[, "sim_start"] %>% max, 
                        min(pull(., "date")) - 1,
                        "day")
      , deaths = NA
      , cases  = NA
    ) %>% 
      mutate(day = as.numeric(date - date_origin))
    , .
  )} else{ # if the data start before the last sim start date, trim the data to the last sim start and add one row of NAs to deal with accumulator variable
    filter(., date >=  variable_params[, "sim_start"] %>% max) %>% 
      {rbind(data.frame(date = min(pull(., "date")) - 1,
                       deaths = NA, 
                       cases = NA) %>% 
               mutate(day = as.numeric(date - date_origin)),
             .)}
  ## Remove dates after one week after movement data ends
  }} %>% dplyr::filter(date < (max(mobility$datestr) + 7))
  
# create mobility covariate table
mob.covtab <- covariate_table(
  sip_prop       = mobility$sip_prop
  , order        = "constant"
  , times        = mobility$day
  , intervention = rep(0, nrow(mobility))
  , iso_mild_level   = NA
  , iso_severe_level = NA
)

covid_mobility <- pomp(
  data         = county.data %>% select(day, cases, deaths)
  , times      = "day"
  , t0         = 1
  , covar      = mob.covtab
  , rprocess   = euler(sir_step_mobility, delta.t = dt)
  , rmeasure   = {if(con_theta){rmeas_multi_logis_con}else{rmeas_multi_logis_ind}}
  , dmeasure   = {if(con_theta){dmeas_multi_logis_con}else{dmeas_multi_logis_ind}}
  , rinit      = sir_init
  , partrans   = {if(con_theta){par_trans_con}else{par_trans_ind}}
  , accumvars  = accum_names
  , paramnames = param_names
  , statenames = state_names
  , globals    = globs
)  

param_array <- array(
  data = 0
, dim  = c(nparams, n.mif_runs, 10))
dimnames(param_array)[[3]] <- c(
  "log_lik"
, "log_lik_es"
, "beta0est"
, "E_init"
, "detect_k"
, "detect_mid" 
, "detect_max"
, "theta"
, "theta2"
, "beta_min"
  )  

startvals <- array(
  data = 0
, dim  = c(8, n.mif_runs, nparams))
dimnames(startvals)[[1]] <- c(
  "beta0"
, "E_init"
, "detect_k"
, "detect_mid"
, "detect_max"
, "theta"
, "theta2"
, "beta_min"  
)

Reff   <- data.frame(paramset = 0, date = 0, Reff = 0)  ; Reff <- Reff[-1, ]
betat  <- data.frame(paramset = 0, date = 0, betat = 0) ; betat <- betat[-1, ]
detect <- data.frame(paramset = 0, date = 0, detect = 0); detect <- detect[-1, ]

for (i in 1:nrow(variable_params)) {

if (variable_params[i, ]$beta0est == 0) {

  ## Fit runs in parallel based on number of cores
checktime <- system.time({
## Run mifs_local one time as a "burn-in" run (using MCMC terms...)
registerDoRNG(610408799)
mifs_local <- foreach(j = 1:n.mif_runs, .combine = c) %dopar%  {
    
library(pomp)
library(dplyr)
  
start_vals <- c(
      beta0        = rlnorm(1, log(3.5), 0.3)
    , E_init       = rpois(1, 2) + 1
    , detect_k     = rlnorm(1, log(0.1), 0.2)
    , detect_mid   = rlnorm(1, log(60), 0.2)
    , detect_max   = runif(1, 0.1, 0.7)
    , theta        = rlnorm(1, log(1), 0.4)
    , theta2       = ifelse(con_theta, runif(1, 0.1, 0.9) , rlnorm(1, log(1), 0.4))
    , beta_min     = runif(1, 0.3, 0.7) 
)

mifs_temp <- covid_mobility %>%
  mif2(
    t0      = as.numeric(variable_params[i, "sim_start"] - date_origin)
  , params  = c(
    c(fixed_params %>% t() %>% 
        as.data.frame() %>% 
        mutate(d = variable_params[i, "alpha"] * lambda_a + 
                       (1 - variable_params[i, "alpha"]) * variable_params[i, "mu"] * (lambda_p + lambda_m) + 
                       (1 -variable_params[i, "alpha"])*(1 - variable_params[i, "mu"])*(lambda_p + lambda_s)
                   ) %>%
        # convert periods to rates
        mutate(
               gamma    = -1/(nE*dt)*log(1-nE*dt/gamma),
               lambda_a = -1/(nIa*dt)*log(1-nIa*dt/lambda_a),
               lambda_p = -1/(nIp*dt)*log(1-nIp*dt/lambda_p), 
               lambda_m = -1/(nIm*dt)*log(1-nIm*dt/lambda_m),
               lambda_s = -1/(nIs*dt)*log(1-nIs*dt/lambda_s),
               rho_r    = -1/dt*log(1-dt/rho_r),
          ) %>% 
        unlist
    , Ca    = variable_params[i, ]$Ca
    , alpha = variable_params[i, ]$alpha
    , delta = variable_params[i, ]$delta
    , mu    = variable_params[i, ]$mu
    , rho_d = {-1/dt*log(1-dt/variable_params[i, "rho_d"])}
    , detect_t0 = min(county.data[!is.na(county.data$cases), ]$day) - 1
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
  , cooling.fraction.50 = 0.50
  , rw.sd  = rw.sd(
      beta0       = n.mif_rw.sd
    , E_init      = ivp(n.mif_rw.sd)
    , detect_max  = n.mif_rw.sd
    , detect_mid  = n.mif_rw.sd
    , detect_k    = n.mif_rw.sd
    , theta       = n.mif_rw.sd
    , theta2      = n.mif_rw.sd
    , beta_min    = n.mif_rw.sd
  )
        ) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.50) %>%
  mif2(Nmif = n.mif_length, cooling.fraction.50 = 0.50)
  
ll <- replicate(10, mifs_temp %>% pfilter(Np = n.mif_particles_LL) %>% logLik())
ll <- logmeanexp(ll, se = TRUE)
return(list(mifs_temp, ll, start_vals))

}
})

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
best.fit   <- which(loglik.est == max(loglik.est, na.rm = T))

variable_params[i, "beta0est"]   <- coef(mifs_local[[best.fit]])["beta0"]
variable_params[i, "E_init"]     <- coef(mifs_local[[best.fit]])["E_init"]
variable_params[i, "detect_max"] <- coef(mifs_local[[best.fit]])["detect_max"]
variable_params[i, "detect_k"]   <- coef(mifs_local[[best.fit]])["detect_k"]
variable_params[i, "detect_mid"] <- coef(mifs_local[[best.fit]])["detect_mid"]  
variable_params[i, "theta"]      <- coef(mifs_local[[best.fit]])["theta"]
variable_params[i, "theta2"]     <- coef(mifs_local[[best.fit]])["theta2"]
variable_params[i, "beta_min"]   <- coef(mifs_local[[best.fit]])["beta_min"]

all.params <- sapply(mifs_local, coef, simplify = "array")

param_array[i, , "beta0est"]   <- all.params[which(dimnames(all.params)[[1]] == "beta0"), ]
param_array[i, , "E_init"]     <- all.params[which(dimnames(all.params)[[1]] == "E_init"), ]
param_array[i, , "detect_k"]   <- all.params[which(dimnames(all.params)[[1]] == "detect_k"), ]
param_array[i, , "detect_mid"] <- all.params[which(dimnames(all.params)[[1]] == "detect_mid"), ]  
param_array[i, , "detect_max"] <- all.params[which(dimnames(all.params)[[1]] == "detect_max"), ]
param_array[i, , "theta"]      <- all.params[which(dimnames(all.params)[[1]] == "theta"), ]
param_array[i, , "theta2"]     <- all.params[which(dimnames(all.params)[[1]] == "theta2"), ]
param_array[i, , "beta_min"]   <- all.params[which(dimnames(all.params)[[1]] == "beta_min"), ]

variable_params[i, "log_lik"]      <- loglik.est[best.fit]
variable_params[i, "log_lik_se"]   <- loglik.se[best.fit]
param_array[i, , "log_lik"]        <- loglik.est

}

## Simulate from fits
SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid_mobility
    , t0           = as.numeric(variable_params[i, "sim_start"] - date_origin)
    , times        = county.data$day
    , params       = coef(mifs_local[[best.fit]])
    , nsim         = nsim
    , format       = "d"
    , include.data = F
    , seed         = 1001))

SEIR.sim <- SEIR.sim %>%
  mutate(
      date     = date_origin + day
    , paramset = variable_params[i, ]$paramset)

if (ci.epidemic) {
  epi_ids <- SEIR.sim %>% 
    group_by(.id) %>% 
    summarise(total_infect = max(D + R)) %>% 
    filter(total_infect > ci.epidemic_cut*ceiling(variable_params[i, "E_init"])) %>% 
    pull(.id)
  print(paste0("limiting to epidemics, including ", length(epi_ids), " simulations"))
  SEIR.sim %<>% filter(.id %in% epi_ids)
}

SEIR.sim.s <- SEIR.sim %>% 
  group_by(day) %>%
  summarize(S = mean(S), D = mean(D))

SEIR.sim %<>% {
  rbind(.,
        group_by(., day) %>%
          dplyr::select(-.id) %>%
          summarise_all(median) %>%
          mutate(.id = "median"))
}

betat.t      <- variable_params[i, "beta0est"] * exp(log(variable_params[i, "beta_min"])*mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "sip_prop"), which(mob.covtab@times >= min(county.data$day))])

Reff.t       <- with(variable_params[i, ], covid_R0(
   beta0est      = betat.t
 , fixed_params  = c(fixed_params, unlist(variable_params[i, ]))
 , sd_strength   = 1
 , prop_S        = SEIR.sim.s %>% 
   mutate(prop_S = S / (fixed_params["N"] - D)) %>% 
   filter(day <= max(mob.covtab@times)) %>% # limit to same times as mob.covtab
   pull(prop_S)
  )
  )

# detect.t <- c(rep(0, min(county.data[!is.na(county.data$cases), ]$day) - 1)
#               , (variable_params[i, "detect_max"] / 
#                    (1 + exp(-variable_params[i, "detect_k"] * (mob.covtab@times - variable_params[i, "detect_mid"])))
#               )[-seq(1, min(county.data[!is.na(county.data$cases), ]$day) - 1)]
# )
# 
betat.t <- data.frame(
  paramset   = variable_params[i, ]$paramset
  , date     = mob.covtab@times[which(mob.covtab@times >= min(county.data$day))] + date_origin
  , betat    = betat.t
)

Reff.t <- data.frame(
  paramset    = variable_params[i, ]$paramset
  , date     = mob.covtab@times[which(mob.covtab@times >= min(county.data$day))] + date_origin
  , Reff     = Reff.t
)

# detect.t <- data.frame(
#   paramset = variable_params[i, ]$paramset
#   , date     = seq(variable_params[i, ]$sim_start
#                    , variable_params[i, ]$sim_start + (length(mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "sip_prop"), ]) - 1), by = 1)
#   , detect   = detect.t
# )
# 
betat  <- rbind(betat, betat.t)
Reff   <- rbind(Reff, Reff.t)
# detect <- rbind(detect, detect.t)

## Stich together output
if (i == 1) {
  SEIR.sim.f <- SEIR.sim  
} else {
  SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim) 
}  

## Keep track of progress
if (((i / 20) %% 1) == 0) {
  print(paste(round(i / nrow(variable_params), 2)*100, "% Complete", sep = ""))
}

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
       , paste(focal.county, fit.minus, {if(con_theta){"con_theta"}else{"ind_theta"}}, Sys.Date(), "cont", "temp", sep = "_")
         , sep = "")
     , "Rds", sep = "."))

}

if (sim.plot) {

{
  SEIR.sim.f.c.a <- SEIR.sim.f %>% dplyr::select(date, day, .id, cases, paramset)
  SEIR.sim.f.d.a <- SEIR.sim.f %>% dplyr::select(date, day, .id, deaths, paramset)
  SEIR.sim.f.D.a <- SEIR.sim.f %>% dplyr::select(date, day, .id, D, paramset)
  
  SEIR.sim.f.c <- SEIR.sim.f %>% 
    group_by(date, paramset) %>% 
    summarize(
      #    lwr_c = quantile(cases, c(0.100))
      #  , mid_c = quantile(cases, c(0.500))
      #  , upr_c = quantile(cases, c(0.900))
        lwr_c = quantile(cases, c(0.050))
      , mid_c = quantile(cases, c(0.500))
      , upr_c = quantile(cases, c(0.950))
    )
  
  SEIR.sim.f.d <- SEIR.sim.f %>% 
    group_by(date, paramset) %>% 
    summarize(
      lwr_d = quantile(deaths, c(0.050))
      , mid_d = quantile(deaths, c(0.500))
      , upr_d = quantile(deaths, c(0.950))
    )
  
  SEIR.sim.f.D <- SEIR.sim.f %>% 
    group_by(date, paramset) %>% 
    summarize(
      lwr_D = quantile(D, c(0.050))
      , mid_D = quantile(D, c(0.500))
      , upr_D = quantile(D, c(0.950))
    )
  
  # unique(SEIR.sim.f.d$paramset)
  
}
SEIR.sim.f.d.a <- SEIR.sim.f.d.a %>% filter(date < "2020-06-18")
SEIR.sim.f.D.a <- SEIR.sim.f.D.a %>% filter(date < "2020-06-18")
SEIR.sim.f.c.a <- SEIR.sim.f.c.a %>% filter(date < "2020-06-18")
Reff         <- Reff %>% filter(date < "2020-06-18")
detect       <- detect %>% filter(date < "2020-06-18")
{
  gg1 <- ggplot(SEIR.sim.f.d.a) + 
    #  geom_ribbon(aes(
    #    x = date
    #  , ymin = lwr_d
    #  , ymax = upr_d
    #  , group = paramset), alpha = 0.05) +
    geom_line(data = (SEIR.sim.f.d.a %>% filter(.id != "median"))
              , aes(
                  x = date
                , y = deaths
                , group = interaction(paramset, .id))
              , alpha = 0.05, lwd = .2) + 
    geom_line(data = (SEIR.sim.f.d.a %>% filter(.id == "median"))
              , aes(
                x = date
                , y = deaths
                , group = paramset)
              , alpha = 1.0, lwd = 1.0) + 
    geom_point(data = county.data
               , aes(
                 x = date
                 , y = deaths), colour = "dodgerblue4", lwd = 2) + 
    scale_y_continuous(trans = "log10") +
    scale_x_date(labels = date_format("%b"), date_breaks = "2 month") +
    theme(
      axis.text.x = element_text(size = 10)
      , legend.title = element_text(size = 12)
      , plot.title = element_text(size = 12)) +
    xlab("Date") + ylab("Daily Deaths")
  
  gg5 <- ggplot(SEIR.sim.f.D) + 
    #  geom_ribbon(aes(
    #    x = date
    #  , ymin = lwr_D
    #  , ymax = upr_D
    #  , group = paramset), alpha = 0.05) +
    geom_line(data = (SEIR.sim.f.D.a %>% filter(.id != "median"))
              , aes(
                x = date
                , y = D
                , group = interaction(paramset, .id))
              , alpha = 0.05, lwd = .2) + 
    geom_line(data = (SEIR.sim.f.D.a %>% filter(.id == "median"))
              , aes(
                x = date
                , y = D
                , group = paramset)
              , alpha = 1.0, lwd = 1.0) + 
    geom_point(data = deaths
               , aes(
                 x = date
                 , y = deaths), colour = "dodgerblue4", lwd = 2) + 
    scale_y_continuous(trans = "log10") +
    scale_x_date(labels = date_format("%b"), date_breaks = "2 month") +
    theme(
      axis.text.x = element_text(size = 10)
      , legend.title = element_text(size = 12)
      , plot.title = element_text(size = 12)) +
    xlab("Date") + ylab("Total Deaths")
  
  gg2 <- ggplot(SEIR.sim.f.c) + 
    #  geom_ribbon(aes(
    #    x = date
    #  , ymin = lwr_c
    #  , ymax = upr_c
    #  , group = paramset), alpha = 0.05) +
    geom_line(data = (SEIR.sim.f.c.a %>% filter(.id != "median"))
              , aes(
                x = date
                , y = cases
                , group = interaction(paramset, .id))
              , alpha = 0.05, lwd = .2) + 
    geom_line(data = (SEIR.sim.f.c.a %>% filter(.id == "median"))
              , aes(
                x = date
                , y = cases
                , group = paramset)
              , alpha = 1.0, lwd = 1.0) + 
    geom_point(data = county.data
               , aes(
                 x = date
                 , y = cases), colour = "firebrick4", lwd = 2) + 
    scale_y_continuous(trans = "log10") +
    scale_x_date(labels = date_format("%b"), date_breaks = "2 month") +
    theme(
      axis.text.x = element_text(size = 10)
      , legend.title = element_text(size = 12)
      , plot.title = element_text(size = 12)) +
    xlab("Date") + ylab("Observed Cases") 
  
  gg3 <- ggplot(data = Reff, aes(x = date, y = Reff)) + 
    geom_line(aes(group = paramset), alpha = 0.5) +
    theme(
      axis.text.x = element_text(size = 10)
      , legend.title = element_text(size = 12)
      , plot.title = element_text(size = 12)) +
    xlab("Date") + ylab("Reff") +
    scale_x_date(labels = date_format("%b"), date_breaks = "2 month") +
    geom_hline(yintercept = 1, linetype = "dashed", lwd = 0.5)
  
  gg4 <- ggplot(data = detect, aes(x = date, y = detect)) + 
    geom_line(aes(group = paramset), alpha = 0.5) +
    theme(
      axis.text.x = element_text(size = 10)
      , legend.title = element_text(size = 12)
      , plot.title = element_text(size = 12)) +
    scale_x_date(labels = date_format("%b"), date_breaks = "2 month") +
    xlab("Date") + ylab("Case detection proportion")
}
gridExtra::grid.arrange(gg1, gg5, gg2, gg3, gg4, ncol = 1)

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
       , paste(focal.county, fit.minus, {if(con_theta){"con_theta"}else{"ind_theta"}}, Sys.Date(), "cont", "final", sep = "_")
         , sep = "")
     , "Rds", sep = "."))

