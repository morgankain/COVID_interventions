##################################################################
## Use pomp fits to simulate dynamics for various interventions ##
##################################################################
## Using continuous human movement data, fitting to cases and deaths and fitting case detection 

## Note: This script is set up to be able to be run from its own from a saved .Rds.
 ## Thus, packages are set to be loaded etc. (these lines do not need to be run if COVID_fit_cont.R was just run)

## For detailed conceivable intervention scenarios and details on how tosimulate them see:
## "potential_intervention_details_and_plots.R"

####
## Parameters
####

set.seed(10001)
fitting            <- FALSE   ## Small change in pomp objects if fitting or simulating
## TRUE if COVID_fit previously run, FALSE if COVID_fit was just run and global environment is still full
use.rds            <- TRUE    
#rds.name           <- "output/Santa Clara_TRUE_FALSE_0_2020-06-02_cont_temp_NZ_sig_exp.Rds"
#rds.name           <- "output/Santa Clara_TRUE_FALSE_0_2020-06-04_cont_temp_NZ_sig_exp_gamma.Rds"
rds.name           <- "output/Miami-Dade_TRUE_FALSE_0_2020-06-05_cont_final_NZ_sig_exp_gamma.Rds"
mobility.file      <- "mobility/unfolded_Jun03.rds"
more.params.uncer  <- FALSE   ## Fit with more (FALSE) or fewer (TRUE) point estimates for a number of parameters
nsim               <- 10     ## Number of epidemic simulations for each parameter set
fit.E0             <- TRUE    ## Was E0 also fit?
usable.cores       <- 1       ## Number of cores to use to fit

int.beta0_sigma    <- 0.16       ## heterogeneity value

## Where to simulate from
sir_init.mid       <- FALSE         ## Starts the epidemic from some non-zero timepoint
sir_init.mid.t     <- "2020-05-28"  ## Date to simulate forward from

## Sim and plotting details
sim_length         <- 200     ## How many days to run the simulation
state.plot         <- "D"     ## State variable for plotting (Hospit [H], Death [D], or Cases [C])
plot.log10         <- TRUE    ## Plot on a log10 scale or not
print.plot         <- FALSE

counter.factual   <- FALSE    ## If true do a special analysis that ignores a lot of these other parameters
cf.type           <- "may1"   ## Specifically modeled counterfactual analyses: no_int, delay, may1 coded for now

####
## Intervention parameters
####
int.movement       <- "post"     ## Shape of human movement after the data ends ::: pre, post, mid
int.type           <- "tail"    ## Extra intervention apart from the movement  ::: none, inf_iso, tail
int.init           <- "2020-06-08"
int.end            <- "2020-08-01"

## Test and Isolate parameters
iso_mild_level     <- 0.2
iso_severe_level   <- 0.2

## Tail chopping parameters
int.beta_catch     <- 1000      ## beta0 values caught by intervention
int.catch_eff      <- 0         ## effectiveness at catching beta0 values above the beta_catch (0 - 1)
# int.beta_red       <- 1       ## new values after those that are caught

# crude way to ensure int.beta_catch > 0
if(int.beta_catch <= 0){
  int.beta_catch <- 0.0001
}

loglik.thresh      <- 2       ## Keep parameter sets with a likelihood within top X loglik units
params.all         <- TRUE    ## Keep all fitted parameters above loglik thresh?...
nparams            <- 50      ## ...if FALSE, pick a random subset for speed

fit.with           <- "D_C"   ## Fit with D (deaths) or H (hospitalizations) -- Need your own data for H -- or both deaths and cases (D_C).
fit_to_sip         <- TRUE    ## Fit beta0 and shelter in place strength simultaneously?
meas.nb            <- TRUE    ## Negative binomial measurement process?
import_cases       <- FALSE   ## Use importation of cases?
fit.minus          <- 0       ## Use data until X days prior to the present

fixed.E0           <- !fit.E0
detect.logis       <- T

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
  , "data.table"
  , "doRNG")

## load packages. Install all packages that return "FALSE"
lapply(needed_packages, require, character.only = TRUE)

## Be very careful here, adjust according to your machine's capabilities
  registerDoParallel(cores = usable.cores)
# registerDoParallel(cores = (Sys.getenv("SLURM_NTASKS_PER_NODE")))
  
## Theme for pretty plots
source("ggplot_theme.R")
## Bring in pomp objects.
source("COVID_pomp_gammabeta_int.R")
  
if (fit.with == "D_C" | fit.with == "D") {
## Scrape death data from the NYT github repo to stay up to date or load a previously saved dataset
#deaths <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
deaths  <- read.csv("us-counties.txt")
deaths  <- deaths %>% mutate(date = as.Date(date)) %>% dplyr::filter(county == focal.county)
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
  
## Load the previously saved fits
 ## If COVID_fit_cont.R was just run, use parameers already stored in the global env 
if (use.rds) {
prev.fit         <- readRDS(rds.name)
variable_params  <- prev.fit[["variable_params"]]
fixed_params     <- prev.fit[["fixed_params"]]
}
  
## drop the rows that have 0s for likelihood (in case exited prematurely) 
 ## and keep only the best fits as defined by loglik
variable_params <- variable_params %>% 
  filter(log_lik != 0) %>% 
 # filter(log_lik == max(log_lik)) ## only used for manuscript dyanmic trajectory plots
 filter(log_lik > (max(log_lik) - loglik.thresh))

## For debug plots just pick one paramset
# variable_params <- variable_params[variable_params$paramset == 2, ]

## A few adjustments for preprint counterfactuals and other scenarios
# 1) Start later
# variable_params <- variable_params %>% mutate(int_start2 = int_start2 + 7) 
# 2) adjustment to social distancing strength
# variable_params <- variable_params %>% mutate(soc_dist_level_sip = 0.2)

## Adjust variable params for the simulation scenario. !!! Some adjustment to be made
# variable_params <- variable_params %>% 

if (!params.all) {
  variable_params <- variable_params[sample(1:nrow(variable_params), nparams), ]
}

####
## Simulate: loop over rows of variable_params
####

Reff   <- data.frame(paramset = 0, date = 0, Reff = 0)  ; Reff <- Reff[-1, ]
betat  <- data.frame(paramset = 0, date = 0, betat = 0) ; betat <- betat[-1, ]
detect <- data.frame(paramset = 0, date = 0, detect = 0); detect <- detect[-1, ]

for (i in 1:nrow(variable_params)) {
  
county.data <- deaths %>% 
  mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) %>% 
  dplyr::filter(day > 0) %>%
  dplyr::select(day, date, deaths, cases) %>% 
  mutate(deaths = deaths - lag(deaths),
         cases = cases - lag(cases)) %>% 
    # three day moving median
  mutate(deaths = mutate(., 
                         deaths_lag = lag(deaths),
                         deaths_lead = lead(deaths)) %>% 
           dplyr::select(contains("deaths")) %>%
           purrr::pmap_dbl(~median(c(...)))) %>% 
  mutate(cases = mutate(., 
                         cases_lag = lag(cases),
                         cases_lead = lead(cases)) %>% 
           dplyr::select(contains("cases")) %>%
           purrr::pmap_dbl(~median(c(...)))) %>%
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

mobility <- readRDS(mobility.file) %>% 
  dplyr::filter(county_name == focal.county & state_abbr == focal.state_abbr)  %>%
  dplyr::select(datestr, sip_prop) %>% 
  dplyr::filter(!is.na(sip_prop)) %>%
  mutate(day = as.numeric(datestr - as.Date("2019-12-31"))) %>% 
  arrange(day) %>% 
  # three day moving median
  mutate(sip_prop = mutate(., 
                         sip_prop_lag = lag(sip_prop),
                         sip_prop_lead = lead(sip_prop)) %>% 
           select(contains("sip_prop")) %>%
           purrr::pmap_dbl(~median(c(...)))) %>% 
  dplyr::filter(datestr >= variable_params[i, ]$sim_start) %>%
  mutate(day = day - as.numeric((variable_params[i, ]$sim_start - as.Date("2019-12-31")))) %>%
  dplyr::filter(!is.na(sip_prop)) %>% 
  dplyr::filter(day > 0)

if (min(mobility$day) > 1) {
mobility <- rbind(
  data.frame(
    datestr  = as.Date(seq(1:(min(mobility$day) - 1)), origin = variable_params[i, ]$sim_start)
  , sip_prop = mean(mobility$sip_prop[1:10])
  , day      = seq(1:(min(mobility$day) - 1))
  )
, mobility
  )  
}
 
mobility <- mobility %>% filter(sip_prop != 0)

## Remove dates after one week after movement data ends
county.data <- county.data %>% dplyr::filter(date < (max(mobility$datestr) + 7))

if (counter.factual) {
  
if (cf.type == "no_int") {
mob.covtab <- covariate_table(
   sip_prop         = rep(mean(head(mobility$sip_prop, 7)), sim_length)
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
 , intervention     = rep(0, sim_length))  
} else if (cf.type == "delay") {
mob.covtab <- covariate_table(
   sip_prop         = c(rep(mean(head(mobility$sip_prop, 7)), 7), mobility$sip_prop
     , rep(mean(tail(mobility$sip_prop, 3)), sim_length - (length(mobility$sip_prop) + 7))) 
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
 , intervention     = rep(0, sim_length))  
} else if (cf.type == "may1") { 
  
int.begin    <- as.numeric((as.Date("2020-05-01") - variable_params[i, ]$sim_start))
int.phase1   <- nrow(mobility) + (int.begin - nrow(mobility))
  
mob.covtab <- covariate_table(
   sip_prop         = c(
      mobility$sip_prop[1:int.begin]
    , rep(mean(head(mobility$sip_prop, 7)), sim_length - int.phase1)
      )
 , order            = "constant"
 , times            = seq(1, sim_length, by = 1)
 , detect_t0        = min(county.data[!is.na(county.data$cases), ]$day) - 1
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
 , intervention     = rep(0, sim_length))  
}
  
} else {
  
####
## Set up intervention scenario
####

int.begin    <- as.numeric((as.Date(int.init) - variable_params[i, ]$sim_start))
int.duration <- as.Date(int.end) - as.Date(int.init)

int.phase1 <- nrow(mobility) + (int.begin - nrow(mobility))
int.phase2 <- as.numeric(int.duration)
int.phase3 <- sim_length - int.phase1 - int.phase2

## Covariate table with all of the intervention scenarios
{
mob.covtab <- covariate_table(
  ## First, can add whatever sip_prop into the future
  
   sip_prop      = {
     if (int.movement == "pre") {
     c(
      mobility$sip_prop
    , rep(mean(tail(mobility$sip_prop, 3)), int.phase1 - length(mobility$sip_prop))
    , rep(mean(head(mobility$sip_prop, 7)), sim_length - int.phase1)
      )       
     } else if (int.movement == "post") {
     c(
      mobility$sip_prop
    , rep(mean(tail(mobility$sip_prop, 3)), int.phase1 - length(mobility$sip_prop))
    , rep(mean(tail(mobility$sip_prop, 3)), sim_length - int.phase1)
      )
     } else if (int.movement == "mid") {
     c(
      mobility$sip_prop
    , rep(mean(tail(mobility$sip_prop, 3)), int.phase1 - length(mobility$sip_prop))
    , rep((mean(head(mobility$sip_prop, 7)) + mean(tail(mobility$sip_prop, 3))) / 2, sim_length - int.phase1)
      )         
     }
    
   }
  
 , order         = "constant"
 , times         = seq(1, sim_length, by = 1)
 , detect_t0     = min(county.data[!is.na(county.data$cases), ]$day) - 1
  
 , iso_mild_level   = iso_mild_level
 , iso_severe_level = iso_severe_level
  
  ## Next can add an intervention that stops large gatherings
  
 , intervention  = {
   if (int.type == "none") {
   c(
    rep(0, int.phase1)
  , rep(0, int.phase2)
  , rep(0, int.phase3)
  )     
   } else if (int.type == "tail") {
   c(
    rep(0, int.phase1)
  , rep(1, int.phase2)
  , rep(0, int.phase3)
   ) 
   } else if (int.type == "inf_iso") {
   c(
    rep(0, int.phase1)
  , rep(2, int.phase2)
  , rep(2, int.phase3)
   )   
  }
   }
    )
}

}

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
 , globals    = trunc_gamma
)

## Second bit of the intervention will affect these parameters:
fixed_params["beta0_sigma"] <- int.beta0_sigma  ## heterogeneity value
fixed_params["beta_catch"]  <- int.beta_catch  ## beta0 values caught by intervention
fixed_params["catch_eff"]   <- int.catch_eff  ## beta0 values caught by intervention
# fixed_params["beta_red"]    <- int.beta_red  ## new values after those that are caught


SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid_mobility
    , times        = mob.covtab@times
    , params = c(fixed_params
      # prev.fit[["fixed_params"]]
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
      
SEIR.sim <- SEIR.sim %>%
  mutate(
      date     = round(as.Date(day, origin = variable_params[i, ]$sim_start))
    , paramset = variable_params[i, ]$paramset)

# SEIR.sim %>% filter(date == sir_init.mid.t, .id == "median")

## Need to think about the most principled way of simulating forward from the present. For now just take the median up until the present as the
 ## starting values, and project from there
if (sir_init.mid) {
  
new.startvals <- SEIR.sim %>% filter(.id == "median" & date == sir_init.mid.t)  
new.days      <- seq(as.numeric((as.Date(sir_init.mid.t) - variable_params[i, ]$sim_start)), sim_length)

mob.covtab.c <- mob.covtab
mob.covtab.c@times <- mob.covtab@times[1:(length(new.days))]
mob.covtab.c@table <- mob.covtab@table[, min(new.days):sim_length]

covid_mobility <- pomp(
   data       = county.data
 , times      = "day"
 , t0         = 1
 , covar      = mob.covtab.c
 , rprocess   = euler(sir_step_mobility, delta.t = 1/6)
 , rmeasure   = rmeas_multi_logis
 , dmeasure   = dmeas_multi_logis
 , rinit      = sir_init_mid
 , partrans   = par_trans
 , accumvars  = accum_names
 , paramnames = c(param_names, mid_init_param_names)
 , statenames = state_names
 , globals    = trunc_gamma
)

SEIR.sim.c <- do.call(
  pomp::simulate
  , list(
    object         = covid_mobility
    , times        = mob.covtab.c@times
    , params = c(prev.fit[["fixed_params"]]
      , c(
        beta0      = variable_params[i, "beta0est"]
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
        )
      , unlist(c(
        E_init     = round(new.startvals$E)
      , S0         = round(new.startvals$S)
      , Ia0        = round(new.startvals$Ia)
      , Ip0        = round(new.startvals$Ip)
      , Is0        = round(new.startvals$Is)
      , Im0        = round(new.startvals$Im)
      , Hr0        = round(new.startvals$Hr)
      , Hd0        = round(new.startvals$Hd)
      , R0         = round(new.startvals$R)
      , D0         = round(new.startvals$D)
      )))
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

SEIR.sim.c <- SEIR.sim.c %>%
  mutate(
      date     = round(as.Date(day, origin = sir_init.mid.t))
    , paramset = variable_params[i, ]$paramset)

# SEIR.sim.c %>% filter(date == as.Date(sir_init.mid.t) + 3, .id == "median")
# SEIR.sim %>% filter(date == as.Date(sir_init.mid.t) + 2, .id == "median")

}

SEIR.sim.s <- SEIR.sim %>% 
  dplyr::filter(.id != "median") %>%
  group_by(day) %>%
  summarize(S = mean(S), D = mean(D))

betat.t      <- variable_params[i, "beta0est"] * exp(log(variable_params[i, "beta_min"])*mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "sip_prop"), ])

Reff.t       <- with(variable_params[i, ], covid_R0(
   beta0est      = betat.t
 , fixed_params  = c(fixed_params, unlist(variable_params[i, ]))
 , sd_strength   = 1
 , prop_S        = SEIR.sim.s$S / (prev.fit[["fixed_params"]]["N"] - SEIR.sim.s$D)
  )
  )

detect.t <- c(rep(0, mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "detect_t0"), 1])
  , (variable_params[i, "detect_max"] / 
    (1 + exp(-variable_params[i, "detect_k"] * (mob.covtab@times - variable_params[i, "detect_mid"])))
    )[-seq(1, mob.covtab@table[which(dimnames(mob.covtab@table)[[1]] == "detect_t0"), 1])]
)

betat.t <- data.frame(
  paramset = variable_params[i, ]$paramset
, date     = seq(variable_params[i, ]$sim_start
, variable_params[i, ]$sim_start + (sim_length - 1), by = 1)
, betat    = betat.t
  )

Reff.t <- data.frame(
  paramset = variable_params[i, ]$paramset
, date     = seq(variable_params[i, ]$sim_start
, variable_params[i, ]$sim_start + (sim_length - 1), by = 1)
, Reff     = Reff.t
  )

detect.t <- data.frame(
  paramset = variable_params[i, ]$paramset
, date     = seq(variable_params[i, ]$sim_start
, variable_params[i, ]$sim_start + (sim_length - 1), by = 1)
, detect   = detect.t
  )

if (sir_init.mid) {
SEIR.sim.c <- left_join(SEIR.sim.c, detect.t[, c(2, 3)], by = "date")
SEIR.sim.c <- SEIR.sim.c %>% mutate(cases = I_new_sympt * detect)
}

betat  <- rbind(betat, betat.t)
Reff   <- rbind(Reff, Reff.t)
detect <- rbind(detect, detect.t)

## Stich together output
if (sir_init.mid) {
if (i == 1) {
 SEIR.sim.f <- SEIR.sim.c 
} else {
 SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim.c) 
}
} else {
if (i == 1) {
 SEIR.sim.f <- SEIR.sim  
} else {
 SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim) 
}  
}

## Keep track of progress
if (((i / 20) %% 1) == 0) {
  print(paste(round(i / nrow(variable_params), 2)*100, "% Complete", sep = ""))
}

}

## Because of the staggered start dates there is some oddity at the end of the simulation, just
 ## remove those few days
SEIR.sim.f <- SEIR.sim.f %>% filter(date < min(variable_params$sim_start + sim_length))

####
## Summary and plotting
####

SEIR.sim.f.c.a <- SEIR.sim.f %>% dplyr::select(date, day, .id, cases, paramset)
SEIR.sim.f.d.a <- SEIR.sim.f %>% dplyr::select(date, day, .id, deaths, paramset)
SEIR.sim.f.D.a <- SEIR.sim.f %>% dplyr::select(date, day, .id, D, paramset)

SEIR.sim.f.c <- SEIR.sim.f %>% 
  group_by(date, paramset) %>% 
  summarize(
#    lwr_c = quantile(cases, c(0.100))
#  , mid_c = quantile(cases, c(0.500))
#  , upr_c = quantile(cases, c(0.900))
    lwr_c = quantile(cases, c(0.025))
  , mid_c = quantile(cases, c(0.500))
  , upr_c = quantile(cases, c(0.975))
  )

SEIR.sim.f.d <- SEIR.sim.f %>% 
  group_by(date, paramset) %>% 
  summarize(
    lwr_d = quantile(deaths, c(0.100))
  , mid_d = quantile(deaths, c(0.500))
  , upr_d = quantile(deaths, c(0.900))
  )

SEIR.sim.f.D <- SEIR.sim.f %>% 
  group_by(date, paramset) %>% 
  summarize(
    lwr_D = quantile(D, c(0.100))
  , mid_D = quantile(D, c(0.500))
  , upr_D = quantile(D, c(0.900))
  )

# unique(SEIR.sim.f.d$paramset)

# SEIR.sim.f.d <- SEIR.sim.f.d %>% filter(date < "2020-03-15")
# SEIR.sim.f.D <- SEIR.sim.f.D %>% filter(date < "2020-03-15")
# SEIR.sim.f.c <- SEIR.sim.f.c %>% filter(date < "2020-03-15")
# Reff         <- Reff %>% filter(date < "2020-03-15")
# detect       <- detect %>% filter(date < "2020-03-15")

## Actual plots
if (print.plot) {
test.paramset <- 2
{
gg1 <- ggplot(SEIR.sim.f.d
  #[SEIR.sim.f.d$paramset == test.paramset, ]
  ) + 
  geom_ribbon(aes(
    x = date
  , ymin = lwr_d
  , ymax = upr_d
  , group = paramset), alpha = 0.05) +
  geom_line(aes(
      x = date
    , y = mid_d
    , group = paramset)
      , alpha = 0.50, lwd = .5) + 
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

gg5 <- ggplot(SEIR.sim.f.D
  #[SEIR.sim.f.D$paramset == test.paramset, ]
  ) + 
  geom_ribbon(aes(
    x = date
  , ymin = lwr_D
  , ymax = upr_D
  , group = paramset), alpha = 0.05) +
  geom_line(aes(
      x = date
    , y = mid_D
    , group = paramset)
      , alpha = 0.50, lwd = .5) + 
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

gg2 <- ggplot(SEIR.sim.f.c
  #[SEIR.sim.f.c$paramset == test.paramset, ]
  ) + 
  geom_ribbon(aes(
    x = date
  , ymin = lwr_c
  , ymax = upr_c
  , group = paramset), alpha = 0.05) +
  geom_line(aes(
      x = date
    , y = mid_c
    , group = paramset)
      , alpha = 0.50, lwd = .5) + 
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

SEIR.sim.f.c.a <- SEIR.sim.f.c.a %>% filter(date < "2020-04-01")

ggplot(SEIR.sim.f.c.a) + 
  geom_line(aes(
      x = date
    , y = cases
    , group = .id)
      , alpha = 0.10, lwd = .2) + 
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

SEIR.sim.f.c <- SEIR.sim.f.c %>% filter(date < "2020-04-01")

ggplot(SEIR.sim.f.c) + 
  geom_ribbon(aes(
    x = date
  , ymin = lwr_c
  , ymax = upr_c
  , group = paramset), alpha = 0.50) +
  geom_line(aes(
      x = date
    , y = mid_c
    , group = paramset)
      , alpha = 0.50, lwd = .5) + 
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
