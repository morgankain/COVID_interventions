##################################################################
## Use pomp fits to simulate dynamics for various interventions ##
##################################################################

## Note: This script is set up to be able to be run from its own from a saved .Rds.
 ## Thus, packages are set to be loaded etc. (these lines do not need to be run if COVID_fit.R was just run)

####
## Parameters
####
set.seed(10001)
focal.county       <- "Santa Clara"
fitting            <- FALSE   ## Small change in pomp objects if fitting or simulating
## TRUE if COVID_fit previously run, FALSE if COVID_fit was just run and global environment is still full
use.rds            <- TRUE    
rds.name           <- "output/Santa Clara_FALSE_FALSE_0_2020-04-25_final.Rds"
more.params.uncer  <- FALSE   ## Fit with more (FALSE) or fewer (TRUE) point estimates for a number of parameters
nsim               <- 200     ## Number of epidemic simulations for each parameter set
fit.E0             <- FALSE   ## Was E0 also fit?
meas.nb            <- TRUE    ## Negative binomial measurement process?
## Intervention scenarios. Only one at a time allowed right now! 
inf_iso            <- FALSE   ## Do we ever reduce from shelter in place to some form of strong/moderate social distancing?
test_and_isolate_s <- 0.2     ## Additional proportional reduction of severe cases under test and isolate
test_and_isolate_m <- 0.2     ## Additional proportional reduction of mild cases under test and isolate
light              <- FALSE   ## Lightswitch method
red_shelt.t        <- 76      ## Time shelter in place changes to a reduced form (inf_iso or light) -- Now June 1 (day 76 of shelter in place orders)
red_shelt.s        <- 0.5     ## New social dist strength after time red_shelt.t
thresh_H.start     <- 15      ## Threshold when lightswtich turns on (when we get higher than this)
thresh_H.end       <- 5       ## Threshold when lightswtich turns off (when we drop from above to this value)
sim_length         <- 400     ## How many days to run the simulation
state.plot         <- "D"     ## State variable for plotting (Hospit [H], Death [D], or Cases [C])

loglik.thresh      <- 2       ## Keep parameter sets with a likelihood within top X loglik units
params.all         <- FALSE   ## Keep all fitted parameters above loglik thresh?...
nparams            <- 50      ## ...if FALSE, pick a random subset for speed
plot.log10         <- FALSE   ## Plot on a log10 scale or not
fit.with           <- "D"     ## Needed again to organize data correctly 
fit.minus          <- 0       ## Use data until X days prior to the present

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
  , "data.table"
  )

## load packages. Install all packages that return "FALSE"
lapply(needed_packages, require, character.only = TRUE)

## Theme for pretty plots
source("ggplot_theme.R")
## Bring in pomp objects.
source("COVID_pomp.R")

if (fit.with == "D") {
  ## Scrape death data from the NYT github repo to stay up to date or load a previously saved dataset
   ## Only used in this script for plotting purposes
#deaths <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
deaths  <- read.csv("us-counties.txt")
deaths  <- deaths %>% mutate(date = as.Date(date)) %>% filter(county == focal.county)
deaths  <- deaths %>% dplyr::filter(date < max(date) - fit.minus)
} else if (fit.with == "H") {
## See warning notes in COVID_fit.R
hospit  <- read.csv("contra_costa/ccc_data.csv")
hospit  <- hospit %>% 
  mutate(date = as.Date(REPORT_DATE)) %>% 
  filter(CURRENT_HOSPITALIZED != "NULL") %>% 
  mutate(ch = as.numeric(as.character(CURRENT_HOSPITALIZED))) %>% 
  dplyr::select(date, ch)
hospit    <- hospit %>% dplyr::filter(date < max(date) - fit.minus)  
}

## Load the previously saved fits
 ## If COVID_fit.R was just run, use parameers already stored in the global env 
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

## A few adjustments for preprint counterfactuals and other scenarios
# 1) Start later
# variable_params <- variable_params %>% mutate(int_start2 = int_start2 + 7) 
# 2) adjustment to social distancing strength
# variable_params <- variable_params %>% mutate(soc_dist_level_sip = 0.2)

## Adjust variable params for the simulation scenario
variable_params <- variable_params %>% 
  mutate(
    int_length2        = red_shelt.t
  , iso_start          = int_start2 + red_shelt.t
  , soc_dist_level_red = red_shelt.s
  ) 

if (!params.all) {
  variable_params <- variable_params[sample(1:nrow(variable_params), nparams), ]
}

####
## Simulate: loop over rows of variable_params
####

for (i in 1:nrow(variable_params)) {
  
## Adjust data for parameter choices and the start date for the specific parameter set
if (fit.with == "D") {
  
## Adjust the data for the current start date
county.data <- deaths %>% 
  mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) %>% 
  filter(day > 0) %>%
  dplyr::select(day, date, deaths) %>% 
  mutate(deaths_cum = deaths) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum)) %>% 
  replace_na(list(deaths = 0)) %>%
  dplyr::select(-deaths_cum)
  
} else {
  
## Adjust the data for the current start date
county.data <- hospit %>% mutate(day = as.numeric(date - variable_params[i, ]$sim_start))   
names(county.data)[2] <- "hosp"  

}

## Create intervention covariate table for the full forecast
intervention.forecast <- with(variable_params[i, ], {

 covariate_table(
  day              = 1:sim_length
, intervention     = c(
      # No intervention until intervention start time
    rep(0, int_start1 - sim_start)                   
      # Intervention style 1
  , rep(1, int_length1)
      # Intervention style 2
  , rep(2, int_length2)
      # Post intervention close
  , {
    if (!inf_iso & !light) {
      rep(0, sim_length - (int_start2 - sim_start) - int_length2)
    } else if (inf_iso & !light) {
      rep(2, sim_length - (int_start2 - sim_start) - int_length2)  
    } else if (!inf_iso & light) {
      rep(3, sim_length - (int_start2 - sim_start) - int_length2)        
    }
    }
  )
, thresh_int_level = rep(soc_dist_level_sip, sim_length)
, back_int_level   = rep(red_shelt.s, sim_length)
, isolation = { 
  
if (!inf_iso) {
  
    rep(0, sim_length)
  
} else {
  
  c(
    rep(0, iso_start - sim_start)  
  , rep(1, sim_length - (iso_start - sim_start))
  )
  
}}
, iso_severe_level = rep(test_and_isolate_s, sim_length)      # % of contats that severe cases maintain
, iso_mild_level   = rep(test_and_isolate_m, sim_length)   # % of contats that mild cases maintain
, soc_dist_level_wfh = rep(soc_dist_level_wfh, sim_length) 
, soc_dist_level_sip = {
if (!inf_iso) {
  rep(soc_dist_level_sip, sim_length) 
} else {
  c(
   rep(soc_dist_level_sip, (int_start1 - sim_start) + int_length1 + int_length2)
 , rep(soc_dist_level_red, sim_length - ((int_start1 - sim_start) + int_length1 + int_length2))
  )
}
}
, thresh_H_start   = rep(thresh_H.start, sim_length) 
, thresh_H_end     = rep(thresh_H.end, sim_length)
, order            = "constant"
, times            = "day"
  )

})

## Need to set up the pomp object, but no fitting, rmeas and dmeas using D or H is irrelevant
covid.fitting <- county.data %>%
  pomp(
    time       = "day"
  , t0         = 1
  , covar      = intervention.forecast
  , rprocess   = euler(sir_step, delta.t = 1/6)
  , rmeasure   = { 
   if (fit.with == "D") { rmeas_deaths } else if (fit.with == "H") { rmeas_hosp }
  }
  , dmeasure   = {
   if (fit.with == "D") { dmeas_deaths } else if (fit.with == "H") { dmeas_hosp }
  }
  , rinit      = sir_init
  , partrans   = par_trans
  , accumvars  = accum_names
  , paramnames = param_names
  , statenames = state_names) 

## Simulate from the beginning and project forward with this beta0. 
  ## Could consider prefer pmcmc for uncertainty in beta0.
if (!more.params.uncer) {

SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid.fitting
    , times        = intervention.forecast@times
    , params = c(fixed_params
    , c(
      beta0              = variable_params[i, "beta0est"]
    , soc_dist_level_sip = variable_params[i, "soc_dist_level_sip"]
    , Ca                 = variable_params[i, ]$Ca
    , alpha              = variable_params[i, ]$alpha
    , delta              = variable_params[i, ]$delta
    , mu                 = variable_params[i, ]$mu
      )
  , if (fit.E0) {
    c(E_init = variable_params[i, "E_init"])
  } else {
    c(E0 = variable_params[i, ]$E0)
  }
      )
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
    object         = covid.fitting
    , times        = intervention.forecast@times
    , params = c(fixed_params
    , c(
      beta0              = variable_params[i, "beta0est"]
    , soc_dist_level_sip = variable_params[i, "soc_dist_level_sip"]
    , Ca                 = variable_params[i, ]$Ca
    , alpha              = variable_params[i, ]$alpha
    , delta              = variable_params[i, ]$delta
    , mu                 = variable_params[i, ]$mu
    , lambda_a           = variable_params[i, ]$lambda_a
    , lambda_s           = variable_params[i, ]$lambda_s
    , lambda_m           = variable_params[i, ]$lambda_m  
      )
  , if (fit.E0) {
    c(E_init = variable_params[i, "E_init"])
  } else {
    c(E0 = variable_params[i, ]$E0)
  }
      )
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

## Plot of a single run
## source("single_run_plot.R")

SEIR.sim <- SEIR.sim %>%
  mutate(date = round(as.Date(day, origin = variable_params[i, ]$sim_start)))

SEIR.sim <- SEIR.sim %>% 
  mutate(
    paramset = variable_params[i, ]$paramset
  )

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

}

## Because of the staggered start dates there is some oddity at the end of the simulation, just
 ## remove those few days
SEIR.sim.f <- SEIR.sim.f %>% filter(date < min(variable_params$sim_start + sim_length))

## Summary of observable cases to plot against case data in the county
SEIR.sim.f <- SEIR.sim.f %>% mutate(C = Is + Im + H)

####
## Summary and plotting
####

state.plot    <- "D"

SEIR.sim.f.m <- SEIR.sim.f %>% 
  dplyr::filter(.id == "median") %>% 
  dplyr::select(state.plot, date, paramset)

## Choose some desired PI for each run
SEIR.sim.f.s <- SEIR.sim.f %>% 
  dplyr::filter(.id != "median") %>%
  dplyr::group_by(date) %>%
  dplyr::summarize(
#   lwr = quantile(get(state.plot), c(0.025))
# , est = quantile(get(state.plot), c(0.50))
# , upr = quantile(get(state.plot), c(0.975))
    lwr = quantile(get(state.plot), c(0.05))
  , est = quantile(get(state.plot), c(0.50))
  , upr = quantile(get(state.plot), c(0.95))
  )

## Attempt at dynamic labels for plots depending on parameter choices
if (inf_iso & !light) {
  plot.title <- paste("Relax shelter in place to", red_shelt.s, "of baseline contacts on", red_shelt.t, sep = " ")
} else if (!inf_iso & light) {
  plot.title <- paste("Switch strong distancing on when H =", thresh_H.start, "Move to weak distancing when H <"
  , thresh_H.end, sep = " ")
} else {
  plot.title <- paste("Lift shelter in place on", red_shelt.t, sep = " ")
}

gg.1 <- ggplot(SEIR.sim.f.s) + 
  geom_ribbon(aes(x = date, ymin = lwr, ymax = upr), colour = NA, fill = "dodgerblue4", alpha = 0.25) +
  geom_line(data = SEIR.sim.f.m, aes(x = date, y = get(state.plot), group = paramset), colour = "grey50", alpha = 0.35) +
  scale_x_date(labels = date_format("%Y-%b"), date_breaks = "1 month") +
  geom_line(aes(x = date, y = est), colour = "dodgerblue4", lwd = 1.5) + 
  xlab("Date") + 
  ylab(state.plot) +
  guides(color = FALSE) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12)
  , legend.title = element_text(size = 12)
  , plot.title = element_text(size = 12)) +
  ggtitle(plot.title) 

if (plot.log10) {
gg.1 <- gg.1 + scale_y_log10()
}

## No hospit data 
if (state.plot == "H") {
  (gg.1 <- gg.1 + geom_point(data = county.data, aes(date, hosp)))
} else if (state.plot == "D") {
  (gg.1 <- gg.1 + geom_point(data = deaths, aes(date, deaths)))
} else if (state.plot == "C") {
  ## If case data can uncomment this
# (gg.1 <- gg.1 + geom_point(data = confirmed, aes(date, cases)))
  gg.1
}

