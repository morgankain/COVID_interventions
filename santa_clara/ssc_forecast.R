############################################################################
## Use pomp to simulate dynamics up until time t and then project forward ##
############################################################################

# for Sherlock
# setwd("/scratch/users/kainm/covid")

## Other scripts used:
 ## scc_pomp_objs.R -- pomp stuff and R0 function
 ## scc_summary.R   -- summarize simulated dynamics
 ## scc_plotting.R  -- plot results

## Notes:
 ## parameters in params.R -- update there
 ## Fill in the "Parameters to adjust below" and then run from the top 

### Parameters to adjust for the given runs
focal.county  <- "Santa Clara"
county.N      <- 1.938e6
nparams       <- 200            ## number of parameter samples (more = longer)
nsim          <- 200            ## number of simulations for each fitted beta0

## Search !! for next steps
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
)

lapply(needed_packages, require, character.only = TRUE)

source("../ggplot_theme.R")

## Be very careful here, adjust according to your machine
 ## Not acutally used in the script right now, but important for expanding pomp fits
registerDoParallel(
  cores = 2
  )

## Bring in pomp objects
source("scc_pomp_objs.R")

## Add in infected isolation? (contact tracing?)
inf_iso <- TRUE

####
## Step 1: Pull the data
####

deaths     <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
deaths     <- deaths %>% mutate(date = as.Date(date)) %>% filter(county == focal.county)

####
## Step 2: Establish a reasonable parameter set
####

params <- read.csv("params.csv", stringsAsFactors = FALSE)
params <- params %>% mutate(Value = sapply(Value, function(x) eval(parse(text = x))))

fixed_params        <- params$Value
names(fixed_params) <- params$Parameter

fixed_params        <- c(fixed_params, N = county.N)

set.seed(10001)
 ## parameters that will vary
variable_params <- sobolDesign(
  lower = c(
    E0          = 1
  , sim_start   = 14
  , int_start1  = 60
  , int_length2 = 60
  , sd_m1       = 0.6
  , sd_m2       = 0.05
  )
, upper = c(
    E0          = 10
  , sim_start   = 28
  , int_start1  = 69
  , int_length2 = 120
  , sd_m1       = 0.9
  , sd_m2       = 0.3
)
, nseq  = nparams
) %>% mutate(
  sim_start   = round(sim_start)
, int_start1  = round(int_start1)
, int_length2 = round(int_length2)
, E0 = round(E0))

variable_params <- variable_params %>% mutate(
    sim_start  = as.Date(sim_start, origin = as.Date("2020-01-01"))
  , int_start1 = as.Date(round(int_start1), origin = as.Date("2020-01-01"))
  , int_start2 = as.Date("2020-03-17")
 ) %>% mutate(
    int_length1  = round(int_start2 - int_start1)
  ## Column for estimated beta0 (!! see lower down for desire to store likelyhood profile to define pmcmc runs)
  , beta0est     = 0
  , paramset     = seq(1, nparams)
) 

variable_params <- variable_params %>% 
  mutate(
    ## Previously had infected isolation starting after 60 days, changing this to when the SD gets lifted a bit
    iso_start = int_start2 + int_length2
  , log_lik   = 0
  ) 

####
## Step 3: Build the pomp object from data and simulation parameters
####

## Run parameters
sim_start  <- variable_params$sim_start
sim_length <- 500
sim_end    <- sim_start + sim_length

## container for results
SEIR.sim.ss.t.ci <- data.frame(
  name     = character(0)
, lwr      = numeric(0)
, est      = numeric(0)
, upr      = numeric(0)
, paramset = numeric(0)
)

for (i in 1:nrow(variable_params)) {
  
## Adjust the data for the current start date
county.data <- deaths %>% 
  mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) %>% 
  select(day, date, deaths) %>% 
  mutate(deaths_cum = deaths) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum)) %>% 
  replace_na(list(deaths = 0)) %>%
  select(-deaths_cum)
  
## Add days from the start of the sim to the first recorded day in the dataset
county.data <- rbind(
  data.frame(
    day    = seq(1:(min(county.data$day) - 1))
  , date   = as.Date(seq(1:(min(county.data$day) - 1)), origin = variable_params[i, ]$sim_start)
  , deaths = 0
  )
, county.data
  )

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
  , rep(1, int_length2)
      # Post intervention close
  , {
    if (!inf_iso) {
      rep(0, sim_length - (int_start2 - sim_start) - int_length2)
    } else {
      rep(1, sim_length - (int_start2 - sim_start) - int_length2)  
    }
    }
  ) 
  
, isolation = { 
  
if (!inf_iso) {
  
  rep(0, sim_length)
  
} else {
  
  c(
    rep(0, iso_start - sim_start)  
  , rep(1, sim_length - (iso_start - sim_start))
  )
  
}}

, iso_severe_level = rep(0, sim_length)      # % of contats that severe cases maintain
, iso_mild_level   = rep(0.05, sim_length)   # % of contats that mild cases maintain
  
, soc_dist_level = {
  
if (!inf_iso) {
   
  c(                       # intensity of the social distancing interventions
  rep(sd_m1, int_start2 - sim_start)
     ## slightly odd way to do this, but should work
, rep(sd_m2, sim_length - (int_start2 - sim_start)))
  
} else {
  
  c(                       # intensity of the social distancing interventions
  rep(sd_m1, int_start2 - sim_start)
     ## slightly odd way to do this, but should work
, rep(sd_m2, iso_start - int_start2)
     ## if infected isolation is going on, increase the background contact rate
, rep(sd_m1, sim_length - (iso_start - sim_start))
  )
  
}}
  
, thresh_H_start   = rep(2, sim_length)  
, thresh_H_end     = rep(10, sim_length)
   
     ## No reason to have a second parameter here, just use the same val that the user picks for social distancing
, thresh_int_level = c(                       # level of social distancing implemented with the threshhold intervention
  rep(sd_m1, int_start2 - sim_start)
     ## slightly odd way to do this, but should work
, rep(sd_m2, sim_length - (int_start2 - sim_start))
)
  
, order            = "constant"
, times            = "day"
  )

})

covid.fitting <- county.data %>%
  pomp(
    time       = "day"
  , t0         = 1
  , covar      = intervention.forecast
  , rprocess   = euler(sir_step, delta.t = 1/6)
  , rmeasure   = rmeas 
  , dmeasure   = dmeas
  , rinit      = sir_init
  , partrans   = par_trans
  , accumvars  = accum_names
  , paramnames = param_names
  , statenames = state_names
  ) 

####
## Step 3: Particle filter to get a fitted beta0
####

## So that the runs can be stopped and picked up again
if (variable_params[i, "beta0est"] == 0) {

timeoneparam <- system.time({
mifs_local <- covid.fitting %>%
  mif2(
    t0     = 1
  , params = c(fixed_params, c(beta0 = 2.5/7, E0 = variable_params[i, ]$E0))
  , Np     = 3000
  , Nmif   = 50
  , cooling.fraction.50 = 0.5
  , rw.sd  = rw.sd(beta0 = 0.02)
        )
})

# source("for_pmcmc.R")

}

variable_params[i, "beta0est"] <- coef(mifs_local)["beta0"]
variable_params[i, "log_lik"]  <- mifs_local@loglik

####
## Step 4: Simulate from the beginning and project forward with this beta0
## Would prefer pmcmc (see for_pmcmc) for uncertainty in beta0. Next step...
####

SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid.fitting
    , times        = intervention.forecast@times
    , params       = c(fixed_params, c(beta0 = variable_params[i, "beta0est"], E0 = variable_params[i, ]$E0))
    , nsim         = nsim
    , format       = "d"
    , include.data = F
    , seed         = 1001)) %>% {
      rbind(.,
         group_by(., day) %>%
           select(-.id) %>%
           summarise_all(median) %>%
                    mutate(.id = "median"))
      } 

## Plot of a single run
## source("scc_sr_plot.R")

####
## Step 5: Calculate summary statistics from these runs
####

## A) Rate of change calculated with hospitalizations
## B) Maximum number of hospitalizations at one point in time
## C) Timing of maximum number of hospitalizations

## Maybe somewhat controversial [?] choice here to remove all sims that don't take off as we
 ## know that this didn't happen in reality (keeping these would really skew our summary values)
SEIR.sim.s  <- SEIR.sim  %>% 
  dplyr::group_by(.id) %>% 
  dplyr::summarize(total_H = sum(H))

SEIR.sim    <- left_join(SEIR.sim, SEIR.sim.s, by = ".id")

SEIR.sim    <- SEIR.sim %>% 
  dplyr::filter(
    total_H > 10
  ) %>% droplevels()

## Calc the summary statistics
SEIR.sim.ss <- SEIR.sim %>% 
  dplyr::filter(.id != "median") %>%
  dplyr::mutate(week   = day %/% 7) %>%
  dplyr::group_by(week, .id) %>%
  dplyr::summarize(
    ## Mean in hospitalizations by week
    sum_H = sum(H_new)
    ) %>%
  dplyr::group_by(.id) %>%  
  dplyr::mutate(
    ## Difference in hospitalizations at the level of the week
    diff_H = c(0, diff(sum_H))
    ) %>%
  dplyr::group_by(.id) %>% 
  dplyr::summarize(
    ## Maximum hospitalizations reached, summarized at the level of the week
    when_max_H = week[which.max(sum_H)]
    ## How many hospitalizations are reached in that week
  , max_H      = max(sum_H)
    ## First week we see a reduction in the number of hospitalizations from a runs _global_ rate peak
  , when_red_H = week[min(which(diff_H == max(diff_H)))]
    )
  
SEIR.sim.ss2 <- SEIR.sim %>%
  filter(.id != "median") %>%
  dplyr::group_by(.id) %>% 
  dplyr::summarize(
    total_D = max(D)
  , total_R = max(R)
      )

SEIR.sim.ss.t    <- left_join(SEIR.sim.ss, SEIR.sim.ss2, by = ".id")

## Get their CI
SEIR.sim.ss.t.s <- SEIR.sim.ss.t %>%
  pivot_longer(cols = -.id) %>%
  dplyr::group_by(name) %>%
  dplyr::summarize(
    lwr = quantile(value, c(0.025))
  , est = quantile(value, c(0.500))
  , upr = quantile(value, c(0.975))
  ) %>% 
  dplyr::mutate(paramset = variable_params[i, ]$paramset)

SEIR.sim.ss.t.s <- rbind(SEIR.sim.ss.t.s
  , data.frame(
      name = "prop_ext"
    , lwr  = nsim - nrow(SEIR.sim.ss.t)
    , est  = nsim - nrow(SEIR.sim.ss.t)
    , upr  = nsim - nrow(SEIR.sim.ss.t)
    , paramset = variable_params[i, ]$paramset))

SEIR.sim.ss.t.ci <- rbind(SEIR.sim.ss.t.ci, SEIR.sim.ss.t.s)

print(i)

}

source("scc_summary.R")
source("scc_plotting.R")

