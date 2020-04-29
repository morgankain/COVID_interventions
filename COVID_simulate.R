##################################################################
## Use pomp fits to simulate dynamics for various interventions ##
##################################################################

set.seed(10001)
fitting        <- FALSE   ## Small change in pomp objects if fitting or simulating
use.rds        <- TRUE
rds.name       <- "output/Contra Costa_TRUE_FALSE_0_2020-04-27_temp.Rds"
nsim           <- 300     ## Number of epidemic simulations for each parameter set
#### ONLY ONE AT A TIME ALLOWED RIGHT NOW
inf_iso        <- FALSE   ## Do we ever reduce from shelter in place to some form of strong/moderate social distancing?
light          <- TRUE    ## Lightswitch method
red_shelt.t    <- 76      ## time shelter in place changes to a reduced form (inf_iso or light) -- Now June 1
red_shelt.s    <- 0.5     ## new social dist strength after time red_shelt.t
thresh_H.start <- 15      ## Threshold when lightswtich turns on (when we get higher than this)
thresh_H.end   <- 5       ## Threshold when lightswtich turns off (when we drop from above to this value)
sim_length     <- 500     ## how many days to run the simulation
state.plot     <- "H"     ## State variable for plotting (Hospit [H], Death [D], or Cases [C])
focal.county   <- "Contra Costa"
#county.N       <- 1.938e6
county.N      <- 1.147e6 ## County population size
loglik.thresh  <- 2       ## Keep runs within top X loglik units
params.all     <- TRUE    ## Keep all fitted parameters above loglik thresh?...
nparams        <- 100     ## ...if FALSE, pick a random subset for speed
plot.log10     <- FALSE   ## log10 scale or not
fit.with       <- "H"
fit.minus      <- 0       ## Use data until X days prior to the present

test_and_isolate_s <- 0.1 ## Additional proportional reduction of severe cases under test and isolate
test_and_isolate_m <- 0.15 ## Additional proportional reduction of mild cases under test and isolate
  
needed_packages <- c(
    "pomp"
  , "plyr"
  , "dplyr"
  , "ggplot2"
  , "magrittr"
  , "scales"
  , "lubridate"
  , "tidyr"
  , "data.table")

lapply(needed_packages, require, character.only = TRUE)

source("ggplot_theme.R")
source("COVID_pomp.R")

#deaths <- read.csv("us-counties.txt")
#deaths <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
#deaths <- deaths %>% mutate(date = as.Date(date)) %>% filter(county == focal.county)

if (fit.with == "D") {
# deaths   <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
deaths     <- read.csv("us-counties.txt")
deaths     <- deaths %>% mutate(date = as.Date(date)) %>% filter(county == focal.county)
deaths     <- deaths %>% dplyr::filter(date < max(date) - fit.minus)
} else if (fit.with == "H") {
## !! Not supported right now for SCC, but placing here for completeness
   ## !! Right now only usable for CCC
hospit     <- read.csv("contra_costa/ccc_data.csv")
hospit     <- hospit %>% 
  mutate(date = as.Date(REPORT_DATE)) %>% 
  filter(CURRENT_HOSPITALIZED != "NULL") %>% 
  mutate(ch = as.numeric(as.character(CURRENT_HOSPITALIZED))) %>% 
  dplyr::select(date, ch)
hospit    <- hospit %>% dplyr::filter(date < max(date) - fit.minus)  
}

## Load the previously saved fits
if (use.rds) {
prev.fit         <- readRDS(rds.name)
variable_params  <- prev.fit[["variable_params"]]
fixed_params     <- prev.fit[["fixed_params"]]
} else {
## use fits that were just run with COVID_fit.R. and do nothing here

#params <- read.csv("params.csv", stringsAsFactors = FALSE)
#params <- params %>% mutate(Value = sapply(est, function(x) eval(parse(text = x))))

#fixed_params        <- params$Value
#names(fixed_params) <- params$Parameter

#fixed_params        <- c(fixed_params, N = county.N) 
  
}

## drop the rows that have 0s for likelihood (didnt' run) and keep only the "best" fits
variable_params <- variable_params %>% 
  filter(log_lik != 0) %>% 
# filter(log_lik == max(log_lik))
  filter(log_lik > (max(log_lik) - loglik.thresh))

## for preprint plot counterfactual
# variable_params <- variable_params %>% mutate(int_start2 = int_start2 + 7) 

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

## debug
# variable_params$soc_dist_level_sip <- 0.2

for (i in 1:nrow(variable_params)) {
  
if (fit.with == "D") {
  
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

####
## Step 3: Simulate from the beginning and project forward with this beta0
## Would prefer pmcmc (see for_pmcmc) for uncertainty in beta0. Next step...
####

SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid.fitting
    , times        = intervention.forecast@times
    , params = c(fixed_params
    , c(
      beta0              = variable_params[i, "beta0est"]
    , soc_dist_level_sip = variable_params[i, "soc_dist_level_sip"]
    , E0                 = variable_params[i, ]$E0
    , Ca                 = variable_params[i, ]$Ca
    , alpha              = variable_params[i, ]$alpha
    , delta              = variable_params[i, ]$delta
    , mu                 = variable_params[i, ]$mu
   , rho_d              = variable_params[i, ]$rho_d
   , rho_r              = variable_params[i, ]$rho_r
      ))
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
## source("single_run_plot.R")

SEIR.sim <- SEIR.sim %>% mutate(date = round(as.Date(day, origin = variable_params[i, ]$sim_start)))

SEIR.sim <- SEIR.sim %>% 
  mutate(
    paramset = variable_params[i, ]$paramset
  )

if (i == 1) {
  
SEIR.sim.f <- SEIR.sim  
  
} else {
  
SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim) 
  
}

if (((i / 20) %% 1) == 0) {
  print(paste(round(i / nrow(variable_params), 2)*100, "% Complete", sep = ""))
}

}

## because of the staggered start dates there is some oddity at the end of the simulation
SEIR.sim.f <- SEIR.sim.f %>% filter(date < min(variable_params$sim_start + sim_length))

## summary of observable cases to plot against case data in the county
SEIR.sim.f <- SEIR.sim.f %>% mutate(C = Is + Im + H)

####
## Summary and plotting
####

state.plot    <- "H"

SEIR.sim.f.m <- SEIR.sim.f %>% 
  dplyr::filter(.id == "median") %>% 
  dplyr::select(state.plot, date, paramset)

SEIR.sim.f.s <- SEIR.sim.f %>% 
  dplyr::filter(.id != "median") %>%
  dplyr::group_by(date) %>%
  dplyr::summarize(
#    lwr = quantile(get(state.plot), c(0.025))
#  , est = quantile(get(state.plot), c(0.50))
#  , upr = quantile(get(state.plot), c(0.975))
    lwr = quantile(get(state.plot), c(0.05))
  , est = quantile(get(state.plot), c(0.50))
  , upr = quantile(get(state.plot), c(0.95))
  )

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
# (gg.1 <- gg.1 + geom_point(data = confirmed, aes(date, cases)))
  gg.1
}
