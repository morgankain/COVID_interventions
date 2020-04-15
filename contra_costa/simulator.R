##############################
### Simulate from the fits ###
##############################

needed_packages <- c("pomp", "plyr", "dplyr", "ggplot2", "magrittr", "scales", "lubridate", "tidyr", "data.table")
lapply(needed_packages, require, character.only = TRUE)

####
## Parameters
####

## Note:
  ## There are 373 stored parameter sets with a fitted beta0 with loglikelihood within 3 units of the max
  ## Using all 373 will be quite slow, so can specify if you want a random subset (set nparams < 373)
nparams       <- 373

focal.county  <- "Contra Costa"
county.N      <- 1.147e6

  ## Do we ever reduce from shelter in place to some form of strong/moderate social distancing?
inf_iso       <- FALSE
 ## time shelter in place changes to a reduced form with infected isolation
  ## ignored if inf_iso = FALSE
red_shelt.t   <- 120
 ## strength of the new contact amount after red_shelt.t time elapsed since the start of shelter in place (proportion of baseline)
red_shelt.s   <- 0.50
 ## number of epidemic simulations for each parameter set
nsim          <- 100
 ## how many days to run the simulation
sim_length    <- 500
 ## State variable for plotting (Hospit, Death, or Cases)
   ## if H, D, or C plot H, D or C from the data on the plot
state.plot    <- "C"
 ## log10 scale or not
plot.log10    <- TRUE

## Load the current Contra Costa H for plotting purposes 
 ## !! Update with new data. Data only through April 10 here currently
hospit     <- read.csv("ccc_data.csv")
hospit     <- hospit %>% 
  mutate(date = as.Date(REPORT_DATE)) %>% 
  filter(CURRENT_HOSPITALIZED != "NULL") %>% 
  mutate(ch = as.numeric(as.character(CURRENT_HOSPITALIZED))) %>% 
  dplyr::select(date, ch)

### Collapsable for convenience
{
####
## Setup stuff
####

source("../ggplot_theme.R")

## Bring in pomp objects
source("../santa_clara/scc_pomp_objs.R")

## Garb new death data from online or load the data in the github
#deaths   <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
deaths    <- read.csv("us-counties.txt")
deaths    <- deaths %>% mutate(date = as.Date(date)) %>% filter(county == focal.county)

## Load the current parameter set
params <- read.csv("params.csv", stringsAsFactors = FALSE)
params <- params %>% mutate(Value = sapply(Value, function(x) eval(parse(text = x))))

fixed_params        <- params$Value
names(fixed_params) <- params$Parameter

fixed_params        <- c(fixed_params, N = county.N)

####
## Previous fits
####

## Load the previous fits
variable_params  <- read.csv("variable_params.csv"); variable_params <- variable_params %>% dplyr::select(-X)

## notes on this data frame:
 ## E0 = starting infecteds (unlikely to be higher than 1-4, but stochasticity is too high if 1). 
   ## Higher number will require a bit of a later start, but should be fine
 ## sim_start = first case
 ## int_stat1 = first suggested "work from home" in the bay area
 ## int_length2 = length of shelter in place. Overwritten by your parameter choices shortly
 ## sd_m1     = social distancing strength for the "work from home" suggestion
 ## sd_m2     = social distancing strength for shelter in place
 ## int_start2 = start date of shelter in place
 ## int_length1 = length of the "work from home" prior to shelter in place
 ## beta0est    = estimated beta0
 ## iso_start   = start of infected isolation. Overwritten shortly. Ignored if inf_iso = FALSE

## Keep only the "best" fits
variable_params <- variable_params %>% 
  filter(
    log_lik > (max(log_lik) - 3)
  )

## Adjust variable params for the scenario
variable_params <- variable_params %>% 
  mutate(
    int_length2 = red_shelt.t
  ) %>% 
  mutate(
    iso_start   = int_start2 + red_shelt.t
  ) 

if (nparams < 373) {
  variable_params <- variable_params[sample(1:373, nparams), ]
}

####
## Simulations
####

for (i in 1:nrow(variable_params)) {
  
## Adjust death data by the start date for this specific parameter set. Not used for fitting here in any way
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
, rep(red_shelt.s, sim_length - (iso_start - sim_start))
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

SEIR.sim <- SEIR.sim %>% filter(.id != "median")
SEIR.sim <- SEIR.sim %>% mutate(date = as.Date(day, origin = variable_params[i, ]$sim_start))

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
  }

####
## Summary and plotting
####

## summary of observable cases to plot against case data in the county
SEIR.sim.f <- SEIR.sim.f %>% mutate(C = Is + Im + H)

SEIR.sim.f.s <- SEIR.sim.f %>% 
  dplyr::group_by(date) %>%
  dplyr::summarize(
    lwr = quantile(get(state.plot), c(0.025))
  , est = quantile(get(state.plot), c(0.50))
  , upr = quantile(get(state.plot), c(0.975)))

data.frame(
  R0  = c("base", "under shelter in place")
, lwr = c(quantile(variable_params$R0_base, c(0.025)), quantile(variable_params$R0_int, c(0.025)))
, est = c(quantile(variable_params$R0_base, c(0.500)), quantile(variable_params$R0_int, c(0.500)))
, upr = c(quantile(variable_params$R0_base, c(0.975)), quantile(variable_params$R0_int, c(0.975)))
)

if (inf_iso) {
  plot.title <- paste("Relax shelter in place to", red_shelt.s, "of baseline contacts on day", red_shelt.t, sep = " ")
} else {
  plot.title <- paste("Lift shelter in place on day", red_shelt.t, sep = " ")
}

gg.1 <- ggplot(SEIR.sim.f.s) + 
  geom_line(aes(x = date, y = est), colour = "dodgerblue4") + 
  geom_ribbon(aes(x = date, ymin = lwr, ymax = upr), colour = NA, fill = "dodgerblue4", alpha = 0.4) +
  scale_x_date(labels = date_format("%Y-%b"), date_breaks = "1 month") +
  xlab("Date") + 
  ylab(state.plot) +
  guides(color = FALSE) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)
  , legend.title = element_text(size = 12)) +
  ggtitle(plot.title) 

if (plot.log10) {
gg.1 <- gg.1 + scale_y_log10()
}

if (state.plot == "H") {
  (gg.1 <- gg.1 + geom_point(data = hospit , aes(date, ch)))
} else if (state.plot == "D") {
  (gg.1 <- gg.1 + geom_point(data = deaths, aes(date, deaths)))
} else if (state.plot == "C") {
  (gg.1 <- gg.1 + geom_point(data = deaths, aes(date, cases)))
}
