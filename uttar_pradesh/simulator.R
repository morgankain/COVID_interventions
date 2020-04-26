#########################################
### Simulate for Uttar Pradesh, India ###
#########################################

## Note: No fits here, just using parameter ranges and simulations

####
## Adjustable parameters
###################################

county.N      <- 232882078     ## pop size
lockdown.eff  <- c(0.05, 0.30) ## expected efficiency range of lockdown. Proportion of baseline contacts
                               ## Interventions. Only one at a time, TRUE + TRUE will break
epidemic_start <- c(40, 50)    ## Days since Jan 1 when we expect the epidemic to have started

inf_iso       <- TRUE          ## Reduce from shelter in place to some form of strong/moderate social distancing?
light         <- FALSE         ## Reduce from shelter in place to lightswitch
red_shelt.t   <- as.Date("2020/05/03") ## Day inf_iso or light begin
red_shelt.s   <- 0.65          ## strength of the new contact amount after red_shelt.t time elapsed since the start of shelter in place (proportion of baseline)
len.on        <- 21            ## length the intervention stays on for the lightswitch method
len.off       <- 21            ## length the intervention turns off during lightswitch
nsim          <- 100           ## number of epidemic simulations for each parameter set
sim_length    <- 500           ## how many days to run the simulation
state.plot    <- "D"           ## State variable for plotting (Hospit [H], Death [D], or Cases [C])
plot.log10    <- TRUE          ## log10 scale or not
nparams       <- 50            ## number of parameter combinations

###################################

####
## Setup stuff
####

needed_packages <- c("pomp", "plyr", "dplyr", "ggplot2", "magrittr", "scales", "lubridate", "tidyr", "data.table")
lapply(needed_packages, require, character.only = TRUE)

source("../ggplot_theme.R")

## Bring in pomp objects
source("pomp_objs.R")

## Load the current Contra Costa H for plotting purposes 
confirmed <- read.csv("confirmed.csv") %>% dplyr::select(UP, date)
confirmed <- confirmed %>% mutate(date = as.Date(as.character(date), format = "%m/%d/%Y"))
names(confirmed)[1] <- "cases"
dead      <- read.csv("dead.csv") %>% dplyr::select(UP, date)
dead      <- dead %>% mutate(date = as.Date(as.character(date), format = "%m/%d/%Y"))
names(dead)[1] <- "deaths"
dead      <- dead %>% mutate(total_deaths = cumsum(deaths))

## Load the current parameter set
params <- read.csv("params.csv", stringsAsFactors = FALSE)
params <- params %>% mutate(Value = sapply(Value, function(x) eval(parse(text = x))))

fixed_params        <- params$Value
names(fixed_params) <- params$Parameter

fixed_params        <- c(fixed_params, N = county.N)

## Load 1000 parameter combinations run for Santa Clara, extract the parameters within the top 1 likelihood unit
 ## and use those parameters to loop over 
variable_params.fit <- readRDS("../santa_clara/output/variable_params_scc_wide_params.Rds")

## Keep only the "best" fits. Here be very stringent 
variable_params.fit <- variable_params.fit %>% 
  filter(
    log_lik > (max(log_lik) - 1)
  )

 ## from here build a parameter set for UP, for now sample across what we think the shelter in place strength
  ## may be as well as what beta0 may be
variable_params <- sobolDesign(
  lower = c(
    E0          = 3
  , sim_start   = epidemic_start[1]
  , sd_m2       = lockdown.eff[1]
  , beta0est    = mean(variable_params.fit$beta0) - sd(variable_params.fit$beta0) 
  )
, upper = c(
    E0          = 3
  , sim_start   = epidemic_start[2]
  , sd_m2       = lockdown.eff[2]
  , beta0est    = mean(variable_params.fit$beta0) + sd(variable_params.fit$beta0) 
)
, nseq  = nparams
)

variable_params <- variable_params %>% mutate(
    sim_start  = as.Date(sim_start, origin = as.Date("2020-01-01"))
  , int_start2 = as.Date("2020-03-24")
  , paramset     = seq(1, nparams)
 )

variable_params <- variable_params %>% mutate(
    int_length2 = red_shelt.t - int_start2
 )

variable_params <- variable_params %>% 
  mutate(
    iso_start = int_start2 + int_length2
  , log_lik   = 0
  ) 

####
## Simulations
####

for (i in 1:nrow(variable_params)) {
  
county.data <- dead %>% mutate(day = round(as.numeric(date - variable_params[i, ]$sim_start)))
  
county.data <- rbind(
  data.frame(
    deaths = 0
  , date   = as.Date(seq(1:(min(county.data$day) - 1)), origin = variable_params[i, ]$sim_start)
  , day    = seq(1:(min(county.data$day) - 1))
  , total_deaths = 0
  )
, county.data
  )
  
## Create intervention covariate table for the full forecast
intervention.forecast <- with(variable_params[i, ], {

 covariate_table(
  day              = 1:sim_length
  
, intervention     = c(
      # No intervention until intervention start time
    rep(0, int_start2 - sim_start)                   
      # Intervention
  , rep(1, int_length2)
      # Post intervention close
  , {
    if (!inf_iso & !light) {
      rep(0, sim_length - (int_start2 - sim_start) - int_length2 + 1)
    } else if (inf_iso & !light) {
      rep(1, sim_length - (int_start2 - sim_start) - int_length2 + 1)  
    } else if (!inf_iso & light) {
      rep(2, sim_length - (int_start2 - sim_start) - int_length2 + 1)        
    }
    }
  ) 
  
, isolation = { 
  
if (!inf_iso) {
  
  rep(0, sim_length)
  
} else {
  
  c(
    rep(0, iso_start - sim_start)  
  , rep(1, sim_length - (iso_start - sim_start) + 1)
  )
  
}}

, iso_severe_level = rep(0, sim_length)      # % of contats that severe cases maintain
, iso_mild_level   = rep(0.05, sim_length)   # % of contats that mild cases maintain
  
, soc_dist_level = {
  
if (!inf_iso & !light) {
  c(                       # intensity of the social distancing interventions
  rep(1, int_start2 - sim_start)
     ## slightly odd way to do this, but should work
, rep(sd_m2, sim_length - (int_start2 - sim_start) + 1)
  )
  
} else if (inf_iso & !light) {
  
  c(                       # intensity of the social distancing interventions
  rep(1, int_start2 - sim_start)
     ## slightly odd way to do this, but should work
, rep(sd_m2, iso_start - int_start2)
     ## if infected isolation is going on, increase the background contact rate
, rep(red_shelt.s, sim_length - (iso_start - sim_start) + 1)
  )
  
} else if (!inf_iso & light) {
  
  c(                       # intensity of the social distancing interventions
  rep(1, int_start2 - sim_start)
     ## slightly odd way to do this, but should work
, rep(sd_m2, iso_start - int_start2)
     ## if infected isolation is going on, increase the background contact rate
, rep(
  c(rep(red_shelt.s, len.off), rep(sd_m2, len.on)), length = sim_length - (iso_start - sim_start) + 1)
  )  
  
}
  
  }
  
, thresh_H_start   = rep(NA, sim_length)  
, thresh_H_end     = rep(NA, sim_length)
   
# level of social distancing implemented with the threshhold intervention
, thresh_int_level = rep(sd_m2, sim_length)
, back_int_level   = rep(red_shelt.s, sim_length)
  
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

SEIR.sim <- SEIR.sim %>% mutate(date = round(as.Date(day, origin = variable_params[i, ]$sim_start)))

SEIR.sim <- SEIR.sim %>% mutate(paramset = variable_params[i, ]$paramset)

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

state.plot    <- "D"

SEIR.sim.f.m <- SEIR.sim.f %>% 
  dplyr::filter(.id == "median") %>% 
  dplyr::select(state.plot, date, paramset)

SEIR.sim.f.s <- SEIR.sim.f %>% 
  dplyr::filter(.id != "median") %>%
  dplyr::group_by(date) %>%
  dplyr::summarize(
    lwr = quantile(get(state.plot), c(0.050))
  , est = quantile(get(state.plot), c(0.50))
  , upr = quantile(get(state.plot), c(0.95)))

if (inf_iso & !light) {
  plot.title <- paste("Relax shelter in place to", red_shelt.s, "of baseline contacts on", red_shelt.t, sep = " ")
} else if (!inf_iso & light) {
  plot.title <- paste("Starting on", red_shelt.t, "every", len.on, "days, switch back and forth from social distancing of"
    , red_shelt.s, "and", round(mean(variable_params$sd_m2), 2), sep = " ")
} else {
  plot.title <- paste("Lift shelter in place on", red_shelt.t, sep = " ")
}

gg.1 <- ggplot(SEIR.sim.f.s) + 
  geom_line(aes(x = date, y = est), colour = "dodgerblue4", lwd = 1.5) + 
  geom_ribbon(aes(x = date, ymin = lwr, ymax = upr), colour = NA, fill = "dodgerblue4", alpha = 0.25) +
  geom_line(data = SEIR.sim.f.m, aes(x = date, y = get(state.plot), group = paramset), colour = "grey50", alpha = 0.35) +
  scale_x_date(labels = date_format("%Y-%b"), date_breaks = "1 month") +
  xlab("Date") + 
  ylab(state.plot) +
  guides(color = FALSE) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12)
  , legend.title = element_text(size = 12)
  , plot.title = element_text(size = 12)) +
  ggtitle(plot.title) 
#  ggtitle(
#"Starting on 2020-05-03 every 21 days, switch back 
#and forth from social distancing of 0.4 and 0.15"  
#  )

if (plot.log10) {
gg.1 <- gg.1 + scale_y_log10()
}

## No hospit data for UP
if (state.plot == "H") {
#  (gg.1 <- gg.1 + geom_point(data = hospit , aes(date, ch)))
  gg.1
} else if (state.plot == "D") {
  (gg.1 <- gg.1 + geom_point(data = dead, aes(date, total_deaths)))
} else if (state.plot == "C") {
  (gg.1 <- gg.1 + geom_point(data = confirmed, aes(date, cases)))
}

