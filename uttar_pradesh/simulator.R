#########################################
### Simulate for Uttar Pradesh, India ###
#########################################

## Note: No fits here, just using parameter ranges and simulations

####
## Adjustable parameters
###################################

county.N      <- 232882078     ## pop size
lockdown.eff  <- c(0.10, 0.30) ## expected efficiency range of lockdown. Proportion of baseline contacts
                               ## Interventions. Only one at a time, TRUE + TRUE will break
epidemic_start <- c(45, 60)    ## Days since Jan 1 when we expect the epidemic to have started
beta0est       <- c(0.55, 0.9) ## Possible range of baseline transmission rate
inf_iso       <- FALSE          ## Reduce from shelter in place to some form of strong/moderate social distancing?
light         <- TRUE         ## Reduce from shelter in place to lightswitch
red_shelt.t   <- as.Date("2020/05/03") ## Day inf_iso or light begin
red_shelt.s   <- 0.40          ## strength of the new contact amount after red_shelt.t time elapsed since the start of shelter in place (proportion of baseline)
len.on        <- 21            ## length the intervention stays on for the lightswitch method
len.off       <- 14            ## length the intervention turns off during lightswitch
nsim          <- 100           ## number of epidemic simulations for each parameter set
sim_length    <- 500           ## how many days to run the simulation
state.plot    <- "D"           ## State variable for plotting (Hospit [H], Death [D], or Cases [C])
plot.log10    <- TRUE          ## log10 scale or not
nparams       <- 20            ## number of parameter combinations

###################################
### Assumed parameter ranges. Adjustable by opening params.csv
###################################

## parameter   est  lwr upr       details
##
## Ca	         0.6	0.4	0.8	      Relative infectiousness of asymptomatic infected individuals
## Cp	         1			            Relative infectiousness of preymptomatic infected individuals
## Cs	         1			            Relative infectiousness of severely symptomatic infected individuals
## Cm	         1			            Relative infectiousness of mildly symptomatic infected individuals
## alpha  	   0.43	0.3	0.5	      Proportion of infections that are asymptomatic
## gamma	     1/3.5			        Incubation period
## lambda_a 	 1/7			          Infectious period for asymptomatic infections
## lambda_s	   1/5.5			        Infectious period for severely symptomatic infections
## lambda_m	   1/5.5		        	Infectious period for mildly symptomatic infections
## lambda_p	   1/1.5		        	Infectious period for presymptomatic infections
## delta	     0.2	0.1	0.3     	Fatality rate among hospitalizations
## mu	         0.950	0.925	0.975	Proportion of symptomatic infections that require hospitalization
## rho_d	     1/13.3		        	1 / duration of hospital stay of those that die of infection
## rho_r	     1/15		          	1 / duration of hospital stay of those that recover from infection
## import_rate 1/4		           	Importation rate of infections

####
## Setup stuff
####

## collapsable
{
needed_packages <- c("pomp", "plyr", "dplyr", "ggplot2", "magrittr", "scales", "lubridate", "tidyr", "data.table", "gridExtra")
lapply(needed_packages, require, character.only = TRUE)

source("../ggplot_theme.R")

## Bring in pomp objects
source("COVID_pomp.R")

## Load the current Contra Costa H for plotting purposes 
confirmed <- read.csv("confirmed.csv") %>% dplyr::select(UP, Date)
confirmed <- confirmed %>% mutate(Date = as.Date(as.character(Date), format = "%m/%d/%Y"))
names(confirmed) <- c("cases", "date")
dead      <- read.csv("dead.csv") %>% dplyr::select(UP, Date)
dead      <- dead %>% mutate(Date = as.Date(as.character(Date), format = "%m/%d/%Y"))
names(dead) <- c("deaths", "date")
dead      <- dead %>% mutate(total_deaths = cumsum(deaths))

## Load the current parameter set
params <- read.csv("params.csv", stringsAsFactors = FALSE)
params <- params %>% mutate(Value = sapply(est, function(x) eval(parse(text = x))))

fixed_params        <- params$Value
names(fixed_params) <- params$Parameter

fixed_params        <- c(fixed_params, N = county.N)

## Parameters to vary across simulations
source("variable_params_less.R")

variable_params <- variable_params %>% 
  mutate(
    int_length2        = (red_shelt.t - int_start2)
  , iso_start          = int_start2 + (red_shelt.t - int_start2)
  , soc_dist_level_red = red_shelt.s
  ) 
}

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
   
, light_on = c( 
 ## On until start of intervention
  rep(1, iso_start - sim_start)
 ## And then flickers according to the params
, rep(
  c(rep(0, len.off), rep(1, len.on)), length = sim_length - (iso_start - sim_start))
  )
   
# level of social distancing implemented with the threshhold intervention
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

, iso_severe_level = rep(0, sim_length)      # % of contats that severe cases maintain
, iso_mild_level   = rep(0.05, sim_length)   # % of contats that mild cases maintain
  
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
   
, soc_dist_level_wfh = rep(soc_dist_level_wfh, sim_length) 
   
, thresh_H_start   = rep(NA, sim_length)  
, thresh_H_end     = rep(NA, sim_length)
   
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
  , rmeasure   = rmeas_deaths 
  , dmeasure   = dmeas_deaths
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
    , params       = c(fixed_params
      , c(
      beta0              = variable_params[i, "beta0est"]
    , soc_dist_level_sip = variable_params[i, "soc_dist_level_sip"]
    , E0                 = variable_params[i, ]$E0
    , Ca                 = variable_params[i, ]$Ca
    , alpha              = variable_params[i, ]$alpha
    , delta              = variable_params[i, ]$delta
    , mu                 = variable_params[i, ]$mu
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
  plot.title <- paste("Starting on", red_shelt.t, len.on, "days of strong SD", len.off, "days of softer SD =", red_shelt.s, sep = " ")
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
  gg.2 <- gg.1 + geom_point(data = confirmed, aes(date, cases * 10)) + 
    ggtitle("Rough rule of thumb that we catch ~ 1/15 cases")
  gg.1 <- gg.1 + geom_point(data = confirmed, aes(date, cases))
 
  grid.arrange(gg.1, gg.2, ncol = 1)
}

