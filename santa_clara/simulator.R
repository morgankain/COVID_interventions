### Simulate for all of the fits 

### My vision is to have 4 scenarios:
 ## 1) Lift shelter in place early
   ## -- FALSE, 60
 ## 2) Never lift shelter in palce
   ## -- FALSE, 120
 ## 3) Reduce shelter in place after 60 days
   ## -- TRUE, 60
 ## 4) Reduce shelter in place after 90 days
   ## -- TRUE, 90
scenario <- c("F60", "F120", "T60", "T90")

  ## Do we ever reduce shelter in place?
inf_iso      <- TRUE
 ## time shelter in place changes to a reduced form with infected isolation
red_shelt.t  <- 90
 ## strength of this new contact amount after red_shelt_t
red_shelt.s  <- 0.50
 ## nsims for each run
nsim         <- 100

## Load the data
SEIR.sim.ss.t.ci <- readRDS("output/SEIR.sim.ss.t.ci_ccc.Rds")
variable_params  <- readRDS("output/variable_params_ccc_H.Rds")

####
## Run ccc_forecast.R up to the variable params section
####

## Keep only the "best" fits
variable_params <- variable_params %>% 
  filter(
    log_lik > (max(log_lik) - 20)
  )

## Adjust variable params for the scenario
variable_params <- variable_params %>% 
  mutate(
    int_length2 = red_shelt.t
  ) %>% 
  mutate(
    iso_start   = int_start2 + red_shelt.t
  ) 

sim_start  <- variable_params$sim_start
sim_length <- 500
sim_end    <- sim_start + sim_length

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
  , scenario = scenario[2])

if (i == 1) {
  
SEIR.sim.f <- SEIR.sim  
  
} else {
  
SEIR.sim.f <- rbind(SEIR.sim.f, SEIR.sim) 
  
}

print(i)
  
}

SEIR.sim.f.s <- SEIR.sim.f %>% 
  dplyr::group_by(date) %>%
  dplyr::summarize(
    lwr = quantile(H, c(0.025))
  , est = quantile(H, c(0.50))
  , upr = quantile(H, c(0.975)))

# SEIR.sim.f.s.f240 <- SEIR.sim.f.s
# SEIR.sim.f.s.f120 <- SEIR.sim.f.s
# SEIR.sim.f.s.f90  <- SEIR.sim.f.s

# SEIR.sim.f.s.t120.35 <- SEIR.sim.f.s
# SEIR.sim.f.s.t90.35  <- SEIR.sim.f.s
# SEIR.sim.f.s.t60.35  <- SEIR.sim.f.s

# SEIR.sim.f.s.t120.50 <- SEIR.sim.f.s
# SEIR.sim.f.s.t90.50  <- SEIR.sim.f.s
# SEIR.sim.f.s.t60.50  <- SEIR.sim.f.s

# SEIR.sim.f.s.t120.65 <- SEIR.sim.f.s
# SEIR.sim.f.s.t90.65  <- SEIR.sim.f.s
# SEIR.sim.f.s.t60.65  <- SEIR.sim.f.s



